//! Core module for detecting open reading frames in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for finding open reading frames (ORFs)
//! in a set of aligned reads.
//!
//! In short, every possible open reading frame (ORF) is detected for every
//! read in the query set. For every potential ORF, learning models and databases
//! are used to determine whether the ORF is a true ORF, a false positive.
//! All the data from each reliable ORF is collected and subjected to another
//! learning model trained with true ORFs and false positives. The process is
//! heavily parallelized to offer fast performance on large datasets.

use dashmap::DashMap;
use genepred::{Bed12, GenePred};
use hashbrown::HashMap;
use memchr::memchr;
use memmap2::Mmap;
use rayon::prelude::*;
use smol_str::{SmolStr, ToSmolStr};

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::from_utf8;
use std::{fmt::Debug, io::Read};

use crate::consts::*;

/// Reads the entire content of a file into a `String`.
///
/// This function provides a basic utility for synchronously reading a file's
/// contents. It's generic over any type that can be converted to a `Path` and
/// is `Debug` printable.
///
/// # Arguments
///
/// * `file` - The path to the file to read.
///
/// # Returns
///
/// A `Result<String, Box<dyn std::error::Error>>` containing the file's
/// contents on success, or an error if the file cannot be opened or read.
///
pub fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

/// Runs a command through `bash -c` and exits with the tool's status unless it succeeds.
///
/// A failed tool must not look like a successful run without results: a DIAMOND
/// run that exits non-zero used to leave an empty hit table, which the blast step
/// then wrote out as "no protein evidence" for every candidate of the chunk. The
/// process exits with the tool's own code rather than panicking (101), so the task
/// fails the way the tool did: 137 after an out-of-memory kill is one of the codes
/// xorf's errorStrategy retries with more memory.
///
/// # Arguments
///
/// * `cmd` - The command line.
/// * `tool` - Name used in the error message.
///
/// # Example
///
/// ```rust, ignore
/// run_command("diamond version", "diamond");
/// ```
pub fn run_command(cmd: &str, tool: &str) {
    let status = std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute {tool} command -> {e}"));

    if let Err(code) = exit_code(&status) {
        eprintln!("ERROR: {tool} exited with {status} -> {cmd}");
        std::process::exit(code);
    }
}

/// `Ok` for a successful run, otherwise the code to exit with: the command's own
/// exit code, or 128 + signal when the shell itself was killed.
fn exit_code(status: &std::process::ExitStatus) -> Result<(), i32> {
    use std::os::unix::process::ExitStatusExt;

    if status.success() {
        return Ok(());
    }
    Err(status
        .code()
        .unwrap_or_else(|| 128 + status.signal().unwrap_or(0)))
}

/// Translates a sequence into amino acids.
///
/// # Arguments
///
/// * `sequence` - The sequence to translate.
///
/// # Returns
///
/// The translated sequence as a string.
///
/// # Panics
///
/// This function will panic if:
/// - A codon is not a valid codon.
///
/// # Example
///
/// ```rust, ignore
/// let sequence = vec![b'A', b'T', b'G', b'C'];
///
/// let aa = translate(&sequence);
/// ```
pub fn translate(sequence: &[u8]) -> String {
    let mut aa = String::new();
    for codon in sequence.chunks(3) {
        // INFO: codon to upper case
        let codon = codon.to_ascii_uppercase();
        let amino_acid = translate_codon(&codon);

        if amino_acid == "X" {
            println!(
                "WARN: codon -> {:?} is not a valid codon from sequence -> {:?}!",
                from_utf8(&codon).unwrap(),
                from_utf8(sequence).unwrap()
            );

            aa.push('X'); // WARN: accepting invalid codons!
        }
        aa.push_str(amino_acid);
    }
    aa
}

/// Translates a codon into an amino acid.
///
/// # Arguments
///
/// * `codon` - The codon to translate.
///
/// # Returns
///
/// The amino acid corresponding to the given codon.
///
/// # Panics
///
/// This function will panic if:
/// - The codon is not a valid codon.
///
/// # Example
///
/// ```rust, ignore
/// let codon = vec![b'A', b'T', b'G'];
///
/// let amino_acid = translate_codon(codon);
/// ```
pub fn translate_codon(codon: &[u8]) -> &'static str {
    for (table_codon, amino_acid) in &CODON_TABLE {
        if codon == *table_codon {
            return amino_acid;
        }
    }

    "X" // INFO: unknown codon
}

/// Scans the sequence for stop codons.
///
/// # Arguments
///
/// * `sequence` - The sequence to scan.
///
/// # Returns
///
/// The number of stop codons found in the sequence.
///
/// # Example
///
/// ```rust, ignore
/// let sequence = vec![b'A', b'T', b'G', b'C'];
///
/// let stops = scan_stops(sequence);
/// ```
pub fn scan_stops(sequence: Vec<u8>) -> usize {
    sequence
        .windows(3)
        .enumerate()
        .filter(|(idx, window)| idx % 3 == 0 && STOP_CODONS_BYTES.contains(window))
        .count()
}

/// Checks if the start codon is valid.
///
/// # Arguments
///
/// * `start_codon` - The start codon to check.
/// * `id` - The ID of the gene prediction.
/// * `orf_start` - The start position of the ORF.
///
/// # Panics
///
/// This function will panic if:
/// - The start codon is not long enough.
/// - The start codon is not ATG.
///
/// # Example
///
/// ```rust, ignore
/// let start_codon = "ATG";
/// let id = "ENSMUST00000176751.2";
/// let orf_start = 128352418;
///
/// __check_start_codon(start_codon, id, orf_start);
/// ```
pub fn __check_start_codon(start_codon: &str, id: &str, orf_start: u64) {
    if start_codon.len() < 3 {
        panic!(
            "ERROR: start codon is not long enough -> {start_codon:?} from {id:?} using {orf_start:?}"
        );
    }

    if start_codon != START_CODON {
        dbg!(
            "WARN: start codon is not ATG -> {:?} from {:?} using {:?}",
            start_codon,
            &id,
            orf_start
        );
    }
}

/// Reformats a FASTA file.
///
/// # Arguments
///
/// * `fasta` - The path to the FASTA file.
/// * `records` - A hash map of gene predictions.
/// * `outdir` - The output directory.
///
/// # Returns
///
/// A tuple containing the path to the formatted FASTA file and a hash map of sequences.
///
/// # Panics
///
/// This function will panic if:
/// - The FASTA file cannot be parsed.
/// - The formatted FASTA file cannot be created.
///
/// # Example
///
/// ```rust, ignore
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::BufWriter;
///
/// let fasta = PathBuf::from("path/to/fasta.fa");
/// let records = HashMap::new();
/// let outdir = PathBuf::from("path/to/outdir");
///
/// let (fmt, seqs) = refmt(&fasta, &records, &outdir);
/// ```
pub fn refmt(
    fasta: &PathBuf,
    records: &HashMap<String, GenePred>,
    outdir: &Path,
) -> (PathBuf, HashMap<SmolStr, Vec<u8>>) {
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let fmt = outdir.join(format!(
        "{}.fmt.fa",
        fasta.file_stem().unwrap().to_str().unwrap()
    ));

    let mut writer = create_fasta(&fmt, "fa")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", fmt));

    let accumulator = DashMap::new();

    records.into_par_iter().for_each(|(header, record)| {
        let start = record.start;
        let end = record.end;
        let name = record
            .name()
            .unwrap_or_else(|| panic!("ERROR: name not found for {header:?}!"));
        let chrom = record.chrom();

        let fmt_header = format!(
            ">{}:{}-{}({})({})(0, 0,)",
            from_utf8(chrom).unwrap_or_else(|e| panic!("ERROR: failed to parse chrom -> {e}")),
            start,
            end,
            record.strand.unwrap(),
            from_utf8(name).unwrap_or_else(|e| panic!("ERROR: failed to parse name -> {e}"))
        );

        let seq = seqs.get(&header.clone().to_smolstr()).unwrap_or_else(|| {
            panic!("ERROR: {header} not found in sequences {:?}!", seqs);
        });

        accumulator
            .entry(seq)
            .or_insert_with(Vec::new)
            .push(fmt_header)
    });

    // INFO: will write the first header on collection and will point other headers to it
    accumulator.into_iter().for_each(|(seq, headers)| {
        for header in headers {
            writeln!(writer, "{}", &header).unwrap();
            writer.write_all(seq).unwrap();
            writeln!(writer).unwrap();
        }
    });

    (fmt.clone(), seqs)
}

/// Parses a FASTA file into a hash map.
///
/// # Arguments
///
/// * `f` - The path to the FASTA file.
///
/// # Returns
///
/// A `HashMap<SmolStr, Vec<u8>>` instance containing the parsed FASTA data.
///
/// # Panics
///
/// This function will panic if:
/// - The FASTA file cannot be read.
///
/// # Example
///
/// ```rust, ignore
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::BufWriter;
///
/// let fasta = PathBuf::from("path/to/fasta.fa");
///
/// let seqs = parse_fa(&fasta);
/// ```
pub fn parse_fa<F: AsRef<Path>>(
    f: F,
) -> Result<HashMap<SmolStr, Vec<u8>>, Box<dyn std::error::Error>> {
    let file = File::open(f).unwrap();
    let mmap = unsafe { Mmap::map(&file).unwrap() };
    let data = mmap.as_ref();

    let mut acc = HashMap::new();
    let mut pos = 0;

    while let Some(start) = memchr(b'>', &data[pos..]) {
        let start = pos + start;
        let end = memchr(b'>', &data[start + 1..]).map_or(data.len(), |e| start + 1 + e);
        let entry = &data[start + 1..end];
        let header_end = memchr(b'\n', entry).unwrap();
        let header = from_utf8(&entry[..header_end])
            .unwrap()
            .trim()
            .split(' ')
            .next()
            .unwrap();
        let record = &entry[header_end + 1..];
        let seq = record
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r') // Remove newdatas and carriage returns
            .cloned()
            .collect::<Vec<u8>>()
            .to_ascii_uppercase(); // Always uppercase to avoid softmasking codon matching issues

        acc.insert(SmolStr::new(header), seq);
        pos = end;
    }

    Ok(acc)
}

/// Creates a FASTA file with the given extension.
///
/// # Arguments
///
/// * `fasta` - The path to the FASTA file.
/// * `extension` - The extension to use for the FASTA file.
///
/// # Returns
///
/// A `BufWriter<File>` instance if the file was created successfully, otherwise `None`.
///
/// # Panics
///
/// This function will panic if:
/// - The file already exists.
/// - The file cannot be created.
///
/// # Example
///
/// ```rust, ignore
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::BufWriter;
///
/// let fasta = PathBuf::from("path/to/fasta.fa");
/// let extension = "fmt.fa";
///
/// let writer = create_fasta(&fasta, extension).unwrap();
/// ```
pub fn create_fasta(fasta: &Path, extension: &str) -> Option<BufWriter<File>> {
    let file = fasta.with_extension(extension);
    if !file.exists() {
        let writer = BufWriter::new(
            File::create(&file).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
        );

        Some(writer)
    } else {
        // INFO: overwrite existing file
        std::fs::remove_file(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot remove existing file -> {e}");
        });

        Some(BufWriter::new(File::create(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot create file -> {e}")
        })))
    }
}

/// Read a bed file and convert it to a hashmap of genepred records
///
/// # Arguments
///
/// * `bed` - The path to the bed file
///
/// # Returns
///
/// * `HashMap<String, GenePred>` - A hashmap of genepred records
///
/// # Example
///
/// ```rust
/// let bed = get_bed(&args.bed);
/// ```
pub fn get_bed(bed: &PathBuf) -> HashMap<String, GenePred> {
    genepred::Reader::<Bed12>::from_mmap(bed)
        .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
        .filter_map(|record| {
            record.ok().map(|record| {
                (
                    from_utf8(record.name().unwrap()).unwrap().to_string(),
                    record,
                )
            })
        })
        .collect::<HashMap<String, genepred::GenePred>>()
}

/// Given a position within the concatenated exonic sequence, this function
/// maps it back to its corresponding genomic coordinate.
///
/// The function iterates through the gene's exons, calculating the length of
/// each exon. It keeps a running total of the length of the exons processed so
/// far (`current_pos`). If the target position (`pos`) falls within the current
/// exon, it calculates the offset from the exon's start and returns the absolute
/// genomic coordinate. If the target position is beyond all exons, it returns `None`.
///
/// # Arguments
///
/// * `pos` - A `u64` representing the position within the concatenated exonic sequence (0-indexed).
///
/// # Returns
///
/// An `Option<u64>` which is:
/// * `Some(genomic_position)` if the position is found within the exons.
/// * `None` if the position is outside the total length of all exons.
///
/// # Panics
///
/// This function does not panic under normal conditions.
///
#[derive(Debug, Clone, Copy)]
enum TranscriptPosition {
    Start,
    End,
}

fn get_pos_in_exons(record: &GenePred, pos: u64) -> Option<u64> {
    get_pos_in_exons_by_kind(record, pos, TranscriptPosition::Start)
}

fn get_end_pos_in_exons(record: &GenePred, pos: u64) -> Option<u64> {
    get_pos_in_exons_by_kind(record, pos, TranscriptPosition::End)
}

fn get_pos_in_exons_by_kind(record: &GenePred, pos: u64, kind: TranscriptPosition) -> Option<u64> {
    let mut current_pos = 0;
    let mut exons = record.exons();

    exons = match record.strand() {
        Some(genepred::Strand::Forward) => exons,
        Some(genepred::Strand::Reverse) => {
            // INFO: we scale the coordinates and then sort them
            // INFO: this is equivalent to reversing the exons and reversing the coordinates
            let mut exons = exons
                .iter()
                .map(|(start, end)| (SCALE - end, SCALE - start))
                .collect::<Vec<(u64, u64)>>();

            exons.sort_by(|a, b| a.0.cmp(&b.0));
            exons
        }
        Some(genepred::Strand::Unknown) | None => {
            panic!("ERROR: unexpected strand value: {:?}", record.strand())
        }
    };

    for (exon_start, exon_end) in exons.iter() {
        let block_len = exon_end - exon_start; // INFO: exon length
        let exon_transcript_end = current_pos + block_len;

        let is_in_exon = match kind {
            TranscriptPosition::Start => pos < exon_transcript_end,
            TranscriptPosition::End => pos <= exon_transcript_end,
        };

        if is_in_exon {
            let offset = pos - current_pos;
            return Some(exon_start + offset);
        }

        current_pos += block_len;
    }

    // INFO: edge case -> dispute between 0-based and 1-based coords [always ends]
    if pos >= current_pos {
        return Some(exons.last().unwrap().1);
    }

    None
}

/// Maps the start and end positions of an Open Reading Frame (ORF) from the
/// concatenated exonic sequence to absolute genomic coordinates.
///
/// This function is crucial for converting ORF coordinates, which are relative to the
/// combined sequence of all exons, into the real genomic coordinates (chromosome, start,
/// end). It calls `get_pos_in_exons` for both the ORF's start and end positions.
///
/// The function also accounts for the gene's strand. For the forward strand, the
/// converted coordinates are returned directly. For the reverse strand, the coordinates
/// are flipped and adjusted based on a `SCALE` value to correctly reflect their position
/// on the reverse complement.
///
/// # Arguments
///
/// * `orf_start` - A `u64` representing the starting position of the ORF within the
///   concatenated exonic sequence.
/// * `orf_end` - A `u64` representing the ending position of the ORF within the
///   concatenated exonic sequence.
///
/// # Returns
///
/// A tuple `(u64, u64)` containing the absolute genomic start and end coordinates of the ORF.
///
/// # Panics
///
/// This function will panic if `get_pos_in_exons` returns `None` for either `orf_start`
/// or `orf_end`, indicating that the ORF coordinates are outside the exonic sequence.
///
pub fn get_cds_from_pos(record: &GenePred, orf_start: u64, orf_end: u64) -> (u64, u64) {
    let cds_start = get_pos_in_exons(record, orf_start).unwrap_or_else(|| {
        panic!("ERROR: could not map ORF_START to exonic coordinates: {orf_start} -> {record:?}")
    });

    let cds_end = get_end_pos_in_exons(record, orf_end).unwrap_or_else(|| {
        panic!("ERROR: could not map ORF_END to exonic coordinates: {orf_end} -> {record:?}")
    });

    match record.strand() {
        Some(genepred::Strand::Forward) => (cds_start, cds_end),
        Some(genepred::Strand::Reverse) => (SCALE - cds_end, SCALE - cds_start),
        Some(genepred::Strand::Unknown) | None => {
            panic!("ERROR: unexpected strand value: {:?}", record.strand())
        }
    }
}

/// Maps ORF transcript coordinates to absolute genomic coordinates (strand-aware).
///
/// This function takes ORF coordinates defined in **dataar transcript space**
/// (as produced by a CDS/ORF predictor) and returns the corresponding genomic
/// start and end coordinates of the coding region, accounting for exon structure
/// and strand orientation. For transcripts on the reverse strand, it assumes
/// coordinates have been **strand-normalized** using a fixed scale (`SCALE`) and
/// reverses them appropriately.
///
/// # Arguments
///
/// * `orf_start` - The start of the ORF in transcript (spliced) coordinates.
/// * `orf_end` - The end of the ORF in transcript (spliced) coordinates. This is exclusive.
///
/// # Returns
///
/// * `(u64, u64)` - A tuple of genomic start and end coordinates corresponding
///   to the CDS region. Returns `(start, end)` if mapping succeeds; otherwise,
///   may return `Error` if the coordinates do not map within the exon structure.
///
/// # Strand Behavior
///
/// * On the `Strand::Forward`, transcript coordinates map left-to-right across exons.
/// * On the `Strand::Reverse`, the exons are reversed and mapped using:
///   `SCALE - (exon_end - offset)`, where `SCALE` is a user-defined upper genomic bound
///   used for coordinate normalization.
///
/// # Example
///
/// ```rust, ignore
/// let mut gene = GenePred { /* initialized */ };
/// let orf_start = 15;
/// let orf_end = 45;
/// if let Some((cds_start, cds_end)) = gene.map_absolute_cds(orf_start, orf_end) {
///     println!("Genomic CDS coordinates: {} - {}", cds_start, cds_end);
/// }
/// ```
pub fn map_absolute_cds(record: &GenePred, orf_start: u64, orf_end: u64) -> (u64, u64) {
    // INFO: for Reverse this should be interpreted from 5' to 3'
    let (cds_start, cds_end) = get_cds_from_pos(record, orf_start, orf_end);

    assert!(
        cds_start > 0,
        "ERROR: {orf_start} is not bigger than 0! -> {record:?}"
    );
    assert!(
        cds_end > 0,
        "ERROR: {orf_end} is not bigget than 0! -> {record:?}"
    );

    assert!(
        cds_start >= record.start,
        "ERROR: {cds_start} is less than {}! Mapped ORF values {orf_start} - {orf_end} -> {record:?}",
        record.start
    );
    assert!(
        cds_end <= record.end,
        "ERROR: {cds_end} is greater than {}! Mapped ORF values {orf_start} - {orf_end} -> {record:?}",
        record.end
    );

    (cds_start, cds_end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use genepred::{Bed12, Reader};

    fn status_of(cmd: &str) -> std::process::ExitStatus {
        std::process::Command::new("bash")
            .arg("-c")
            .arg(cmd)
            .status()
            .unwrap()
    }

    #[test]
    fn test_run_command_accepts_success() {
        run_command("true", "true");
    }

    #[test]
    fn test_exit_code_keeps_the_command_code() {
        assert_eq!(exit_code(&status_of("true")), Ok(()));
        assert_eq!(exit_code(&status_of("exit 3")), Err(3));
        // a tool killed by SIGKILL (the OOM killer) -> 137, which xorf retries
        assert_eq!(exit_code(&status_of("sh -c 'kill -9 $$'")), Err(137));
    }

    #[test]
    fn test_exit_code_of_a_killed_shell() {
        assert_eq!(exit_code(&status_of("kill -9 $$")), Err(137));
    }

    #[test]
    fn test_absolute_cds_mapping_forward() {
        let data = "chr6\t8259278\t8593709\tR441_chr6__FC23#TC31#PA0#PR0#IY1000\t60\t+\t8259278\t8593709\t43,118,219\t10\t215,109,62,215,152,87,117,153,211,2165\t0,11194,114627,167692,278538,299221,313903,320308,323301,332266";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 237;
        let orf_end = 420;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start, orf_end);

        assert_eq!(predicted_cds_start, 8270494);
        assert_eq!(predicted_cds_end, 8427004);
    }

    #[test]
    fn test_absolute_cds_mapping_reverse() {
        let data = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1284;
        let orf_end = 3321;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start, orf_end);
        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 45503698);
        assert_eq!(predicted_cds_end, 45515991);
    }

    #[test]
    fn test_get_cds_from_pos_forward() {
        let data = "chr12\t102521030\t102531285\tR909465_chr12__FC30#TC0#PA164#PR164#IY896\t60\t+\t102521030\t102531285\t43,118,219\t8\t399,47,94,69,159,417,476,417\t0,1113,3694,4517,6573,6895,7826,9838";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1530;
        let orf_end = 1674;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        assert_eq!(predicted_cds_start, 102529201);
        assert_eq!(predicted_cds_end, 102530881);
    }

    #[test]
    fn test_get_pos_in_exons_forward() {
        let data = "chr12\t102521030\t102531285\tR909465_chr12__FC30#TC0#PA164#PR164#IY896\t60\t+\t102521030\t102531285\t43,118,219\t8\t399,47,94,69,159,417,476,417\t0,1113,3694,4517,6573,6895,7826,9838";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1530;
        let orf_end = 1674;

        let predicted_cds_start = get_pos_in_exons(&record, orf_start);
        let predicted_cds_end = get_pos_in_exons(&record, orf_end);

        assert_eq!(predicted_cds_start, Some(102529201));
        assert_eq!(predicted_cds_end, Some(102530881));
    }

    #[test]
    fn test_get_cds_from_pos_reverse() {
        let data = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 358;
        let orf_end = 553;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 45519839);
        assert_eq!(predicted_cds_end, 45520034);
    }

    #[test]
    fn test_get_pos_in_exons_reverse() {
        let data = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1284;
        let orf_end = 3321;

        let predicted_cds_end = get_pos_in_exons(&record, orf_start).unwrap();
        let predicted_cds_start = get_pos_in_exons(&record, orf_end);

        dbg!(
            SCALE - predicted_cds_start.unwrap(),
            SCALE - predicted_cds_end
        );

        assert_eq!(SCALE - predicted_cds_start.unwrap(), 45503698);
        assert_eq!(SCALE - predicted_cds_end, 45515991);
    }

    #[test]
    fn test_get_pos_in_exons_forward_refseq_tai_prediction() {
        let data = "chr1\t58713285\t58733227\tNM_009805.4\t0\t+\t58726436\t58732362\t0,0,0\t5\t374,427,106,136,975,\t0,13020,15770,17866,18967,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 505;
        let orf_end = 1150;

        let predicted_cds_start = get_pos_in_exons(&record, orf_start).unwrap();
        let predicted_cds_end = get_pos_in_exons(&record, orf_end).unwrap();

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 58726436);
        assert_eq!(predicted_cds_end + 3, 58732362);
    }

    // #[test]
    // fn test_get_pos_in_exons_reverse_refseq_tai_prediction() {
    //     let data = "chr2\t73092800\t73214447\tNM_025942.2\t0\t-\t73093622\t73212960\t0,0,0\t11\t924,123,97,141,98,81,176,128,144,101,162,\t0,4344,6491,7273,49097,49507,63937,106600,110524,120059,121485,";
    //
    //     let mut reader: Reader<Bed12> =
    //         Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
    //     let record = reader.next().unwrap().unwrap();
    //
    //     let orf_start = 162;
    //     let orf_end = 1350;
    //
    //     let predicted_cds_end = get_pos_in_exons(&record, orf_start).unwrap();
    //     let predicted_cds_start = get_pos_in_exons(&record, orf_end);
    //
    //     dbg!(
    //         SCALE - predicted_cds_start.unwrap(),
    //         SCALE - predicted_cds_end
    //     );
    //
    //     assert_eq!(SCALE - predicted_cds_start.unwrap() - 3, 73093622);
    //     assert_eq!(SCALE - predicted_cds_end, 73212960);
    // }
    #[test]
    fn test_get_pos_in_exons_reverse_refseq_tai_prediction_additional() {
        // 2025-09-10T11:33:47.207Z WARN  [orf::tai] WARN: translationAi predicted a non-stop ORF: "304,484,0.11812351644039154,0.1313595026731491" for GenePred { nam
        // e: "9249", chrom: "chr7", strand: Reverse, start: 99940345148, end: 99940364371, cds_start: 99940345148, cds_end: 99940364371, exons: [(99940345148, 999403
        // 45217), (99940350839, 99940350935), (99940354698, 99940354847), (99940361505, 99940361651), (99940362689, 99940362841), (99940363369, 99940363481), (999403
        // 64269, 99940364371)], introns: [(99940345218, 99940350838), (99940350936, 99940354697), (99940354848, 99940361504), (99940361652, 99940362688), (9994036284
        // 2, 99940363368), (99940363482, 99940364268)], exon_len: 826, exon_count: 7, data: "chr7\t59635629\t59654852\t9249\t60\t-\t59635629\t59654852\t43,118,219\t7
        // \t102,112,152,146,149,96,69\t0,890,1530,2720,9524,13436,19154", is_ref: false }

        let data = "chr7\t59635629\t59654852\t9249\t60\t-\t59635629\t59654852\t43,118,219\t7\t102,112,152,146,149,96,69\t0,890,1530,2720,9524,13436,19154";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 317;
        let orf_end = 437;

        let (predicted_cds_end, predicted_cds_start) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        if predicted_cds_start - 3 < record.start {
            dbg!("WARN: translationAi predicted a non-stop ORF: {orf_start:?} for {?}");
        }

        if predicted_cds_end + 3 > record.end {
            dbg!("WARN: translationAi predicted a non-stop ORF: {orf_start:?} for {?}");
        }

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start - 3, 59638489);
        assert_eq!(predicted_cds_end, 59638372);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_bug_prediction_additional_2() {
        let data = "scaffold_10\t60508574\t60652099\tR11900_scaffold_10__FC0#TC0#PA29#PR30#IY979#SG\t60\t-\t60508574\t60652099\t51,153,255\t4\t499,38,310,77\t0,4353,112181,143448";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 28;
        let orf_end = 748;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 60508750);
        assert_eq!(predicted_cds_end, 60652071);

        let orf_start = 52;
        let orf_end = 748;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 60508750);
        assert_eq!(predicted_cds_end, 60652047);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_bug_prediction_additional_3() {
        let data = "scaffold_17\t80102749\t80103945\tR16258_scaffold_17__FC0#TC0#PA2#PR28#IY971\t0\t+\t80102749\t80103945\t43,118,219\t1\t1196\t0";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 84;
        let orf_end = 483;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 80102833);
        assert_eq!(predicted_cds_end, 80103232);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_refseq_orfipy_prediction() {
        let data = "chr5\t151561660\t151586924\tNM_001102582.1\t0\t-\t151561660\t151586906\t0,0,0\t5\t911,124,228,836,286,\t0,11135,14005,22893,24978,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_end = 2385;

        let predicted_cds_start = get_pos_in_exons(&record, orf_end);

        dbg!(SCALE - predicted_cds_start.unwrap());

        assert_eq!(SCALE - predicted_cds_start.unwrap(), 151561660);
    }

    #[test]
    fn test_get_pos_in_exons_forward_refseq_orfipy_prediction() {
        let data = "chr1\t83116765\t83119167\tNM_001159738.1\t0\t+\t83116823\t83118717\t0,0,0\t4\t137,112,78,472,\t0,1033,1233,1930,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 46;
        let orf_end = 349;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 83116811);
        assert_eq!(predicted_cds_end, 83118717);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow() {
        let data = "NC_037328	102267973	102388511	XM_005201715.5	0	+	102272281	102388043	0,0,0	48	96,118,95,118,110,152,172,100,113,126,132,120,114,85,118,172,117,155,85,57,125,89,50,171,156,204,155,169,97,113,126,129,111,63,117,88,139,134,35,117,149,85,57,125,89,50,168,537	0,4308,5644,7568,7981,9021,10662,12452,13630,15870,18181,18399,20348,29163,32745,34854,38566,40588,41701,42678,43435,45403,48668,50006,58153,68412,70214,71796,73546,73750,73952,75545,77369,79915,81921,83444,84463,89302,90352,90474,93758,96861,97028,97687,100948,116510,117018,120001";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 96;
        let orf_end = 5535;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 102272281);
        assert_eq!(predicted_cds_end, 102388043);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow_additional() {
        let data = "NC_037328	127530720	127609658	XM_059888940	0	-	127533976	127537591	0,0,0	7	6871,105,11028,446,66,75,309	0,44010,46138,58325,59275,77446,78629";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 12029;
        let orf_end = 15644;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 127533976);
        assert_eq!(predicted_cds_end, 127537591);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow_additional_2() {
        let data = "NC_037328	69329582	69357506	XM_002684796.7	0	-	69332533	69357372	0,0,0	12	923,196,85,38,121,126,164,56,107,120,1074,186	0,2951,5733,7446,8130,8564,11264,14986,18831,20186,24458,27738";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 134;
        let orf_end = 2273;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 69332533);
        assert_eq!(predicted_cds_end, 69357372);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow_additional_3() {
        let data = "NC_037357	35186707	35245555	XM_024988460.2	0	+	35201872	35242092	0,0,0	22	158,79,66,42,78,90,87,84,42,87,105,99,96,90,195,195,111,81,72,103,99,2534	0,12066,12624,14433,14577,15164,16087,16869,17192,23711,26945,28746,31210,32094,38359,43644,46231,47344,49432,50523,55286,56314";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 424;
        let orf_end = 2059;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 35201872);
        assert_eq!(predicted_cds_end, 35242092);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow_additional_4() {
        let data = "NC_037335	111666383	112087409	XM_024995855.2	0	+	111926245	112085362	0,0,0	25	240,88,106,143,154,55,34,63,350,975,135,91,108,212,251,234,122,69,63,84,222,92,104,144,2074	0,38344,141165,225665,255268,259862,260086,260311,282916,297251,303053,307317,308931,312650,319990,322543,324983,325986,385391,391222,392566,403274,413083,417769,418952";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 731;
        let orf_end = 4166;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 111926245);
        assert_eq!(predicted_cds_end, 112085362);
    }
}
