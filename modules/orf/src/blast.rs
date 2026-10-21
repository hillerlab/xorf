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

use genepred::GenePred;
use hashbrown::{hash_map::Entry, HashMap};
use log::{debug, info, warn};
use memchr::memchr;
use memmap2::Mmap;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::str::{from_utf8, FromStr};
use std::sync::Arc;

use crate::{
    cli::{BlastArgs, Engine},
    consts::*,
    utils::*,
};

const ESMFOLD2_FAST: &str = include_str!("esmfold2_fast.py");
const ESMFOLD2_MAX_LEN: usize = 2048;

/// Run diamont on the a nested orfipy input + translationAi fasta
///
/// # Arguments
///
/// * `args` - The arguments for running the pipeline
///
/// # Returns
///
/// * `None`
///
/// # Example
///
/// ```rust
/// let args = Args::parse();
/// run_blast(args);
/// ```
pub fn run_blast(args: BlastArgs, threads: usize) {
    let dir = args.outdir.join(format!("{}.orf", args.prefix));
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let fasta = match args.upstream_flank > 0 || args.downstream_flank > 0 {
        true => cut_fasta(
            &args.fasta,
            args.upstream_flank,
            args.downstream_flank,
            &dir,
        ),
        false => args.fasta,
    };

    let pep = run_orfipy(&fasta, &dir);
    let bed = get_bed(&args.bed);

    let orfs = deduplicate(
        &pep,
        true,
        args.orf_min_len,
        args.orf_min_percent,
        "M".as_bytes(),
    );

    let table = get_table(
        orfs,
        bed,
        &dir,
        args.tai,
        args.net,
        args.nmd_distance,
        args.weak_nmd_distance,
        args.atg_distance,
        args.big_exon_dist_to_ej,
        args.prefix.clone(),
    );

    __run_diamond_and_psauron(
        table,
        args.database,
        &dir,
        args.prefix,
        args.engine,
        args.esm,
        threads,
    );

    if !args.keep_temp {
        println!("INFO: removing temporary files in {dir:?}");
        // INFO: remove everything under outdir that does not end with '.result'
        for entry in std::fs::read_dir(&dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_file()
                && !path
                    .file_name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .ends_with(".result")
            {
                std::fs::remove_file(path).unwrap();
            }
        }
    } else {
        println!("INFO: keeping temporary files in {dir:?}");
    }
}

/// Cuts the input fasta file by the given flanks
///
/// # Arguments
///
/// * `fasta` - The path to the input fasta file
/// * `upstream_flank` - The number of bases upstream of the ORF to include in the context sequence
/// * `downstream_flank` - The number of bases downstream of the ORF to include in the context sequence
/// * `outdir` - The output directory
///
/// # Returns
///
/// A `PathBuf` containing the path to the cut fasta file
///
/// # Example
///
/// ```rust, ignore
/// let fasta = cut_fasta(&args.fasta, args.upstream_flank, args.downstream_flank, &dir);
/// ```
fn cut_fasta(
    fasta: &PathBuf,
    upstream_flank: usize,
    downstream_flank: usize,
    outdir: &Path,
) -> PathBuf {
    let file = File::open(fasta).unwrap();
    let mmap = unsafe { Mmap::map(&file).unwrap() };
    let data = mmap.as_ref();

    let mut pos = 0;

    let out_fasta = outdir.join("tmp.fa");
    let mut writer = BufWriter::new(
        File::create(&out_fasta).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
    );

    while let Some(start) = memchr(b'>', &data[pos..]) {
        let start = pos + start;
        let end = memchr(b'>', &data[start + 1..]).map_or(data.len(), |e| start + 1 + e);
        let entry = &data[start + 1..end];
        let header_end = memchr(b'\n', entry).unwrap();
        let header = from_utf8(&entry[..header_end]).unwrap().trim();
        let record = &entry[header_end + 1..];
        let seq = record
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r') // Remove newlines and carriage returns
            .cloned()
            .collect::<Vec<u8>>();

        if seq.len() < upstream_flank + downstream_flank {
            panic!(
                "ERROR: sequence is too short for upstream_flank + downstream_flank -> {} + {} = {} > {} for {}",
                upstream_flank,
                downstream_flank,
                upstream_flank + downstream_flank,
                seq.len(),
                header
            );
        }

        let chopped = seq[upstream_flank..seq.len() - downstream_flank].to_vec();

        writer
            .write_all(format!(">{}\n", header).as_bytes())
            .unwrap_or_else(|e| {
                panic!("ERROR: failed to write to file -> {e}");
            });
        writer.write_all(&chopped).unwrap_or_else(|e| {
            panic!("ERROR: failed to write to file -> {e}");
        });
        writer.write_all(b"\n").unwrap_or_else(|e| {
            panic!("ERROR: failed to write to file -> {e}");
        });

        pos = end;
    }

    out_fasta
}

/// Run diamond on the orfipy output
///
/// # Arguments
///
/// * `table` - The orfipy output
/// * `database` - The path to the translationAi fasta
/// * `outdir` - The output directory
///
/// # Returns
///
/// * `None`
///
/// # Example
///
/// ```rust
/// let table = get_table(orfs, bed, &dir, args.tai);
/// __run_diamond(table, DATABASE, &dir);
/// ```
fn __run_diamond_and_psauron(
    mut table: HashMap<usize, Vec<String>>,
    database: PathBuf,
    outdir: &Path,
    prefix: String,
    engine: Engine,
    esm: bool,
    threads: usize,
) {
    let blasting = match engine {
        Engine::Diamond => outdir.join(format!("{}.orf.diamond", prefix)),
        Engine::Mmseqs2 => outdir.join(format!("{}.orf.mmseqs2", prefix)),
    };
    let orfs = outdir.join(format!("{}.orf.dedup.pep", prefix));
    let psauron = outdir.join(format!("{}.orf.psauron", prefix));
    let mmseqs_tmp = outdir.join(format!("{}.mmseqs.tmp", prefix));

    if engine == Engine::Mmseqs2
        && database
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("dmnd"))
    {
        panic!(
            "ERROR: mmseqs2 needs a FASTA or mmseqs database, not a diamond .dmnd file -> {database:?}"
        );
    }

    if engine == Engine::Mmseqs2 {
        std::fs::create_dir_all(&mmseqs_tmp).unwrap_or_else(|e| {
            panic!("ERROR: could not create mmseqs tmp dir {mmseqs_tmp:?} -> {e}")
        });
    }

    let cmd = match engine {
        Engine::Diamond => format!(
            "diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid pident qlen slen length qstart qend sstart send evalue sseqid --threads {} --sensitive -e 1e-10",
            orfs.display(),
            database.display(),
            blasting.display(),
            threads
        ),
        Engine::Mmseqs2 => format!(
            "mmseqs easy-search {} {} {} {} --threads {} --search-type 1 -e 1e-10 --sort-results 1 --remove-tmp-files 1 --format-output query,pident,qlen,tlen,alnlen,qstart,qend,tstart,tend,evalue,target",
            orfs.display(),
            database.display(),
            blasting.display(),
            mmseqs_tmp.display(),
            threads
        ),
    };

    println!("INFO: Executing -> {}", cmd);

    run_command(&cmd, &engine.to_string());

    if engine == Engine::Mmseqs2 {
        let _ = std::fs::remove_dir_all(&mmseqs_tmp);
    }

    let cmd = format!("psauron -i {} -o {} -p", orfs.display(), psauron.display());
    println!("INFO: Executing -> {}", cmd);

    run_command(&cmd, "psauron");

    let mut seen = HashSet::new();

    // WARN: breaking point -> if BLAST is empty we throw an error; very plausible to happen
    // WARN: we now avoid panicking and instead return an empty table
    let predictions = reader(&blasting).unwrap_or_default();
    if predictions.is_empty() {
        println!(
            "WARN: no predictions found for {blasting:?}. Will proceed with empty BLAST table"
        );
    }

    let psauron_predictions = read_psauron(psauron);
    let esm_predictions = esm.then(|| run_esmfold2_fast(&orfs, outdir, &prefix, threads, &table));

    let mut writer = BufWriter::new(
        File::create(blasting.with_extension(RESULT)).unwrap_or_else(|e| {
            panic!("ERROR: failed to create output file for blast results -> {e}");
        }),
    );

    for line in predictions.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        let header = parts[0];
        let u_header = header.parse::<usize>().unwrap_or_else(|_| {
            panic!("ERROR: failed to parse header as usize -> {header} in parts: {parts:?}")
        });

        // INFO: only storing first row -> sorted by default
        if seen.contains(header) {
            continue;
        }

        seen.insert(header);

        let lines = table.get(&u_header).unwrap_or_else(|| {
            panic!("ERROR: could not find header {header} in table -> {parts:?}!")
        });
        let psauron_score = psauron_predictions.get(&u_header).unwrap_or_else(|| {
            panic!("ERROR: could not find header {header} in psauron predictions -> {parts:?}")
        });

        for line in lines {
            let blast = BlastRecord::from_parts(&parts);
            let output = match &esm_predictions {
                Some(scores) => format!(
                    "{}\t{:.2}\t{:e}\t{}\t{}\t{}\t{}\t{:.4}\t{}\n",
                    line,
                    blast.blast_pid,
                    blast.blast_e_value,
                    blast.blast_offset,
                    blast.blast_alignment_len,
                    blast.percent_aligned,
                    psauron_score,
                    scores[&u_header],
                    blast.reference_id
                ),
                None => format!(
                    "{}\t{:.2}\t{:e}\t{}\t{}\t{}\t{}\t{}\n",
                    line,
                    blast.blast_pid,
                    blast.blast_e_value,
                    blast.blast_offset,
                    blast.blast_alignment_len,
                    blast.percent_aligned,
                    psauron_score,
                    blast.reference_id
                ),
            };

            writer.write_all(output.as_bytes()).unwrap_or_else(|e| {
                panic!("ERROR: failed to write to file -> {e}");
            });
        }

        // INFO: table records with no blast hits remain
        table.remove(&u_header);
    }

    for (u_header, lines) in table.iter() {
        if let Some(psauron_score) = psauron_predictions.get(u_header) {
            for line in lines {
                let output = match &esm_predictions {
                    Some(scores) => format!(
                        "{}\t0\t1\t0\t0\t0\t{}\t{:.4}\tUNKNOWN\n",
                        line, psauron_score, scores[u_header]
                    ),
                    None => format!("{}\t0\t1\t0\t0\t0\t{}\tUNKNOWN\n", line, psauron_score),
                };

                writer.write_all(output.as_bytes()).unwrap_or_else(|e| {
                    panic!("ERROR: failed to write to file -> {e}");
                });
            }
        } else {
            for line in lines {
                let output = match &esm_predictions {
                    Some(scores) => format!(
                        "{}\t0\t1\t0\t0\t0\t0.0\t{:.4}\tUNKNOWN\n",
                        line, scores[u_header]
                    ),
                    None => format!("{}\t0\t1\t0\t0\t0\t0.0\tUNKNOWN\n", line),
                };

                writer.write_all(output.as_bytes()).unwrap_or_else(|e| {
                    panic!("ERROR: failed to write to file -> {e}");
                });
            }
        }
    }
}

fn sanitize_esmfold2_sequence(sequence: &str) -> String {
    sequence
        .trim_end_matches('*')
        .chars()
        .map(|aa| if aa == '*' { 'X' } else { aa })
        .collect()
}

fn parse_numeric_fasta(contents: &str) -> Result<Vec<(usize, String)>, String> {
    let mut records = Vec::new();
    let mut current_id = None;
    let mut sequence = String::new();

    for raw_line in contents.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(seq_idx) = current_id.take() {
                if sequence.is_empty() {
                    return Err(format!("empty sequence for FASTA ID {seq_idx}"));
                }
                records.push((seq_idx, std::mem::take(&mut sequence)));
            }
            current_id = Some(
                header
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| "empty FASTA header".to_string())?
                    .parse::<usize>()
                    .map_err(|_| format!("non-numeric FASTA header: {header}"))?,
            );
        } else if current_id.is_some() {
            sequence.push_str(line);
        } else {
            return Err("sequence before first FASTA header".to_string());
        }
    }

    if let Some(seq_idx) = current_id {
        if sequence.is_empty() {
            return Err(format!("empty sequence for FASTA ID {seq_idx}"));
        }
        records.push((seq_idx, sequence));
    }

    Ok(records)
}

fn parse_esmfold2_scores(contents: &str) -> Result<HashMap<usize, f32>, String> {
    let mut scores = HashMap::new();
    for line in contents.lines().filter(|line| !line.trim().is_empty()) {
        let mut fields = line.split('\t');
        let seq_idx = fields
            .next()
            .ok_or_else(|| format!("missing sequence ID in score row: {line}"))?
            .parse::<usize>()
            .map_err(|_| format!("non-numeric sequence ID in score row: {line}"))?;
        let score = fields
            .next()
            .ok_or_else(|| format!("missing pLDDT in score row: {line}"))?
            .parse::<f32>()
            .map_err(|_| format!("non-numeric pLDDT in score row: {line}"))?;
        if fields.next().is_some() {
            return Err(format!("extra columns in score row: {line}"));
        }
        if !score.is_finite() || !(0.0..=100.0).contains(&score) {
            return Err(format!("pLDDT outside 0..100 in score row: {line}"));
        }
        if scores.insert(seq_idx, score).is_some() {
            return Err(format!("duplicate sequence ID in score row: {seq_idx}"));
        }
    }
    Ok(scores)
}

fn run_esmfold2_fast(
    orfs: &Path,
    outdir: &Path,
    prefix: &str,
    threads: usize,
    table: &HashMap<usize, Vec<String>>,
) -> HashMap<usize, f32> {
    let input = outdir.join(format!("{prefix}.orf.esmfold2.pep"));
    let output = outdir.join(format!("{prefix}.orf.esmfold2.tsv"));
    let script = outdir.join(format!("{prefix}.orf.esmfold2_fast.py"));
    let records = parse_numeric_fasta(
        &reader(orfs).unwrap_or_else(|e| panic!("ERROR: failed to read ESMFold2 input -> {e}")),
    )
    .unwrap_or_else(|e| panic!("ERROR: invalid deduplicated peptide FASTA -> {e}"));
    let mut writer = BufWriter::new(
        File::create(&input)
            .unwrap_or_else(|e| panic!("ERROR: failed to create ESMFold2 input -> {e}")),
    );
    let mut scores = HashMap::new();
    let mut foldable = 0usize;
    let mut residues = 0usize;

    for (seq_idx, sequence) in &records {
        let sequence = sanitize_esmfold2_sequence(sequence);
        if sequence.len() > ESMFOLD2_MAX_LEN {
            let representative = table
                .get(seq_idx)
                .and_then(|lines| lines.first())
                .and_then(|line| line.split('\t').nth(3))
                .unwrap_or("UNKNOWN");
            warn!(
                "ESMFold2-Fast skipped {representative} (sequence {seq_idx}, {} aa > {ESMFOLD2_MAX_LEN}); assigning 0",
                sequence.len()
            );
            scores.insert(*seq_idx, 0.0);
            continue;
        }
        if sequence.is_empty() {
            panic!("ERROR: ESMFold2 sequence {seq_idx} is empty after sanitization");
        }
        writeln!(writer, ">{seq_idx}\n{sequence}")
            .unwrap_or_else(|e| panic!("ERROR: failed to write ESMFold2 input -> {e}"));
        foldable += 1;
        residues += sequence.len();
    }
    drop(writer);

    let started = std::time::Instant::now();
    if foldable > 0 {
        std::fs::write(&script, ESMFOLD2_FAST)
            .unwrap_or_else(|e| panic!("ERROR: failed to write ESMFold2 helper -> {e}"));
        let status = std::process::Command::new("python")
            .arg(&script)
            .arg(&input)
            .arg(&output)
            .arg(threads.to_string())
            .status()
            .unwrap_or_else(|e| panic!("ERROR: failed to start ESMFold2-Fast -> {e}"));
        if !status.success() {
            panic!("ERROR: ESMFold2-Fast exited with status {status}");
        }

        let inferred = parse_esmfold2_scores(
            &reader(&output)
                .unwrap_or_else(|e| panic!("ERROR: failed to read ESMFold2 scores -> {e}")),
        )
        .unwrap_or_else(|e| panic!("ERROR: invalid ESMFold2 score output -> {e}"));
        for (seq_idx, score) in inferred {
            if scores.insert(seq_idx, score).is_some() {
                panic!("ERROR: duplicate ESMFold2 score for sequence {seq_idx}");
            }
        }
    }

    for seq_idx in table.keys() {
        if !scores.contains_key(seq_idx) {
            panic!("ERROR: missing ESMFold2 score for sequence {seq_idx}");
        }
    }
    if scores.len() != table.len() {
        panic!(
            "ERROR: ESMFold2 returned {} scores for {} deduplicated sequences",
            scores.len(),
            table.len()
        );
    }

    let elapsed = started.elapsed().as_secs_f64();
    let rate = if elapsed > 0.0 {
        foldable as f64 / elapsed
    } else {
        0.0
    };
    let residue_rate = if elapsed > 0.0 {
        residues as f64 / elapsed
    } else {
        0.0
    };
    println!(
        "INFO: ESMFold2-Fast summary: {} deduplicated, {} folded, {} skipped, {:.2}s, {:.2} proteins/s, {:.2} residues/s",
        records.len(),
        foldable,
        records.len() - foldable,
        elapsed,
        rate,
        residue_rate
    );
    scores
}

/// Reads the psauron output and converts it to a hashmap of `PsauronRecord`s
///
/// # Arguments
///
/// * `psauron` - The path to the psauron output
///
/// # Returns
///
/// A `HashMap` of `String` and their corresponding `PsauronRecord`
///
/// # Example
///
/// ```rust
/// let psauron = read_psauron(psauron).unwrap();
/// ```
pub fn read_psauron<P: AsRef<Path> + std::fmt::Debug>(psauron: P) -> HashMap<usize, f32> {
    let predictions = reader(&psauron).unwrap_or_default();
    // .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

    let mut accumulator = HashMap::new();

    // INFO: format below
    // /psauron/bin/psauron -i orf.dedup.pep -o test.csv -p
    // psauron score: 70.6
    // description,psauron_is_protein,in-frame_score
    // 64,True,0.99869
    // 19,True,0.99216
    if predictions.is_empty() {
        log::warn!("WARN: no predictions found for PSAURON in {psauron:?}");
        return accumulator;
    } else {
        for line in predictions.lines().skip(3) {
            let record = PsauronRecord::parse(line);
            accumulator.insert(record.u_id, record.psauron_score);
        }
    }

    accumulator
}

/// A struct representing a psauron record
#[derive(Debug, Clone)]
pub struct PsauronRecord {
    pub u_id: usize,
    pub psauron_score: f32,
}

impl PsauronRecord {
    /// Parses a psauron record from a string
    ///
    /// # Arguments
    ///
    /// * `line` - The string to parse
    ///
    /// # Returns
    ///
    /// A `PsauronRecord` struct
    ///
    /// # Panics
    /// This function will panic if the line is not in the correct format
    ///
    /// # Example
    ///
    /// ```rust
    /// let line = "64,True,0.99869";
    /// let record = PsauronRecord::parse(line);
    /// ```
    pub fn parse(line: &str) -> Self {
        let mut parts = line.split(',');

        let u_id = parts
            .next()
            .unwrap_or_else(|| {
                panic!("ERROR: failed to parse u_id from line: {}", line);
            })
            .parse::<usize>()
            .unwrap_or_else(|_| {
                panic!("ERROR: failed to parse header as usize in parts: {parts:?}")
            });

        parts.next();

        let psauron_score = parts
            .next()
            .unwrap_or_else(|| {
                panic!("ERROR: failed to parse psauron_score from line: {}", line);
            })
            .parse::<f32>()
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: failed to parse psauron_score from line: {} -> {e}",
                    line
                );
            });

        Self {
            u_id,
            psauron_score,
        }
    }
}

/// Represents a single BLAST alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BlastRecord {
    pub blast_pid: f32,           // Percentage of identical matches
    pub blast_e_value: f64,       // E-value of the match
    pub blast_offset: i32,        // Offset in the query sequence where the match starts
    pub blast_alignment_len: u32, // Length of the alignment
    pub percent_aligned: f32,     // Percentage of the query sequence that is aligned
    pub reference_id: String,     // Reference ID of the match
}

impl BlastRecord {
    /// Creates a new `BlastRecord` from a slice of string parts, typically
    /// obtained by splitting a line from a DIAMOND BLAST output.
    ///
    /// The expected format of `parts` corresponds to `diamond blastp --outfmt 6`
    /// output: `qseqid pident qlen slen length qstart qend sstart send evalue`.
    ///
    /// # Arguments
    ///
    /// * `parts` - A slice of string slices, where each element represents a column from the BLAST output.
    ///
    /// # Returns
    ///
    /// A `BlastRecord` instance populated with the parsed data. The `blast_id`
    /// field is initially empty and is expected to be set later.
    ///
    /// # Panics
    ///
    /// This function will panic if:
    /// - `parts` does not contain at least 10 elements.
    /// - Any of the numeric fields (`blast_idx_id`, `blast_pid`, `blast_e_value`,
    ///   `blast_offset` components, `blast_alignment_len`, query length for `percent_aligned`)
    ///   cannot be successfully parsed into their respective types.
    ///
    /// # Example
    ///
    /// Follows this format:
    ///
    /// qseqid pident  qlen    slen   length qstart    qend   sstart   send     evalue     sseqid
    ///  17      97.2    142     357     141     1       141     217     357     5.09e-93  sp|O55165|KIF3C_RAT
    ///
    /// ```rust
    /// let parts = ["1", "99.0", "500", "0", "100", "1", "100", "1", "100", "1e-10"];
    /// let record = BlastRecord::from_parts(&parts);
    /// ```
    pub fn from_parts(parts: &[&str]) -> Self {
        if parts.len() < 10 {
            panic!("ERROR: not enough parts to create BlastRecord -> {parts:?}");
        }

        let blast_pid = parts[1].parse::<f32>().unwrap_or_else(|_| {
            panic!("ERROR: failed to parse blast PID from parts: {:?}", parts);
        });

        let blast_e_value = parts[9].parse::<f64>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast E-value from parts: {:?}",
                parts
            );
        });
        // INFO: if parsed to zero, but string was not "0.0", it's subnormal
        let blast_e_value = if blast_e_value == 0.0 {
            // INFO: represent it with the minimum positive value
            f64::MIN_POSITIVE // INFO: ~2.225074e-308
        } else {
            blast_e_value
        };

        let qstart = parts[5].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            );
        });
        let qend = parts[6].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            );
        });

        let blast_offset = parts[7].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        }) - qstart;

        let blast_alignment_len = parts[4].parse::<u32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        });

        let percent_aligned = (qend - qstart + 1) as f32
            / parts[2].parse::<u32>().unwrap_or_else(|_| {
                panic!(
                    "ERROR: failed to parse blast length from parts: {:?}",
                    parts
                );
            }) as f32
            * 100.0;

        let reference_id = parts[10]
            .split("|")
            .last()
            .unwrap_or_else(|| {
                panic!("ERROR: failed to parse reference id from blast line: {parts:?}",);
            })
            .to_string();

        Self {
            blast_pid,
            blast_e_value,
            blast_offset,
            blast_alignment_len,
            percent_aligned,
            reference_id,
        }
    }
}

/// Run orfipy on the input fasta
///
/// # Arguments
///
/// * `fasta` - The path to the input fasta file
/// * `dir` - The output directory
///
/// # Returns
///
/// * `PathBuf` - The path to the orfipy output
///
/// # Example
///
/// ```rust
/// let dir = args.outdir.join("orf");
/// std::fs::create_dir_all(&dir)
///     .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}"));
///
/// let pep = run_orfipy(&args.fasta, &dir);
/// ```
fn run_orfipy(fasta: &Path, dir: &Path) -> PathBuf {
    let cmd = format!(
        "orfipy {} --pep {} --partial-5 --partial-3 --include-stop --strand f --min 100 --ignore-case --outdir {} --single-mode",
        fasta.display(),
        ORF_PEP,
        &dir.display()
    );

    run_command(&cmd, "orfipy");

    dir.join(ORF_PEP)
}

/// Read a file into a string
///
/// # Arguments
///
/// * `file` - The path to the file
///
/// # Returns
///
/// * `Result<String, Box<dyn std::error::Error>>` - The contents of the file
///
/// # Example
///
/// ```rust
/// let pep = reader(&args.fasta);
/// ```
pub fn reader<P: AsRef<Path> + std::fmt::Debug>(
    file: P,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

/// Represents a single ORF record from orfipy
///
/// # Fields
///
/// * `id` - The ID of the ORF
/// * `start` - The start position of the ORF
/// * `end` - The end position of the ORF
/// * `orf_type` - The type of the ORF
/// * `frame` - The frame of the ORF
/// * `start_codon` - The start codon of the ORF
/// * `stop_codon` - The stop codon of the ORF
/// * `strand` - The strand of the ORF
#[derive(Debug, Clone, Eq, Hash, PartialEq)]
pub struct OrfRecord {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub orf_type: OrfType,
    pub frame: i8,
    pub start_codon: String,
    pub stop_codon: String,
    pub strand: char,
    pub seq_idx: usize,
}

/// Splits a header string into its components
///
/// # Arguments
///
/// * `header` - The header string
///
/// # Returns
///
/// * `Option<(usize, usize, char)>` - The start, end, and strand of the ORF
///
/// # Example
///
/// ```rust
/// let (start, end, strand) = split_coords(coords).unwrap_or_else(|| {
///     panic!("ERROR: failed to parse coords from line: {}", line);
/// });
/// ```
pub fn split_coords(header: &str) -> Option<(usize, usize, char)> {
    let mut parts = header.split('-');

    let start = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF start from header: {}", header);
        })
        .strip_prefix('[')
        .unwrap_or_else(|| {
            panic!(
                "ERROR: failed to strip prefix from ORF start from header: {}",
                header
            );
        })
        .parse::<usize>()
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: failed to parse ORF start from header: {} -> {e}",
                header
            );
        });

    parts = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF end from header: {}", header);
        })
        .split(']');

    let end = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF end from header: {}", header);
        })
        .parse::<usize>()
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: failed to parse ORF end from header: {} -> {e}",
                header
            );
        });

    let strand = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF strand from header: {}", header);
        })
        .char_indices()
        .nth(1)
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF strand from header: {}", header);
        })
        .1;

    Some((start, end, strand))
}

impl OrfRecord {
    /// Parses an ORF record from a string
    ///
    /// # Arguments
    ///
    /// * `line` - The string to parse
    ///
    /// # Returns
    ///
    /// * `Self` - The parsed ORF record
    ///
    /// # Example
    ///
    /// ```rust
    /// let record = OrfRecord::parse(line);
    /// ```
    pub fn parse(line: &str) -> Self {
        let mut parts = line.split(' ');

        // WARN: if headers have spaces splitting will fail!
        // INFO: >ENSMUST00000167914.2_ORF.1 [69-330](+) type:complete length:258 frame:1 start:ATG stop:TGA
        let (id, coords, orf_type, _, frame, start_codon, stop_codon) = (
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF ID from line: {}", line);
                })
                .to_string(),
            parts.next().unwrap_or_else(|| {
                panic!(
                    "ERROR: failed to parse ORF coords from line: {} -> {:?}",
                    line, parts
                );
            }),
            OrfType::from_str(parts.next().unwrap_or_else(|| {
                panic!("ERROR: failed to parse ORF type from line: {}", line);
            }))
            .unwrap_or_else(|e| {
                panic!("ERROR: failed to parse ORF type from line: {} -> {e}", line);
            }),
            parts.next(), // INFO: length
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF frame from line: {}", line);
                })
                .strip_prefix("frame:")
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to strip prefix from ORF frame from line: {}",
                        line
                    );
                })
                .parse::<i8>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse ORF frame from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF start codon from line: {}", line);
                })
                .strip_prefix("start:")
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to strip prefix from ORF start codon from line: {}",
                        line
                    );
                })
                .to_string()
                .to_uppercase(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF stop codon from line: {}", line);
                })
                .strip_prefix("stop:")
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to strip prefix from ORF stop codon from line: {}",
                        line
                    );
                })
                .to_string()
                .to_uppercase(),
        );

        let (start, end, strand) = split_coords(coords).unwrap_or_else(|| {
            panic!("ERROR: failed to parse coords from line: {}", line);
        });

        Self {
            id,
            start,
            end,
            orf_type,
            frame,
            start_codon,
            stop_codon,
            strand,
            seq_idx: 0, // INFO: placeholder
        }
    }
}

/// A HashHead is a tuple of (idx, seq) where idx is the index of the hash
#[derive(Debug, Clone, Eq, Hash, PartialEq)]
pub struct HashHead {
    seq: Arc<Vec<u8>>,
}

impl HashHead {
    /// Creates a new `HashHead` from an index and a sequence
    ///
    /// # Arguments
    ///
    /// * `idx` - The index of the hash
    /// * `seq` - The sequence of the hash
    ///
    /// # Returns
    ///
    /// A new `HashHead` instance
    ///
    /// # Example
    ///
    /// ```rust
    /// let hash_head = HashHead::new(0, Arc::from(vec![b'A', b'T', b'G']));
    /// ```
    pub fn new(seq: Arc<Vec<u8>>) -> Self {
        Self { seq }
    }
}

/// An enum representing the type of an ORF
#[derive(Debug, Clone, Copy, Eq, Hash, PartialEq)]
pub enum OrfType {
    Complete,
    ThreePartial,
    FivePartial,
    FivePartialNested,
    ThreePartialNested,
    CompleteNested,
}

impl std::str::FromStr for OrfType {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "type:complete" => Ok(OrfType::Complete),
            "type:3-prime-partial" => Ok(OrfType::ThreePartial),
            "type:5-prime-partial" => Ok(OrfType::FivePartial),
            "type:5-prime-partial-nested" => Ok(OrfType::FivePartialNested),
            "type:3-prime-partial-nested" => Ok(OrfType::ThreePartialNested),
            "type:complete-nested" => Ok(OrfType::CompleteNested),
            _ => Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Invalid ORF type: {}", s),
            ))),
        }
    }
}

impl std::fmt::Display for OrfType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OrfType::Complete => write!(f, "CO"),
            OrfType::ThreePartial => write!(f, "TP"),
            OrfType::FivePartial => write!(f, "FP"),
            OrfType::FivePartialNested => write!(f, "FN"),
            OrfType::ThreePartialNested => write!(f, "TN"),
            OrfType::CompleteNested => write!(f, "CN"),
        }
    }
}

/// Parses a FASTA file into a hashmap of ORF records and their sequences
///
/// # Arguments
///
/// * `f` - The path to the FASTA file
///
/// # Returns
///
/// A `Result` containing a `HashMap` of `OrfRecord` and their sequences
///
/// # Example
///
/// ```rust
/// let orf_records = parse_fa("path/to/file.fa").unwrap();
/// ```
pub fn parse_fa<F: AsRef<Path>>(
    f: F,
) -> Result<HashMap<OrfRecord, Vec<u8>>, Box<dyn std::error::Error>> {
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
        let header = from_utf8(&entry[..header_end]).unwrap().trim();
        let record = &entry[header_end + 1..];
        let seq = record
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r') // Remove newlines and carriage returns
            .cloned()
            .collect::<Vec<u8>>();

        acc.insert(OrfRecord::parse(header), seq);
        pos = end;
    }

    Ok(acc)
}

/// Creates a new file with the given extension if it doesn't exist
///
/// # Arguments
///
/// * `file` - The path to the file
/// * `extension` - The extension of the file
///
/// # Returns
///
/// An `Option` containing a `BufWriter` for the file
///
/// # Example
///
/// ```rust
/// let writer = create_file("path/to/file.txt", "txt").unwrap();
/// ```
pub fn create_file(file: &Path, extension: &str) -> Option<BufWriter<File>> {
    let file = file.with_extension(extension);
    if !file.exists() {
        let writer = BufWriter::new(
            File::create(&file).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
        );

        Some(writer)
    } else {
        // INFO: overwrite existing fil
        dbg!("WARN: overwriting existing file -> {:?}", &file);
        std::fs::remove_file(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot remove existing file -> {e}");
        });

        Some(BufWriter::new(File::create(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot create file -> {e}")
        })))
    }
}

/// Deduplicates ORF records from a FASTA file
///
/// # Arguments
///
/// * `fasta` - The path to the FASTA file
/// * `do_nesting` - Whether to do nesting or not
/// * `min_len` - The minimum length of an ORF
/// * `min_percent` - The minimum percentage of an ORF
/// * `pattern` - The pattern to match for ORFs
///
/// # Returns
///
/// A `HashMap` of `HashHead` and their corresponding `Vec<OrfRecord>`
///
/// # Example
///
/// ```rust
/// let mapper = deduplicate("path/to/file.fa", true, 50, 0.25, &[b'M']).unwrap();
/// ```
pub fn deduplicate(
    fasta: &PathBuf,
    do_nesting: bool,
    min_len: usize,
    min_percent: f32,
    pattern: &[u8],
) -> HashMap<HashHead, Vec<OrfRecord>> {
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let mut mapper: HashMap<HashHead, Vec<OrfRecord>> = HashMap::new(); // INFO: seq -> records
    let mut idx = 1; // INFO: stores deduplicated index of sequence -> 1-based because 0 is placeholder

    // INFO: loops through sequences and populates mapper
    for (mut header, seq) in seqs.into_iter() {
        // INFO: sequence used as key in mapper -> hhead
        let len = seq.len();
        let mut key = Vec::with_capacity(len);
        for b in &seq {
            if *b != b'\n' {
                key.push(*b);
            }
        }

        let hhead = HashHead::new(Arc::from(key));

        // INFO: if hhead exists -> header inherits seq_idx
        match mapper.entry(hhead) {
            Entry::Occupied(mut entry) => {
                // INFO: seq exists -> push header to vec
                // INFO: seq_idx is inherited from first header
                entry.get_mut().push(header.clone());
            }
            Entry::Vacant(entry) => {
                header.seq_idx = idx;
                entry.insert(vec![header.clone()]);
                idx += 1; // INFO: increment index count
            }
        }

        if do_nesting {
            __split_record(
                &header,
                &seq,
                len,
                min_len,
                min_percent,
                &mut mapper,
                pattern,
                &mut idx,
            )
        }
    }

    mapper
}

/// Splits an ORF record into nested ORF records
///
/// # Arguments
///
/// * `header` - The header of the ORF record
/// * `seq` - The sequence of the ORF record
/// * `seq_length` - The length of the sequence
/// * `min_len` - The minimum length of the ORF
/// * `min_percent` - The minimum percentage of the ORF
/// * `mapper` - The mapper to store the nested ORF records
/// * `needle` - The needle to match for the nested ORF
/// * `idx_count` - The index count to increment
///
/// # Returns
///
/// None
///
/// # Example
///
/// ```rust
/// let mapper = HashMap::new();
/// let needle = b"M";
/// let idx_count = 0;
/// __split_record(&header, &seq, seq_length, min_len, min_percent, &mut mapper, &needle, &mut idx_count);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn __split_record(
    header: &OrfRecord,
    seq: &[u8],
    seq_length: usize,
    min_len: usize,
    min_percent: f32,
    mapper: &mut HashMap<HashHead, Vec<OrfRecord>>,
    needle: &[u8], // b"ATG" or b"M"
    idx_count: &mut usize,
) {
    let mut orf_count = 0;

    for (pos, &aa) in seq.iter().enumerate().skip(1) {
        // INFO: default needle is 'M' -> start signal
        if aa == needle[0] {
            let len_remaining = seq_length - pos;
            let percent = len_remaining as f32 / seq_length as f32;

            // INFO: we have found a nested ORF -> create a new record
            // WARN: updating start codons only in non-ATG cases because we are only looking for inner ATGs
            if len_remaining >= min_len && percent >= min_percent {
                orf_count += 1;

                let sub_seq = &seq[pos..];
                let mut nested_header = header.clone();

                nested_header.start += pos * 3;
                nested_header.id += format!("@{}", orf_count).as_str();

                // INFO: only updating start codon if it is not ATG
                if header.start_codon != "ATG" {
                    nested_header.start_codon = "ATG".to_string();
                }

                nested_header.orf_type = match header.orf_type {
                    OrfType::Complete => OrfType::CompleteNested,
                    OrfType::FivePartial => OrfType::FivePartialNested,
                    OrfType::ThreePartial => OrfType::ThreePartialNested,
                    _ => unreachable!(),
                };

                // INFO: subsequence is now a new key in the mapper
                let mut inner_key = Vec::with_capacity(sub_seq.len());
                for &b in sub_seq {
                    if b != b'\n' {
                        inner_key.push(b);
                    }
                }

                let hhead = HashHead::new(inner_key.into());
                match mapper.entry(hhead) {
                    Entry::Occupied(mut entry) => {
                        // INFO: seq exists -> push header to vec
                        // INFO: seq_idx is inherited from first header
                        entry.get_mut().push(nested_header.clone());
                    }
                    Entry::Vacant(entry) => {
                        nested_header.seq_idx = *idx_count;
                        entry.insert(vec![nested_header.clone()]);
                        *idx_count += 1; // INFO: increment index count
                    }
                }
            }
        }
    }
}

/// Gets the table of ORF records and their corresponding lines
///
/// # Arguments
///
/// * `mapper` - The mapper of ORF records and their sequences
/// * `bed` - The bed file
/// * `outdir` - The output directory
/// * `tai` - The path to the translationAi output
///
/// # Returns
///
/// A `HashMap` of `usize` and their corresponding `Vec<String>`
///
/// # Example
///
/// ```rust
/// let table = get_table(mapper, bed, outdir, tai).unwrap();
/// ```
#[allow(clippy::too_many_arguments)]
pub fn get_table(
    mut mapper: HashMap<HashHead, Vec<OrfRecord>>,
    mut bed: HashMap<String, GenePred>,
    outdir: &Path,
    tai: Option<PathBuf>,
    net: Option<PathBuf>,
    nmd_distance: u64,
    weak_nmd_distance: i64,
    atg_distance: u64,
    big_exon_dist_to_ej: u64,
    prefix: String,
) -> HashMap<usize, Vec<String>> {
    // INFO: sequence index -> vec of constructed lines
    let mut table = HashMap::new();

    // INFO: read tai values if provided
    let mut tai = if let Some(tai) = tai {
        info!("INFO: reading tai file: {tai:?}");
        read_tai_table(tai)
    } else {
        warn!("WARN: no tai file provided");
        HashMap::new()
    };

    // INFO: read tai values if provided
    let mut nets = if let Some(net) = net {
        info!("INFO: reading net file: {net:?}");
        read_nets(net)
    } else {
        warn!("WARN: no net file provided");
        HashMap::new()
    };

    let filename = &outdir.join(prefix);
    let mut pep = create_file(filename, "orf.dedup.pep")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", filename));

    let mut unique_tai_idx = 0;

    // INFO: iterate over mapper -> write index + mapper to records and then tai if provided
    for (head, records) in mapper.iter_mut() {
        let idx = records[0].seq_idx;

        // INFO: table holds seq idx -> vec of ORF records
        table.entry(idx).or_insert(Vec::new());

        for record in records.iter_mut() {
            unique_tai_idx += 1;

            let record_cannonical_id = record.id.split("_ORF").next().unwrap_or_else(|| {
                panic!(
                    "ERROR: failed to parse cannonical ID from header: {}",
                    record.id
                );
            });

            let gp = bed.get_mut(record_cannonical_id).unwrap_or_else(|| {
                panic!(
                    "ERROR: could not find gene prediction for ID: {}",
                    record_cannonical_id
                );
            });

            let (orf_start, orf_end) = map_absolute_cds(gp, record.start as u64, record.end as u64);
            let record_key = format!("{}:{}-{}", record_cannonical_id, orf_start, orf_end);

            // INFO: only important thing to inherit from orfipy is orf_type; strand is inherited from bed!
            if tai.contains_key(&record_key) {
                // INFO: ORF predicted by translationAi -> merge with orfipy and remove from tai
                let tai_prediction = tai.get(&record_key).unwrap_or_else(|| {
                    panic!(
                        "ERROR: could not find tai prediction for key: {}. This is a bug!",
                        record_key
                    );
                });

                if tai_prediction.relative_orf_start != record.start
                    || tai_prediction.relative_orf_end != record.end
                    || tai_prediction.start as u64 != orf_start
                    || tai_prediction.end as u64 != orf_end
                    || (record.start_codon != "NA"
                        && tai_prediction.start_codon != record.start_codon)
                    || (record.stop_codon != "NA" && tai_prediction.stop_codon != record.stop_codon)
                {
                    panic!(
                        "ERROR: tai prediction does not match orfipy prediction for key: {}. This is a bug -> {:?} vs {:?}! Key ({:?}) unlocked the following gp: {:?}",
                        record_key, tai_prediction, record, record_key, gp
                    );
                }

                let nmd_type = detect_nmd(
                    gp,
                    orf_start,
                    orf_end,
                    nmd_distance,
                    weak_nmd_distance,
                    atg_distance,
                    big_exon_dist_to_ej,
                );

                let mut line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    from_utf8(&gp.chrom).unwrap(),
                    orf_start,
                    orf_end,
                    record.id,
                    gp.strand.unwrap(),
                    tai_prediction.start_score,
                    tai_prediction.stop_score,
                    record.start,
                    record.end,
                    record.start_codon,
                    tai_prediction.stop_codon, // INFO: on purpose to avoid orfipy NA stop
                    tai_prediction.inner_stops,
                    record.orf_type,
                    nmd_type
                );

                if nets.contains_key(&record_key) {
                    debug!("DEBUG: net record for key: {record_key}");
                    let additional = format!(
                        "\t{}\t{}\t{}\t{}",
                        nets[&record_key].netstart_atg_score,
                        nets[&record_key].transaid_start_score,
                        nets[&record_key].transaid_stop_score,
                        nets[&record_key].transaid_integrated_score,
                    );

                    line.push_str(&additional);
                    debug!("DEBUG: line for {record_key} -> {line}");
                } else {
                    warn!("WARN: no net record for key: {record_key}");
                    line.push_str("\t-1\t-1\t-1\t-1");
                }

                table.entry(idx).or_insert(Vec::new()).push(line);

                // INFO: remove from tai -> tai retains tai-only
                // INFO: remove from net -> net retains net-only
                tai.remove(&record_key);
                nets.remove(&record_key);
            } else {
                // INFO: ORF not predicted by translationAi -> only skipped if it is not complete
                if record.orf_type != OrfType::Complete
                    && record.orf_type != OrfType::CompleteNested
                {
                    continue;
                };

                let nmd_type = detect_nmd(
                    gp,
                    orf_start,
                    orf_end,
                    nmd_distance,
                    weak_nmd_distance,
                    atg_distance,
                    big_exon_dist_to_ej,
                );

                // INFO: orfipy complete or complete-nested -> follow tai fmt + orf_type + strand
                let mut line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    from_utf8(&gp.chrom).unwrap(),
                    orf_start,
                    orf_end,
                    record.id,
                    gp.strand.unwrap(),
                    0.0, // INFO: no tai start score
                    0.0, // INFO: no tai stop score
                    record.start,
                    record.end,
                    record.start_codon,
                    record.stop_codon,
                    1, // INFO: complete or complete-nested only have 1 stop codon, otherwise would not be complete
                    record.orf_type,
                    nmd_type
                );

                if nets.contains_key(&record_key) {
                    debug!("DEBUG: net record for key: {record_key}");
                    let additional = format!(
                        "\t{}\t{}\t{}\t{}",
                        nets[&record_key].netstart_atg_score,
                        nets[&record_key].transaid_start_score,
                        nets[&record_key].transaid_stop_score,
                        nets[&record_key].transaid_integrated_score,
                    );
                    line.push_str(&additional);
                } else {
                    warn!("WARN: no net record for key: {record_key}");
                    line.push_str("\t-1\t-1\t-1\t-1");
                }

                table.entry(idx).or_insert(Vec::new()).push(line);
                nets.remove(&record_key);
            }
        }

        let handle = table.get(&idx).unwrap_or_else(|| {
            panic!("ERROR: could not find table for head -> {:?}", head);
        });

        if handle.is_empty() {
            println!("WARN: empty table for head -> {:?}", idx);
            table.remove(&idx);
        } else {
            let hh = format!(
                ">{}\n{}",
                idx,
                from_utf8(&head.seq).unwrap_or_else(|_| {
                    panic!("ERROR: failed to convert key to UTF-8 -> {:?}", head.seq);
                })
            );
            writeln!(pep, "{}", hh).unwrap_or_else(|e| {
                panic!("ERROR: failed to write record to file -> {e} -> {:?}", hh);
            });
        }
    }

    if !tai.is_empty() {
        for (_, prediction) in tai.iter() {
            unique_tai_idx += 1;

            let rc = format!(">{}\n{}", unique_tai_idx, prediction.aa);
            writeln!(pep, "{}", rc).unwrap_or_else(|e| {
                panic!("ERROR: failed to write record to file -> {e} -> {:?}", rc);
            });

            let record_cannonical_id = prediction.id.split(".p").next().unwrap_or_else(|| {
                panic!(
                    "ERROR: failed to parse cannonical ID from header: {}",
                    prediction.id
                );
            });

            let gp = bed.get_mut(record_cannonical_id).unwrap_or_else(|| {
                panic!(
                    "ERROR: could not find gene prediction for ID: {}",
                    record_cannonical_id
                );
            });
            let record_key = format!(
                "{}:{}-{}",
                record_cannonical_id, prediction.start, prediction.end
            );

            let nmd_type = detect_nmd(
                gp,
                prediction.start as u64,
                prediction.end as u64,
                nmd_distance,
                weak_nmd_distance,
                atg_distance,
                big_exon_dist_to_ej,
            );

            let mut line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                prediction.chr,
                prediction.start,
                prediction.end,
                prediction.id,
                prediction.strand,
                prediction.start_score,
                prediction.stop_score,
                prediction.relative_orf_start,
                prediction.relative_orf_end,
                prediction.start_codon,
                prediction.stop_codon,
                prediction.inner_stops,
                "UN", // INFO: unknown orf_type -> can be guessed with codons but its helpful for tai-only
                nmd_type,
            );

            if nets.contains_key(&record_key) {
                debug!("DEBUG: net record for key in a tai-only table: {record_key}");
                let additional = format!(
                    "\t{}\t{}\t{}\t{}",
                    nets[&record_key].netstart_atg_score,
                    nets[&record_key].transaid_start_score,
                    nets[&record_key].transaid_stop_score,
                    nets[&record_key].transaid_integrated_score,
                );

                line.push_str(&additional);
            } else {
                warn!("WARN: no net record for key in a tai-only table: {record_key}");
                line.push_str("\t-1\t-1\t-1\t-1");
            }

            table.entry(unique_tai_idx).or_insert(Vec::new()).push(line);
            nets.remove(&record_key);
        }
    }

    if !nets.is_empty() {
        info!(
            "INFO: there are {} net records remaining in nets",
            nets.len()
        );
        warn!("WARN: net records remaining in nets: {nets:?}")
    }

    table
}

/// An enum representing the NMD type of an ORF
pub enum NMDType {
    NN, // INFO: no_nmd
    WN, // INFO: weak_nmd
    SN, // INFO: strong_nmd
    UN, // INFO: unknown -> likely an error
}

impl std::fmt::Display for NMDType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NMDType::NN => write!(f, "NN"),
            NMDType::WN => write!(f, "WN"),
            NMDType::SN => write!(f, "SN"),
            NMDType::UN => write!(f, "UN"),
        }
    }
}

/// Detects the NMD type of an ORF
/// NMD types are:
/// 1. NN -> no_nmd
/// 2. WN -> weak_nmd
/// 3. SN -> strong_nmd
///
/// # Arguments
/// * `gp` - The genepred record
/// * `orf_start` - The start of the ORF
/// * `orf_end` - The end of the ORFs
/// * `nmd_distance` - The minimum 3' UTR length for a strong NMD classification.
/// * `weak_nmd_distance` - The maximum distance from the stop codon to the last exon-exon junction
///   for a weak NMD classification.
/// * `atg_distance` - The maximum CDS length for a weak NMD classification.
/// * `big_exon_dist_to_ej` - The maximum distance from the stop codon to the next splice junction
///   for a big exon test, used in weak NMD classification.
///
/// # Returns
/// A `String` representing the NMD type
///
/// # Example
/// ```rust, ignore
/// let nmd_type = detect_nmd(gp, orf_start, orf_end);
/// ```
fn detect_nmd(
    gp: &GenePred,
    orf_start: u64,
    orf_end: u64,
    nmd_distance: u64,
    weak_nmd_distance: i64,
    atg_distance: u64,
    big_exon_dist_to_ej: u64,
) -> NMDType {
    let mut cds_start = orf_start;
    let mut cds_end = orf_end;

    // INFO: noncoding transcripts + errors
    if cds_start == cds_end {
        return NMDType::UN;
    }

    let mut nmd_count: i64 = -1;
    let mut _ex_ex_junction_utr: i64 = -1;
    let mut dist_stop_to_next_sj = 0; // INFO: for big exon test
    let mut in_utr = false;
    let mut utr_len = 0;
    let mut cds_len = 0;
    let mut bp_utr_to_last_ex_ex_jct = 0;

    let mut exons = gp.exons();
    exons = match gp.strand() {
        Some(genepred::Strand::Forward) => exons,
        Some(genepred::Strand::Reverse) => {
            // INFO: we scale the coordinates and then sort them
            // INFO: this is equivalent to reversing the exons and reversing the coordinates
            let mut exons = exons
                .iter()
                .map(|(start, end)| (SCALE - end, SCALE - start))
                .collect::<Vec<(u64, u64)>>();

            exons.sort_by(|a, b| a.0.cmp(&b.0));

            cds_end = SCALE - orf_start;
            cds_start = SCALE - orf_end;

            exons
        }
        None | Some(genepred::Strand::Unknown) => {
            panic!(
                "ERROR: failed to parse strand from genepred record -> {:?}",
                gp
            )
        }
    };

    for (i, exon) in exons.iter().enumerate() {
        let exon_start = exon.0;
        let exon_end = exon.1;

        // INFO: Count EEJs in 3'UTR
        if exon_end >= cds_end {
            _ex_ex_junction_utr += 1;

            // INFO: first exon containing stop codon
            if dist_stop_to_next_sj == 0 {
                dist_stop_to_next_sj = exon_end - cds_end;
            }

            if !in_utr {
                utr_len += exon_end - cds_end;
                in_utr = true;
            } else {
                utr_len += exon_end - exon_start;
            }

            if utr_len >= nmd_distance {
                nmd_count += 1;
            }

            // If last exon, compute bpUTRtoLastEEJ
            if i == exons.len() - 1 {
                bp_utr_to_last_ex_ex_jct = utr_len as i64 - (exon_end as i64 - exon_start as i64);
            }
        }

        // INFO: CDS length accumulation
        if exon_end < cds_start || exon_start > cds_end {
            continue; // INFO: skip pure UTR exons
        }

        // INFO: first coding exon
        if exon_end >= cds_start && cds_start >= exon_start {
            if exon_end >= cds_end {
                cds_len += cds_end - cds_start;
            } else {
                cds_len += exon_end - cds_start;
            }
        }
        // INFO: internal coding exon
        else if exon_start > cds_start && exon_end < cds_end {
            cds_len += exon_end - exon_start;
        }
        // INFO: last coding exon
        else if exon_start > cds_start && exon_end >= cds_end {
            cds_len += cds_end - exon_start;
        }
    }

    // INFO: final classification -> tag [NN: no_nmd, SN: strong_nmd, WN: weak_nmd]
    if nmd_count == 0 || nmd_count == -1 {
        NMDType::NN
    } else if bp_utr_to_last_ex_ex_jct <= weak_nmd_distance
        || dist_stop_to_next_sj >= big_exon_dist_to_ej
        || cds_len <= atg_distance
    {
        NMDType::WN
    } else {
        NMDType::SN
    }
}

/// Reads the translationAi output and converts it to a hashmap of `TaiRecord`s
///
/// # Arguments
///
/// * `tai` - The path to the translationAi output
///
/// # Returns
///
/// A `HashMap` of `String` and their corresponding `TaiRecord`
///
/// # Example
///
/// ```rust
/// let tai = read_tai_table(tai).unwrap();
/// ```
fn read_tai_table(tai: PathBuf) -> HashMap<String, TaiRecord> {
    let predictions = reader(tai)
        .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

    let mut accumulator = HashMap::new();

    // INFO: chr17 26054102 26060105 ENSMUST00000176751.2.p3 - 0.011934959 0.9905367 242 2126 ATG TAG 8
    for line in predictions.lines() {
        let record = TaiRecord::parse(line);
        let cannonical_id =
            record.id.split(".p").next().unwrap_or_else(|| {
                panic!("ERRROR: failed to parse cannonical ID from line: {}", line)
            }); // INFO: ENSMUST00000176751.2.p3 -> ENSMUST00000176751.2
        let key = format!("{}:{}-{}", cannonical_id, record.start, record.end); // INFO: ENSMUST00000176751.2:26054102-26060105
        accumulator.insert(key, record);
    }

    accumulator
}

/// A struct representing a translationAi record
#[derive(Debug, Clone)]
pub struct TaiRecord {
    pub chr: String,
    pub start: usize,
    pub end: usize,
    pub id: String,
    pub strand: char,
    pub start_score: f32,
    pub stop_score: f32,
    pub relative_orf_start: usize,
    pub relative_orf_end: usize,
    pub start_codon: String,
    pub stop_codon: String,
    pub inner_stops: usize,
    pub aa: String,
}

impl TaiRecord {
    pub fn parse(line: &str) -> Self {
        let mut parts = line.split('\t');

        let (
            chr,
            start,
            end,
            id,
            strand,
            start_score,
            stop_score,
            relative_orf_start,
            relative_orf_end,
            start_codon,
            stop_codon,
            inner_stops,
            aa,
        ) = (
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai chr from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai start from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse tai start from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai end from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to parse tai end from line: {} -> {e}", line);
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai id from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai strand from line: {}", line);
                })
                .chars()
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai strand from line: {}", line);
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai start score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse tai start score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai stop score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse tai stop score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse relative orf start from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse relative orf start from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse relative orf end from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse relative orf end from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse start codon from line: {}", line);
                })
                .to_string()
                .to_uppercase(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse stop codon from line: {}", line);
                })
                .to_string()
                .to_uppercase(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse inner stops from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse inner stops from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse aa from line: {}", line);
                })
                .to_string(),
        );

        Self {
            chr,
            start,
            end,
            id,
            strand,
            start_score,
            stop_score,
            relative_orf_start,
            relative_orf_end,
            start_codon,
            stop_codon,
            inner_stops,
            aa,
        }
    }
}

/// Net representation of netstart and transaid
#[derive(Debug, Clone)]
#[allow(unused)]
struct NetRecord {
    sequence_id: String,
    orf_absolute_start: usize,
    orf_absolute_end: usize,
    orf_relative_start: usize,
    orf_relative_end: usize,
    strand: genepred::Strand,
    netstart_atg_score: f32,
    transaid_start_score: f32,
    transaid_stop_score: f32,
    transaid_integrated_score: f32,
}

impl NetRecord {
    /// Parse netstart line
    fn parse(line: &str) -> Self {
        let mut parts = line.split('\t');

        let (
            sequence_id,
            orf_absolute_start,
            orf_absolute_end,
            orf_relative_start,
            orf_relative_end,
            strand,
            netstart_atg_score,
            transaid_start_score,
            transaid_stop_score,
            transaid_integrated_score,
        ) = (
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse sequence_id from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse orf_absolute_start from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse orf_absolute_start from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse orf_absolute_end from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse orf_absolute_end from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse orf_relative_start from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse orf_relative_start from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse orf_relative_end from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse orf_relative_end from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse strand from line: {}", line);
                })
                .chars()
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse strand from line: {}", line);
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse netstart atg score from line: {}",
                        line
                    );
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse netstart atg score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse transaid start score from line: {}",
                        line
                    );
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse transaid start score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse transaid stop score from line: {}",
                        line
                    );
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse transaid stop score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse transaid integrated score from line: {}",
                        line
                    );
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse transaid integrated score from line: {} -> {e}",
                        line
                    );
                }),
        );

        let strand = match strand {
            '+' => genepred::Strand::Forward,
            '-' => genepred::Strand::Reverse,
            _ => panic!("ERROR: failed to parse strand from line: {}", line),
        };

        Self {
            sequence_id,
            orf_absolute_start,
            orf_absolute_end,
            orf_relative_start,
            orf_relative_end,
            strand,
            netstart_atg_score,
            transaid_start_score,
            transaid_stop_score,
            transaid_integrated_score,
        }
    }
}

/// Reads the net output and converts it to a hashmap of `NetRecord`s
///
/// # Arguments
///
/// * `path` - The path to the net output
///
/// # Returns
///
/// A `HashMap` of `String` and their corresponding `NetRecord`
///
/// # Example
///
/// ```rust
/// let nets = read_nets(net).unwrap();
/// ```
fn read_nets<P: AsRef<Path> + std::fmt::Debug>(path: P) -> HashMap<String, NetRecord> {
    let contents =
        reader(&path).unwrap_or_else(|e| panic!("ERROR: failed to read net file {path:?} -> {e}"));
    let mut accumulator = HashMap::new();

    contents.lines().for_each(|line| {
        let record = NetRecord::parse(line);
        let key = format!(
            "{}:{}-{}",
            record.sequence_id, record.orf_absolute_start, record.orf_absolute_end
        );

        accumulator.entry(key).or_insert(record);
    });

    accumulator
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_parts_reads_shared_diamond_mmseqs_layout() {
        let parts = [
            "17",
            "97.2",
            "142",
            "357",
            "141",
            "1",
            "141",
            "217",
            "357",
            "5.09e-93",
            "sp|O55165|KIF3C_RAT",
        ];
        let record = BlastRecord::from_parts(&parts);
        assert_eq!(record.blast_pid, 97.2);
        assert!((record.blast_e_value - 5.09e-93).abs() < 1e-100);
        assert_eq!(record.blast_offset, 216);
        assert_eq!(record.blast_alignment_len, 141);
        assert!((record.percent_aligned - (141.0 / 142.0 * 100.0)).abs() < 0.01);
        assert_eq!(record.reference_id, "KIF3C_RAT");
    }

    #[test]
    fn esmfold2_input_and_score_protocol() {
        for (input, expected) in [
            ("ABC*", "ABC"),
            ("AB*C", "ABXC"),
            ("AB*C*", "ABXC"),
            ("A**BC*", "AXXBC"),
        ] {
            assert_eq!(sanitize_esmfold2_sequence(input), expected);
        }

        let boundary = format!(">7\n{}*\n>8\n{}\n", "A".repeat(2048), "A".repeat(2049));
        let records = parse_numeric_fasta(&boundary).unwrap();
        assert_eq!(sanitize_esmfold2_sequence(&records[0].1).len(), 2048);
        assert_eq!(sanitize_esmfold2_sequence(&records[1].1).len(), 2049);

        let scores = parse_esmfold2_scores("7\t81.2500\n8\t0.0000\n").unwrap();
        assert_eq!(scores[&7], 81.25);
        assert_eq!(scores[&8], 0.0);
        assert!(parse_esmfold2_scores("7\tNaN\n").is_err());
        assert!(parse_esmfold2_scores("7\t50\n7\t51\n").is_err());
    }
}
