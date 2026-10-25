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

use dashmap::DashSet;
use genepred::{Bed12, GenePred, Strand};
use hashbrown::HashMap;
use smol_str::{SmolStr, ToSmolStr};

use std::fs::File;
use std::io::{BufWriter, Write};
use std::str::from_utf8;

use crate::{cli::TaiArgs, consts::*, utils::*};

/// Runs the translationAi pipeline.
///
/// # Arguments
///
/// * `args` - The command-line arguments.
///
/// # Example
///
/// ```rust
/// use clap::Parser;
/// use std::path::PathBuf;
///
/// #[derive(Debug, Parser)]
/// #[command(name = "tai", about = "In-house translationAi caller", version = env!("CARGO_PKG_VERSION"))]
/// pub struct Args {
///     #[arg(
///         short = 'f',
///         long = "fasta",
///         required = true,
///         help = "Path to .fa file"
///     )]
///     pub fasta: PathBuf,
///
///     #[arg(short = 'b', long = "bed", required = true, help = "Path to .bed file")]
///     pub bed: PathBuf,
///
///     #[arg(
///         short = 'o',
///         long = "outdir",
///         required = false,
///         help = "Path to outdir",
///         default_value = "."
///     )]
///     pub outdir: PathBuf,
/// }
///
/// fn run() {
///     let args = Args::parse();
///     run_tai(args);
/// }
/// ```
pub fn run_tai(args: TaiArgs) {
    let dir = args.outdir.join("tai");
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let records = genepred::Reader::<Bed12>::from_mmap(args.bed)
        .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
        .filter_map(|record| {
            record.ok().map(|record| {
                (
                    from_utf8(record.name().unwrap()).unwrap().to_string(),
                    record,
                )
            })
        })
        .collect::<HashMap<String, genepred::GenePred>>();

    let (fasta, sequences) = refmt(&args.fasta, &records, &dir);

    let cmd = format!(
        "translationai -I {} -t 0.01,0.01 -O {}",
        fasta.display(),
        fasta.with_extension(PREDICTIONS).display()
    );

    run_command(&cmd, "translationai");

    let predictions = reader(fasta.with_extension(PREDICTIONS))
        .unwrap_or_else(|e| panic!("ERROR: failed to read predictions file -> {e}"));

    let writer = BufWriter::new(
        File::create(fasta.with_extension("result"))
            .unwrap_or_else(|e| panic!("ERROR: cannot create index from sequences -> {e}")),
    );

    tai(
        predictions,
        records,
        writer,
        sequences,
        args.upstream_flank,
        args.downstream_flank,
    );
}

/// Runs the translationAi pipeline.
///
/// # Arguments
///
/// * `predictions` - The path to the predictions file.
/// * `records` - A hash map of gene predictions.
/// * `accumulator` - A set of unique IDs.
/// * `writer` - A writer for the output file.
/// * `sequences` - A hash map of sequences.
///
/// # Example
///
/// ```rust
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::BufWriter;
///
/// let predictions = PathBuf::from("path/to/predictions");
/// let records = HashMap::new();
/// let accumulator = DashSet::new();
/// let writer = BufWriter::new(File::create("path/to/output").unwrap());
/// let sequences = HashMap::new();
///
/// tai(predictions, records, accumulator, writer, sequences);
/// ```
#[allow(unused_assignments)]
fn tai(
    predictions: String,
    mut records: HashMap<String, GenePred>,
    mut writer: BufWriter<File>,
    sequences: HashMap<SmolStr, Vec<u8>>,
    upstream_flank: usize,
    downstream_flank: usize,
) {
    let accumulator = DashSet::new();

    for line in predictions.lines() {
        let parts: Vec<&str> = line.split('\t').collect();

        // INFO: >chr6:128352418-128362107(+)(ENSMUST00000100926.4)(0, 0,)
        let mut id = parts
            .first()
            .unwrap_or_else(|| {
                panic!("ERROR: ID not found in line: {}", line);
            })
            .to_string();

        // INFO: ENSMUST00000100926.4
        id = id.split("(").collect::<Vec<&str>>()[2]
            .split(")")
            .collect::<Vec<&str>>()[0]
            .to_string();

        for (orf_idx, orf) in parts.iter().skip(1).enumerate() {
            // INFO: 13190,13667,0.6413461491465569,0.921319767832756
            let orf_parts: Vec<&str> = orf.split(',').collect();
            if orf_parts.len() < 2 {
                panic!("ERROR: ORF does not have enough parts to parse: {}", orf);
            }

            let start = orf_parts[0].parse::<u64>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse start position from ORF: {}", orf);
            });
            let mut stop = orf_parts[1].parse::<u64>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse stop position from ORF: {}", orf);
            }); // INFO: stop is inclusive, so we add 3 to include the stop codon
            let start_score = orf_parts[2].parse::<f32>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse start position from ORF: {}", orf);
            });
            let stop_score = orf_parts[3].parse::<f32>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse stop position from ORF: {}", orf);
            });

            if start > stop {
                panic!(
                    "ERROR: start position is greater than stop position in ORF: {}",
                    orf
                );
            }

            // INFO: since context is included, any predicted ORF start < upstream_flank
            // is considered unreliable because has been predicted in the upstream context
            if start < upstream_flank as u64 || stop < upstream_flank as u64 {
                println!("WARN: skipping unreliable ORF predicted in context flanks -> {orf:?}");
                continue;
            }

            // INFO: retrieving the reference gene prediction record
            let record = records.get_mut(&id).unwrap_or_else(|| {
                panic!("ERROR: id not found in BED, this is a bug -> {}!", id);
            });

            if stop - downstream_flank as u64 > record.exonic_length() {
                println!(
                    "WARN: predicted ORF is out-of-bounds -> {orf:?} for {record:?} with exonic_length -> {:?}",
                    record.exonic_length()
                );
                continue;
            };

            let strand = record.strand;

            let (mut orf_start, mut orf_end) = map_absolute_cds(
                record,
                start - upstream_flank as u64,
                stop - downstream_flank as u64,
            );

            // WARN: skipping unreliable ORFs for the current alignment
            if (orf_start == 0 && orf_end == 0) || orf_start > orf_end || orf_end - orf_start < 3 {
                println!(
                    "WARN: skipping unreliable ORF -> {orf:?} with mapped start -> {orf_start:?}, end -> {orf_end:?} for start -> {start:?}, stop -> {stop:?}"
                );
                continue;
            }

            let mut stop_codon = String::new();

            // INFO: before adding up we need to check the stop codon to see if its a real one or not
            let sequence = sequences.get(&id.to_smolstr()).unwrap_or_else(|| {
                panic!("ERROR: sequence not found for {id:?}");
            });
            let mut orf_sequence = sequence.clone();

            let start_codon =
                from_utf8(&sequence[start as usize..(start + 3) as usize])
                    .unwrap_or_else(|e| {
                        panic!("ERROR: failed to parse start codon -> {e} -> {sequence:?} from {id:?} using {orf_start:?}");
                    }).to_uppercase();

            __check_start_codon(&start_codon, &id, start);

            // INFO: stop is inclusive, so we add 3 to include the stop codon
            match record.strand() {
                // WARN: some weird cases where the tool predicts a non-stopped ORF:
                // WARN: R146001_manual_scaffold_1.p1    102     2199
                // WARN: sizes -> 144,117,332,112,120,117,138,78,66,270,246,137,79,156,87 = 2199
                // WARN: record -> manual_scaffold_1 189532046 189543938 R146001 60 + 189532148 189543941
                // WARN: would be out-of-bounds if we add +3 in this case
                Some(Strand::Forward) => {
                    if orf_end + 3 > record.end {
                        println!(
                            "WARN: translationAi predicted a non-stop ORF: {orf:?} for {record:?}"
                        );
                        // orf_end = record.end
                        stop_codon =
                            from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                .unwrap_or_else(|e| {
                                    panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                }).to_uppercase();
                        println!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);

                        orf_sequence = sequence[start as usize..stop as usize].to_vec();
                    } else {
                        stop_codon =
                            from_utf8(&sequence[(stop) as usize..(stop + 3) as usize])
                                .unwrap_or_else(|e| {
                                    panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                }).to_uppercase();

                        if !STOP_CODONS.contains(&stop_codon.as_str()) {
                            // WARN: if stop_codon is not cannonical this is probably a case where the tool is wrong
                            println!(
                                "WARN: stop codon is not TAA, TAG, or TGA -> {:?} from {:?} using {:?}",
                                stop_codon, &id, stop
                            );

                            // INFO: taking stop as last base and going back 2 nt to capture last codon
                            stop_codon =
                                from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                    }).to_uppercase();

                            println!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                            orf_sequence = sequence[start as usize..stop as usize].to_vec();
                        } else {
                            // INFO: stop_codon is cannonical so we can safely add 3 to the end
                            orf_end += 3;
                            orf_sequence = sequence[start as usize..(stop + 3) as usize].to_vec();
                            stop += 3;
                        }
                    }
                }
                Some(Strand::Reverse) => {
                    // INFO: not orf_start - 3, which underflows for a CDS starting within 3 bp of the scaffold start
                    if orf_start < record.start + 3 {
                        println!(
                            "WARN: translationAi predicted a non-stop ORF: {orf:?} for {record:?}"
                        );
                        // orf_start = config::SCALE - record.end
                        stop_codon =
                            from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                .unwrap_or_else(|e| {
                                    panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                }).to_uppercase();
                        println!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                        orf_sequence = sequence[start as usize..stop as usize].to_vec();
                    } else {
                        stop_codon =
                            from_utf8(&sequence[(stop) as usize..(stop + 3) as usize])
                                .unwrap_or_else(|e| {
                                    panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                }).to_uppercase();

                        if !STOP_CODONS.contains(&stop_codon.as_str()) {
                            // WARN: if stop_codon is not cannonical this is probably a case where the tool is wrong
                            println!(
                                "WARN: stop codon is not TAA, TAG, or TGA -> {:?} from {:?} using {:?}",
                                stop_codon, &id, stop
                            );

                            // INFO: taking stop as last base and going back 2 nt to capture last codon
                            stop_codon =
                                from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                    }).to_uppercase();

                            println!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                            orf_sequence = sequence[start as usize..stop as usize].to_vec();
                        } else {
                            // INFO: stop_codon is cannonical so we can safely add 3 to the end
                            orf_start -= 3;
                            orf_sequence = sequence[start as usize..(stop + 3) as usize].to_vec();
                            stop += 3;
                        }
                    }
                }
                _ => panic!("ERROR: unexpected strand value: {:?}", strand),
            }

            let pep = translate(&orf_sequence);
            let inner_stops = scan_stops(orf_sequence);

            // INFO: retrieving the reference gene prediction record
            // INFO: since indexing groups exact similar records
            // INFO: we safely assume ref gp record could be applied to all queries
            let ref_id = format!("{}.p{}", id, orf_idx + 1);

            let ref_line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                from_utf8(&record.chrom)
                    .unwrap_or_else(|e| panic!("ERROR: failed to parse chrom -> {e}")),
                orf_start,
                orf_end,
                ref_id,
                strand.unwrap_or_else(|| panic!("ERROR: strand not found for {id:?}")),
                start_score,
                stop_score,
                start - upstream_flank as u64,
                stop - downstream_flank as u64,
                start_codon,
                stop_codon,
                inner_stops,
                pep
            );

            accumulator.insert(ref_line);
        }
    }

    accumulator.into_iter().for_each(|line| {
        writeln!(writer, "{}", line).unwrap();
    });
}
