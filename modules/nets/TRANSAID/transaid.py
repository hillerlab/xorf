#!/usr/bin/env python
# predict.py

import argparse
import sys
from pathlib import Path
import logging
from src.predictor import TranslationPredictor
from src.utils import save_results, setup_logger, find_perfect_orfs, cleanup_temp_file

MODEL = "model/TRANSAID_Embedding_batch4_best_model.pth"
MODEL_PATH = Path(__file__).resolve().parent / MODEL


__version__ = "0.0.2-orf"


def get_args():
    parser = argparse.ArgumentParser(
        description="Predict translation products from RNA sequences"
    )
    parser.add_argument("--input", type=str, required=True, help="Input FASTA file")
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Base path for output files (CSV and FASTA)",
    )
    parser.add_argument(
        "--model_path",
        type=str,
        required=False,
        help="Path to pretrained model",
        default=MODEL_PATH,
    )
    parser.add_argument("--gpu", type=int, default=-1, help="GPU device ID, -1 for CPU")
    parser.add_argument(
        "--integrated_cutoff",
        type=float,
        default=0.65,
        help="Cutoff value for integrated score",
    )
    parser.add_argument(
        "--filter_mode",
        choices=["all", "best"],
        default="best",
        help="Keep all ORFs above cutoff or only the best one",
    )
    parser.add_argument(
        "--save_raw_predictions",
        action="store_true",
        help="Whether to save raw prediction results in PKL format",
    )
    parser.add_argument(
        "--batch_size", type=int, default=32, help="Batch size for prediction"
    )
    parser.add_argument(
        "--max_seq_len",
        type=int,
        default=27109,
        help="Maximum sequence length to process",
    )
    parser.add_argument(
        "--pad_margin",
        type=int,
        default=1024,
        help="Pad each length-sorted batch to its longest sequence plus this margin "
        "(rounded up to 256, capped at --max_seq_len); >= 1009 leaves scores unchanged; "
        "negative pads every sequence to --max_seq_len",
    )
    parser.add_argument(
        "--max_tokens",
        type=int,
        default=0,
        help="Cap batch size x padded length per batch (0 = no cap)",
    )
    parser.add_argument(
        "--orf_length_cutoff",
        type=int,
        default=30,
        help="Minimum ORF length in amino acids",
    )
    parser.add_argument(
        "--kozak_cutoff", type=float, default=0.05, help="Minimum Kozak sequence score"
    )
    parser.add_argument(
        "--tis_cutoff", type=float, default=0.1, help="Minimum TIS prediction score"
    )
    parser.add_argument(
        "--tts_cutoff", type=float, default=0.1, help="Minimum TTS prediction score"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Version: {__version__}",
    )
    return parser.parse_args()


def main():
    args = get_args()

    output_base = Path(args.output)
    output_dir = output_base.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    log_file = output_base.with_suffix(".log")
    if log_file.exists():
        with open(log_file, "w") as f:
            f.truncate(0)
    logger = setup_logger(name="translation_predictor", log_file=str(log_file))

    remaining_fasta_path_to_clean = None
    try:
        input_path = Path(args.input)
        if not (
            input_path.is_file()
            and input_path.suffix.lower() in {".fasta", ".fa", ".fna"}
        ):
            raise FileNotFoundError(
                f"Input file not found or is not a FASTA file: {args.input}"
            )

        # Step 1: Pre-screening for perfect ORFs
        perfect_orf_results, remaining_fasta_path = find_perfect_orfs(input_path)
        remaining_fasta_path_to_clean = remaining_fasta_path

        model_results = []
        probs = None
        if remaining_fasta_path:
            logger.info("Initializing predictor for remaining sequences...")
            predictor = TranslationPredictor(
                model_path=args.model_path,
                device=args.gpu,
                batch_size=args.batch_size,
                sequence_length=args.max_seq_len,
                tis_cutoff=args.tis_cutoff,
                tts_cutoff=args.tts_cutoff,
                orf_length_cutoff=args.orf_length_cutoff,
            )

            logger.info("Starting model prediction...")
            model_results, probs = predictor.predict_batch(
                remaining_fasta_path,
                pad_margin=args.pad_margin,
                max_tokens=args.max_tokens,
            )
        else:
            logger.info("No sequences remain for model prediction. Skipping.")

        all_results = perfect_orf_results + model_results

        logger.info("Applying filters and marking results...")
        for r in all_results:
            if r.passed_filter:
                continue  # Skip already passed results (e.g., perfect ORFs)

            reasons = []
            if r.integrated_score <= args.integrated_cutoff:
                reasons.append(
                    f"integrated_score({r.integrated_score:.3f}) <= {args.integrated_cutoff}"
                )
            if len(r.protein_sequence) < args.orf_length_cutoff:
                reasons.append(
                    f"protein_length({len(r.protein_sequence)}) < {args.orf_length_cutoff}"
                )

            if not reasons:
                r.passed_filter = True
            else:
                r.passed_filter = False
                r.filter_reason = "; ".join(reasons)

        if args.filter_mode == "best":
            transcript_groups = {}
            for r in all_results:
                transcript_groups.setdefault(r.sequence_id, []).append(r)

            for transcript_id, transcript_results in transcript_groups.items():
                passed_results = [r for r in transcript_results if r.passed_filter]
                if len(passed_results) > 1:
                    best_result = max(passed_results, key=lambda x: x.integrated_score)
                    for r in passed_results:
                        if r != best_result:
                            r.passed_filter = False
                            r.filter_reason = (
                                r.filter_reason + "; " if r.filter_reason else ""
                            ) + "not the best result"

        logger.info("Saving results...")
        sorted_results = sorted(
            all_results, key=lambda x: x.integrated_score, reverse=True
        )
        csv_path = output_base.with_suffix(".csv")
        faa_path = output_base.with_suffix(".faa")

        save_results(
            probs=probs,
            results=sorted_results,
            csv_path=csv_path,
            faa_path=faa_path,
            save_raw=args.save_raw_predictions,
            save_only_passed=True,
        )

        logger.info("Prediction completed successfully!")

    except Exception as e:
        logger.error(f"An unexpected error occurred: {str(e)}", exc_info=True)
        sys.exit(1)

    finally:
        if remaining_fasta_path_to_clean:
            logger.info(f"Cleaning up temporary file: {remaining_fasta_path_to_clean}")
            cleanup_temp_file(remaining_fasta_path_to_clean)


if __name__ == "__main__":
    main()
