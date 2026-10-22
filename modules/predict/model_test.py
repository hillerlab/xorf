#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import predict

DATA = Path(__file__).resolve().parent / "fixtures"
SCRIPT = Path(__file__).resolve().parent / "predict.py"


class BundledModelTests(unittest.TestCase):
    """The bundled models load and score a real chunk as the training code does.

    fixtures/ is one PREDICT chunk (two HLratNor8 coding transcripts and one
    lncRNA, every candidate of each) with the probabilities the training code
    computed for both bundled models.
    """

    def test_bundles_load(self) -> None:
        for path, n_features in ((predict.MODEL_PATH, 63), (predict.MODEL_NO_TRANSAID_PATH, 59)):
            model, bundle = predict.load_model(path)
            self.assertEqual(len(bundle["features"]), n_features)
            self.assertEqual(model.n_features_in_, n_features)
            self.assertTrue(0 < bundle["threshold"] < 1)

    def test_no_transaid_model_has_no_transaid_features(self) -> None:
        _, bundle = predict.load_model(predict.MODEL_NO_TRANSAID_PATH)
        self.assertEqual([f for f in bundle["features"] if "transaid" in f], [])

    def test_scores_match_training_code(self) -> None:
        expected = pd.read_csv(DATA / "expected.tsv", sep="\t")
        for path, column in (
            (predict.MODEL_PATH, "prob_coding"),
            (predict.MODEL_NO_TRANSAID_PATH, "prob_coding_no_transaid"),
        ):
            model, bundle = predict.load_model(path)
            table = predict.predict(DATA / "chunk.blast.tsv", DATA / "chunk.samba.tsv", model, bundle)
            scores = table.set_index("id")["prob_coding"].reindex(expected["id"]).to_numpy()
            np.testing.assert_allclose(scores, expected[column].to_numpy(), rtol=0, atol=1e-5)

    def test_cli_keeps_calls_above_the_model_cutoff(self) -> None:
        expected = pd.read_csv(DATA / "expected.tsv", sep="\t")
        for flags, column, path in (
            ([], "prob_coding", predict.MODEL_PATH),
            (["--no-transaid"], "prob_coding_no_transaid", predict.MODEL_NO_TRANSAID_PATH),
        ):
            cutoff = predict.load_model(path)[1]["threshold"]
            with tempfile.TemporaryDirectory() as directory:
                run = subprocess.run(
                    [sys.executable, str(SCRIPT), "-b", str(DATA / "chunk.blast.tsv"),
                     "-s", str(DATA / "chunk.samba.tsv"), "-a", str(DATA / "chunk.bed"),
                     "-o", directory, "-p", "chunk", "-mp", "1", *flags],
                    check=True, capture_output=True, text=True,
                )
                kept = pd.read_csv(Path(directory) / "chunk.predictions.tsv", sep="\t")
            # the fixture's best scores cannot tell the cutoff from 0.50, the log line can
            self.assertIn(f"min score to keep={cutoff:.4f}", run.stderr)
            above = expected[expected[column] >= cutoff]
            best = above.loc[above.groupby(above["id"].str.split("_ORF").str[0])[column].idxmax(), "id"]
            self.assertEqual(sorted(kept["id"]), sorted(best))

    def test_model_and_no_transaid_are_exclusive(self) -> None:
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "-b", "x", "-s", "x", "-a", "x",
             "--model", str(predict.MODEL_PATH), "--no-transaid"],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)


if __name__ == "__main__":
    unittest.main()
