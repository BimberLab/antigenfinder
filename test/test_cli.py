import unittest

from antigenfinder.__main__ import build_parser

class TestCommandLine(unittest.TestCase):
    def test_skip_processing_defaults_to_false(self):
        args = build_parser().parse_args([
            "prepare_gtf",
            "--gtf-file", "input.gtf",
            "--fasta-file", "genome.fa",
            "--cache-dir", "cache",
        ])

        self.assertFalse(args.skip_processing_if_cached)

    def test_skip_processing_flag(self):
        args = build_parser().parse_args([
            "prepare_gtf",
            "--gtf-file", "input.gtf",
            "--fasta-file", "genome.fa",
            "--cache-dir", "cache",
            "--skip-processing-if-cached",
        ])

        self.assertTrue(args.skip_processing_if_cached)

    def test_required_argument(self):
        with self.assertRaises(SystemExit) as error:
            build_parser().parse_args(["prepare_gtf"])

        self.assertEqual(error.exception.code, 2)