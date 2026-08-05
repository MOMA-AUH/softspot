from __future__ import annotations

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


REPOSITORY = Path(__file__).resolve().parents[1]
_OUTPUT_HEADER = "fraction\tpositive_pairs\teligible_pairs\n"


def _expected_output(score: float, numerator: int, denominator: int) -> str:
    return f"{_OUTPUT_HEADER}{score:.6f}\t{numerator}\t{denominator}\n"


def _write_alignment(
    path: Path,
    records: tuple[tuple[str, bool, int, str, str], ...] | None = None,
    sa_tags: dict[tuple[str, bool], str] | None = None,
    reverse_records: set[tuple[str, bool]] | None = None,
    chromosome: str = "chr1",
    record_flags: dict[tuple[str, bool, int], int] | None = None,
) -> None:
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": chromosome, "LN": 1_000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as alignment:
        records = records or (
            ("rna_t", True, 100, "50M100N50M", "T"),
            ("rna_t", False, 400, "100M", "T"),
            ("rna_not_both_t", True, 100, "50M100N50M", "T"),
            ("rna_not_both_t", False, 400, "100M", "A"),
            ("ordinary_t_pair", True, 400, "100M", "T"),
            ("ordinary_t_pair", False, 500, "100M", "T"),
            ("incomplete_rna_pair", True, 100, "50M100N50M", "T"),
        )
        for name, is_read1, start, cigar, first_base in records:
            read = pysam.AlignedSegment(alignment.header)
            read.query_name = name
            reverse_flag = (
                0x10
                if reverse_records and (name, is_read1) in reverse_records
                else 0
            )
            additional_flags = (record_flags or {}).get((name, is_read1, start), 0)
            read.flag = (
                0x1
                | (0x40 if is_read1 else 0x80)
                | reverse_flag
                | additional_flags
            )
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = cigar
            query_length = sum(
                length
                for operation, length in read.cigartuples
                if operation
                in {
                    pysam.CMATCH,
                    pysam.CINS,
                    pysam.CSOFT_CLIP,
                    pysam.CEQUAL,
                    pysam.CDIFF,
                }
            )
            forward_sequence = first_base + ("C" * (query_length - 1))
            if read.is_reverse:
                complement = str.maketrans("ACGTN", "TGCAN")
                read.query_sequence = forward_sequence.translate(complement)[::-1]
            else:
                read.query_sequence = forward_sequence
            read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
            if sa_tags and (name, is_read1) in sa_tags:
                read.set_tag("SA", sa_tags[(name, is_read1)])
            alignment.write(read)


def _write_gtf(path: Path, chromosome: str = "chr1") -> None:
    path.write_text(
        f'{chromosome}\ttest\texon\t101\t150\t.\t+\t.\ttranscript_id "tx1";\n'
        f'{chromosome}\ttest\texon\t251\t300\t.\t+\t.\ttranscript_id "tx1";\n'
    )


def _run_softspot(
    alignment: Path,
    annotation: Path,
    *options: str,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(REPOSITORY / "src")
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "softspot",
            str(alignment),
            str(annotation),
            *options,
        ],
        cwd=REPOSITORY,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


class SoftspotCliTests(unittest.TestCase):
    def test_reports_fraction_of_splice_supporting_pairs_with_two_t_starts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(alignment)
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(0.5, 1, 2)),
            msg=completed.stderr,
        )

    def test_large_deletion_at_annotated_intron_supports_rna_origin(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("rna_deletion", True, 100, "50M100D50M", "T"),
                    ("rna_deletion", False, 400, "100M", "T"),
                ),
            )
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )

    def test_softclips_at_annotated_exon_boundaries_support_rna_origin(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("left_boundary", True, 100, "5H5S50M", "T"),
                    ("left_boundary", False, 400, "100M", "T"),
                    ("right_boundary", True, 100, "50M5S5H", "T"),
                    ("right_boundary", False, 400, "100M", "T"),
                ),
            )
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 2, 2)),
            msg=completed.stderr,
        )

    def test_split_alignment_at_annotated_junction_supports_rna_origin(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("rna_split", True, 100, "50M", "T"),
                    ("rna_split", False, 400, "100M", "T"),
                ),
                sa_tags={("rna_split", True): "chr1,251,+,50M,60,0;"},
            )
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )

    def test_pair_suffixes_and_reverse_orientation_preserve_t_start(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("rna_suffix/1", True, 100, "50M100N50M", "T"),
                    ("rna_suffix/2", False, 400, "100M", "T"),
                ),
                reverse_records={("rna_suffix/1", True)},
            )
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )

    def test_chr_prefix_difference_between_gtf_and_alignment_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("rna_alias", True, 100, "50M100N50M", "T"),
                    ("rna_alias", False, 400, "100M", "T"),
                ),
                chromosome="1",
            )
            _write_gtf(annotation, chromosome="chr1")

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )

    def test_gzipped_gtf_annotation_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf.gz"
            plain_annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("rna_gzip", True, 100, "50M100N50M", "T"),
                    ("rna_gzip", False, 400, "100M", "T"),
                ),
            )
            _write_gtf(plain_annotation)
            with gzip.open(annotation, "wt") as compressed:
                compressed.write(plain_annotation.read_text())

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )

    def test_short_deletion_at_annotated_junction_supports_rna_origin(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("rna_short_deletion", True, 100, "50M50D50M", "T"),
                    ("rna_short_deletion", False, 400, "100M", "T"),
                ),
            )
            annotation.write_text(
                'chr1\ttest\texon\t101\t150\t.\t+\t.\ttranscript_id "tx1";\n'
                'chr1\ttest\texon\t201\t250\t.\t+\t.\ttranscript_id "tx1";\n'
            )

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )

    def test_near_miss_splice_signals_do_not_enter_the_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            records = (
                ("near_n", True, 100, "50M99N50M", "T"),
                ("near_n", False, 400, "100M", "T"),
                ("short_d", True, 100, "50M99D50M", "T"),
                ("short_d", False, 400, "100M", "T"),
                ("off_boundary", True, 101, "5S50M", "T"),
                ("off_boundary", False, 400, "100M", "T"),
                ("near_sa", True, 100, "50M", "T"),
                ("near_sa", False, 400, "100M", "T"),
            )
            _write_alignment(
                alignment,
                records=records,
                sa_tags={("near_sa", True): "chr1,252,+,50M,60,0;"},
            )
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(0, 0, 0)),
            msg=completed.stderr,
        )

    def test_splice_and_t_evidence_accumulate_across_alignments(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(
                alignment,
                records=(
                    ("multi_alignment", True, 400, "100M", "T"),
                    ("multi_alignment", True, 100, "50M100N50M", "A"),
                    ("multi_alignment", False, 500, "100M", "T"),
                ),
                record_flags={("multi_alignment", True, 100): 0x800},
            )
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )

    def test_default_output_is_a_headered_tsv_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(alignment)
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (completed.returncode, completed.stdout, completed.stderr),
            (0, _expected_output(0.5, 1, 2), ""),
        )

    def test_verbose_option_is_not_available(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(alignment)
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation, "--verbose")

        self.assertEqual(completed.returncode, 2)
        self.assertIn("No such option '--verbose'", completed.stderr)

    def test_min_intron_option_is_not_available(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            _write_alignment(alignment)
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation, "--min-intron", "60")

        self.assertEqual(completed.returncode, 2)
        self.assertIn("No such option '--min-intron'", completed.stderr)

    def test_invalid_alignment_fails_without_a_python_traceback(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "broken.bam"
            annotation = Path(tmpdir) / "genes.gtf"
            alignment.write_text("not an alignment")
            _write_gtf(annotation)

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (
                completed.returncode,
                completed.stdout,
                "Error: Unable to read alignment" in completed.stderr,
                "Traceback" in completed.stderr,
            ),
            (1, "", True, False),
        )

    def test_truncated_gzip_annotation_fails_without_aborting(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / "sample.bam"
            annotation = Path(tmpdir) / "broken.gtf.gz"
            _write_alignment(alignment)
            annotation.write_bytes(gzip.compress(b"gtf contents")[:-8])

            completed = _run_softspot(alignment, annotation)

        self.assertEqual(
            (
                completed.returncode,
                completed.stdout,
                "Error: Unable to read alignment" in completed.stderr,
                "Aborted!" in completed.stderr,
            ),
            (1, "", True, False),
        )

    def test_cram_input_uses_the_explicit_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = Path(tmpdir) / "sample.bam"
            cram = Path(tmpdir) / "sample.cram"
            annotation = Path(tmpdir) / "genes.gtf"
            reference = Path(tmpdir) / "reference.fa"
            _write_alignment(
                bam,
                records=(
                    ("rna_cram", True, 100, "50M100N50M", "T"),
                    ("rna_cram", False, 400, "100M", "T"),
                ),
            )
            _write_gtf(annotation)
            reference.write_text(">chr1\n" + ("C" * 1_000) + "\n")
            pysam.faidx(str(reference))
            with pysam.AlignmentFile(bam, "rb") as source:
                with pysam.AlignmentFile(
                    cram,
                    "wc",
                    template=source,
                    reference_filename=str(reference),
                ) as destination:
                    for read in source.fetch(until_eof=True):
                        destination.write(read)

            completed = _run_softspot(
                cram,
                annotation,
                "--reference",
                str(reference),
            )

        self.assertEqual(
            (completed.returncode, completed.stdout),
            (0, _expected_output(1, 1, 1)),
            msg=completed.stderr,
        )


if __name__ == "__main__":
    unittest.main()
