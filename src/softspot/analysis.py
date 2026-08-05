from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import TypeVar

import pysam


_READ1 = 0x01
_READ2 = 0x02
_T_READ1 = 0x04
_T_READ2 = 0x08
_REFERENCE_CIGAR_SPAN = re.compile(r"(\d+)([MDN=X])")
_EMPTY = frozenset()
_AnnotationValue = TypeVar("_AnnotationValue")


def _pair_name(query_name: str) -> str:
    if query_name.endswith(("/1", "/2")):
        return query_name[:-2]
    return query_name


@dataclass(frozen=True)
class AnalysisResult:
    splice_supporting_pairs: int
    t_starting_splice_supporting_pairs: int

    @property
    def score(self) -> float:
        if self.splice_supporting_pairs == 0:
            return 0.0
        return self.t_starting_splice_supporting_pairs / self.splice_supporting_pairs


@dataclass(frozen=True)
class _SpliceAnnotation:
    junctions: dict[str, frozenset[tuple[int, int]]]
    exon_starts: dict[str, frozenset[int]]
    exon_ends: dict[str, frozenset[int]]


def _with_chromosome_aliases(
    values: dict[str, set[_AnnotationValue]],
) -> dict[str, frozenset[_AnnotationValue]]:
    result: dict[str, frozenset[_AnnotationValue]] = {}
    for chromosome in values:
        alias = (
            chromosome[3:]
            if chromosome.startswith("chr")
            else f"chr{chromosome}"
        )
        merged = frozenset(values.get(chromosome, set())) | frozenset(
            values.get(alias, set())
        )
        result[chromosome] = merged
        result[alias] = merged
    return result


def _is_eligible(read: pysam.AlignedSegment) -> bool:
    return (
        read.is_paired
        and (read.is_read1 or read.is_read2)
        and not read.is_unmapped
        and bool(read.cigartuples)
    )


def _load_annotation(annotation: str | Path) -> _SpliceAnnotation:
    exons_by_transcript: dict[
        tuple[str, str], list[tuple[int, int]]
    ] = defaultdict(list)
    annotation_path = Path(annotation)
    opener = gzip.open if annotation_path.name.lower().endswith(".gz") else open
    with opener(annotation_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            transcript_field = next(
                (
                    item.strip()
                    for item in fields[8].split(";")
                    if item.strip().startswith("transcript_id ")
                ),
                None,
            )
            if transcript_field is None:
                continue
            transcript = transcript_field.split(maxsplit=1)[1].strip('"')
            start = int(fields[3]) - 1
            end = int(fields[4])
            exons_by_transcript[(fields[0], transcript)].append((start, end))

    junctions: dict[str, set[tuple[int, int]]] = defaultdict(set)
    exon_starts: dict[str, set[int]] = defaultdict(set)
    exon_ends: dict[str, set[int]] = defaultdict(set)
    for (chromosome, _), exons in exons_by_transcript.items():
        ordered = sorted(exons)
        for start, end in ordered:
            exon_starts[chromosome].add(start)
            exon_ends[chromosome].add(end)
        for left, right in zip(ordered, ordered[1:]):
            if right[0] > left[1]:
                junctions[chromosome].add((left[1], right[0]))
    return _SpliceAnnotation(
        junctions=_with_chromosome_aliases(junctions),
        exon_starts=_with_chromosome_aliases(exon_starts),
        exon_ends=_with_chromosome_aliases(exon_ends),
    )


def _supports_splice(
    read: pysam.AlignedSegment,
    junctions: frozenset[tuple[int, int]],
    exon_starts: frozenset[int],
    exon_ends: frozenset[int],
) -> bool:
    position = read.reference_start
    for operation, length in read.cigartuples or ():
        if operation == pysam.CREF_SKIP and (position, position + length) in junctions:
            return True
        if operation == pysam.CDEL and (position, position + length) in junctions:
            return True
        if operation in {
            pysam.CMATCH,
            pysam.CDEL,
            pysam.CREF_SKIP,
            pysam.CEQUAL,
            pysam.CDIFF,
        }:
            position += length

    operations = list(read.cigartuples or ())
    while operations and operations[0][0] == pysam.CHARD_CLIP:
        operations.pop(0)
    while operations and operations[-1][0] == pysam.CHARD_CLIP:
        operations.pop()
    if operations:
        if (
            operations[0][0] == pysam.CSOFT_CLIP
            and read.reference_start in exon_starts
        ):
            return True
        if operations[-1][0] == pysam.CSOFT_CLIP and read.reference_end in exon_ends:
            return True

    if read.has_tag("SA"):
        for entry in read.get_tag("SA").rstrip(";").split(";"):
            fields = entry.split(",")
            if len(fields) < 4 or fields[0] != read.reference_name:
                continue
            try:
                split_start = int(fields[1]) - 1
            except ValueError:
                continue
            split_end = split_start + sum(
                int(length)
                for length, _ in _REFERENCE_CIGAR_SPAN.findall(fields[3])
            )
            if (
                (read.reference_end, split_start) in junctions
                or (split_end, read.reference_start) in junctions
            ):
                return True
    return False


def analyze_alignment(
    alignment: str | Path,
    annotation: str | Path,
    *,
    reference: str | Path | None = None,
    threads: int = 1,
) -> AnalysisResult:
    splice_annotation = _load_annotation(annotation)
    open_options: dict[str, object] = {"threads": threads}
    if reference is not None:
        open_options["reference_filename"] = str(reference)
    mode = "rb" if str(alignment).lower().endswith(".bam") else "rc"

    splice_supporting_names: set[str] = set()
    with pysam.AlignmentFile(str(alignment), mode, **open_options) as records:
        for read in records.fetch(until_eof=True):
            if not _is_eligible(read):
                continue

            chromosome = read.reference_name
            if _supports_splice(
                read,
                splice_annotation.junctions.get(chromosome, _EMPTY),
                splice_annotation.exon_starts.get(chromosome, _EMPTY),
                splice_annotation.exon_ends.get(chromosome, _EMPTY),
            ):
                splice_supporting_names.add(_pair_name(read.query_name))

    pair_flags = {name: 0 for name in splice_supporting_names}
    del splice_supporting_names
    del splice_annotation
    with pysam.AlignmentFile(str(alignment), mode, **open_options) as records:
        for read in records.fetch(until_eof=True):
            if not _is_eligible(read):
                continue
            pair_name = _pair_name(read.query_name)
            if pair_name not in pair_flags:
                continue

            flags = pair_flags.get(pair_name, 0)
            sequence = read.get_forward_sequence()
            starts_with_t = bool(sequence) and sequence[0].upper() == "T"
            if read.is_read1:
                flags |= _READ1
                if starts_with_t:
                    flags |= _T_READ1
            else:
                flags |= _READ2
                if starts_with_t:
                    flags |= _T_READ2
            pair_flags[pair_name] = flags

    complete = _READ1 | _READ2
    both_t = _T_READ1 | _T_READ2
    splice_pairs = 0
    final_pairs = 0
    for flags in pair_flags.values():
        if flags & complete != complete:
            continue
        splice_pairs += 1
        if flags & both_t == both_t:
            final_pairs += 1

    return AnalysisResult(
        splice_supporting_pairs=splice_pairs,
        t_starting_splice_supporting_pairs=final_pairs,
    )
