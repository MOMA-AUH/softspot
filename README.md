# softspot

[![Conda Version](https://img.shields.io/conda/vn/MOMA-AUH/softspot?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/softspot) [![Conda Downloads](https://img.shields.io/conda/dn/MOMA-AUH/softspot?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/softspot)

The `softspot` tool measures a prepared RNA-library signature in paired-end DNA sequencing
data. It combines two independent observations:

1. A read pair has alignment evidence supporting an annotated splice event.
2. Both R1 and R2 start with `T`, as expected from an RNA-library chemistry with
   T-overhangs.

The command reports one number: the fraction of complete, splice-supporting read
pairs for which both mates start with `T`.

`softspot` does not apply a calibrated threshold and does not emit a
contaminated/not-contaminated verdict. The result is an evidence score that must
be interpreted in the context of the library preparation and a locally validated
workflow.

## Installation

The recommended installation is from the `MOMA-AUH` conda channel:

```bash
conda create -n softspot MOMA-AUH::softspot
conda activate softspot
```

`softspot` requires Python 3.10 or newer. It can also be installed from a source
checkout with pip:

```bash
python -m pip install .
```

## Usage

```text
softspot [OPTIONS] ALIGNMENT ANNOTATION
```

`ALIGNMENT` is a paired-end BAM or CRAM file. `ANNOTATION` is a GTF file, which
may be plain text or gzip-compressed (`.gtf.gz`). For example:

```bash
softspot sample.bam genes.gtf
```

Supply the matching reference FASTA when decoding CRAM:

```bash
softspot sample.cram genes.gtf.gz \
  --reference reference.fa \
  --threads 4
```

The available options are:

| Option | Meaning |
|:--|:--|
| `-T, --reference PATH` | Reference FASTA used to decode CRAM input. |
| `--threads INTEGER` | Alignment decompression threads; must be at least 1. Default: `1`. |
| `--version` | Print the installed version and exit. |
| `--help` | Show command help and exit. |

### Output

On success, standard output is a headered tab-separated table. The `fraction` is
formatted with six decimal places; `positive_pairs` is the numerator and
`eligible_pairs` is the denominator. This is both human-readable and
straightforward to parse in a workflow:

```console
$ softspot sample.bam genes.gtf
fraction        positive_pairs  eligible_pairs
0.125000        100             800
```

For example, extract the fraction with `awk 'NR == 2 { print $1 }'`.

## What the score means

The score is calculated as:

```text
        complete splice-supporting pairs whose R1 and R2 both start with T
score = ------------------------------------------------------------------
                     all complete splice-supporting pairs
```

A *complete* pair is one for which eligible R1 and R2 records were both observed.
A pair is *splice-supporting* when at least one of its eligible alignments matches
one of the four annotation-based splice heuristics described below.

For each mate, `softspot` checks the first base returned by pysam's
`get_forward_sequence()`. This restores the original, as-sequenced orientation
before testing for `T`, including for alignments mapped to the reverse strand.
Both R1 and R2 must pass this test for the pair to enter the numerator.

If the denominator is zero, `softspot` reports a score of `0.000000` and zero
counts. The counts distinguish "zero observed signal" from "no splice-supporting
complete pairs were available."

The score is conditional on splice support. It is not the fraction of all read
pairs, and it is not a calibrated estimate of the proportion of contaminating
library molecules.

## Splice-support heuristics

`softspot` derives exon boundaries and introns from consecutive exons belonging
to the same transcript in the GTF. It then recognises exactly four forms of splice
support:

1. **Reference skip (`N`)**: a CIGAR `N` operation spans an annotated junction
   exactly.
2. **Deletion (`D`)**: a CIGAR `D` operation spans an annotated junction exactly.
3. **Exon-boundary soft clip (`S`)**: after ignoring terminal hard clips, a
   leading soft clip starts at an annotated exon start or a trailing soft clip
   ends at an annotated exon end.
4. **Supplementary alignment (`SA`)**: a same-chromosome split alignment encoded
   in the `SA` tag has primary and supplementary endpoints matching an annotated
   junction.

The soft-clipped sequence itself is not tested for `T`. Soft clipping supplies
splice evidence; the separate first-base test on both mates supplies the
T-overhang signature.

Only paired, mapped R1/R2 records with a CIGAR are eligible. Pairing uses the
exact query name after removing a terminal `/1` or `/2`, if present. Secondary,
supplementary, duplicate and vendor-QC-failed records are not filtered out;
qualifying evidence is accumulated by pair name across all eligible records.

GTF parsing uses `exon` records with a `transcript_id` attribute. Chromosome names
with and without a `chr` prefix are treated as aliases.

## Performance

The alignment is scanned sequentially twice:

1. The first pass stores the exact names of pairs with splice-supporting
   alignments.
2. The second pass records R1/R2 presence and T-start state only for those names.

Using exact names avoids hash-collision errors. The two-pass design also keeps
pair-state memory proportional to the number of splice-supporting names rather
than every pair in the alignment. Its trade-off is reading and decompressing the
alignment twice. `--threads` controls decompression threading for both passes.
The sequential whole-file scan does not require a BAM or CRAM index.

## Interpretation and limitations

- The method is intended for contamination by already prepared RNA-library
  molecules whose library chemistry produces T-overhangs. It is not a general
  assay for RNA molecules present before library preparation.
- The relevant RNA kit and expected background should be established for the
  local workflow. Kits without this both-mates-T signature are outside the
  method's scope.
- There is no universal cutoff in `softspot`. Comparisons and decision thresholds
  require independently validated controls generated with a representative
  protocol.
- A low denominator makes the score unstable. Retain the reported counts
  alongside the score.
- Splice evidence is annotation-dependent. Unannotated splice events are missed,
  and an incomplete or mismatched GTF changes the denominator.
- The `N` and deletion heuristics require an exact annotated junction match. The
  `SA` heuristic considers same-chromosome splits only.
- Incomplete pairs, unmapped records and records without CIGAR operations do not
  contribute. No mapping-quality, duplicate, vendor-QC-fail, secondary or
  supplementary filter is applied.
- The result is a signature score, not a contamination percentage or a clinical
  classification.

## Development and testing

Create an editable development environment from a source checkout:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --editable . coverage
```

Run the unit tests and coverage report with the standard library test runner:

```bash
coverage run -m unittest discover -s tests -v
coverage report
```

The tests exercise the four splice heuristics, forward/reverse read orientation,
pair-name suffixes, chromosome aliases, compressed GTF input, CLI output and
error handling. Continuous integration runs them on Python 3.10 through 3.14 and
also builds the wheel and source distribution and smoke-tests the wheel.

To build the Python distributions locally:

```bash
python -m pip install build
python -m build
```

## Releases

Pushing a tag named `v<version>` starts `.github/workflows/conda-release.yml`.
The workflow builds and smoke-tests the conda package on a GitHub-hosted runner,
then publishes the artifact to the `MOMA-AUH` Anaconda channel. Publishing
requires the repository secret `ANACONDA_API_TOKEN`; existing artifacts are
skipped.

Prepare version changes on a release branch and open a pull request to `master`:

```bash
git switch master
git pull --ff-only origin master
git switch -c release/v0.2.0
python tools/set_version.py 0.2.0
coverage run -m unittest discover -s tests -v
git add pyproject.toml devtools/conda-recipe/meta.yaml src/softspot/__init__.py
git commit -m "Release v0.2.0"
git push -u origin release/v0.2.0
```

After the required checks pass, squash-merge the release pull request. Tag the
resulting commit on `master` and push only the release tag:

```bash
git switch master
git pull --ff-only origin master
git tag v0.2.0
git push origin v0.2.0
```

Releases are published to <https://anaconda.org/MOMA-AUH/softspot>.

## License

`softspot` is distributed under the MIT License. See `LICENSE`.
