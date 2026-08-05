from __future__ import annotations

from pathlib import Path

import click

from softspot import __version__
from softspot.analysis import analyze_alignment


@click.command()
@click.argument(
    "alignment",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "annotation",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--reference",
    "-T",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Reference FASTA used to decode CRAM input.",
)
@click.option(
    "--threads",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Alignment decompression threads.",
)
@click.version_option(version=__version__)
def main(
    alignment: Path,
    annotation: Path,
    reference: Path | None,
    threads: int,
) -> None:
    """Report the prepared RNA-library signature in ALIGNMENT using ANNOTATION."""
    try:
        result = analyze_alignment(
            alignment,
            annotation,
            reference=reference,
            threads=threads,
        )
    except (EOFError, OSError, ValueError) as error:
        raise click.ClickException(
            f"Unable to read alignment or annotation: {error}"
        ) from error
    click.echo(
        "fraction\tpositive_pairs\teligible_pairs"
    )
    click.echo(
        f"{result.score:.6f}\t{result.t_starting_splice_supporting_pairs}"
        f"\t{result.splice_supporting_pairs}"
    )
