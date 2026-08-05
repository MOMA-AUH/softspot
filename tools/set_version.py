#!/usr/bin/env python3
"""Update the softspot version in all package metadata.

Usage:
    python tools/set_version.py 0.2.0
"""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
VERSION_PATTERN = re.compile(r"^\d+\.\d+\.\d+(?:[a-zA-Z0-9_.+-]*)?$")

FILES_AND_PATTERNS = (
    (
        ROOT / "pyproject.toml",
        re.compile(r'^(version\s*=\s*")([^"]+)(")$', flags=re.MULTILINE),
    ),
    (
        ROOT / "devtools" / "conda-recipe" / "meta.yaml",
        re.compile(
            r'^(\{\%\s*set\s+version\s*=\s*")([^"]+)("\s*\%\})$',
            flags=re.MULTILINE,
        ),
    ),
    (
        ROOT / "src" / "softspot" / "__init__.py",
        re.compile(r'^(__version__\s*=\s*")([^"]+)(")$', flags=re.MULTILINE),
    ),
)


def _replace_version_once(
    text: str,
    pattern: re.Pattern[str],
    new_version: str,
    path: Path,
) -> tuple[str, str]:
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one version match in {path}, found {len(matches)}"
        )

    old_version = matches[0].group(2)
    updated = pattern.sub(rf"\g<1>{new_version}\g<3>", text, count=1)
    return updated, old_version


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("Usage: python tools/set_version.py <new-version>")
        return 2

    new_version = argv[1].strip()
    if not VERSION_PATTERN.fullmatch(new_version):
        print(f"Invalid version: {new_version}")
        return 2

    updates: list[tuple[Path, str]] = []
    old_versions: set[str] = set()

    # Validate every target before writing any of them, avoiding partial updates.
    for path, pattern in FILES_AND_PATTERNS:
        text = path.read_text(encoding="utf-8")
        updated, old_version = _replace_version_once(
            text, pattern, new_version, path
        )
        updates.append((path, updated))
        old_versions.add(old_version)

    if len(old_versions) != 1:
        versions = ", ".join(sorted(old_versions))
        print(f"Refusing to update inconsistent existing versions: {versions}")
        return 1

    for path, updated in updates:
        path.write_text(updated, encoding="utf-8")

    old_version = next(iter(old_versions))
    print(f"Updated version {old_version} -> {new_version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
