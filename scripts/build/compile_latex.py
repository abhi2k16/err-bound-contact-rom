"""Compile a LaTeX source file to PDF with latexmk."""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    default_source = Path(__file__).resolve().parents[2] / "manuscript_master_v2.tex"
    parser = argparse.ArgumentParser(
        description="Compile a LaTeX source file to PDF using latexmk."
    )
    parser.add_argument(
        "source",
        nargs="?",
        type=Path,
        default=default_source,
        help="LaTeX source file (default: manuscript_master_v2.tex in the repository root)",
    )
    args = parser.parse_args()
    source = args.source.expanduser().resolve()

    if not source.is_file():
        parser.error("LaTeX source file not found: {}".format(source))

    latexmk = shutil.which("latexmk")
    if latexmk is None:
        print(
            "Error: latexmk was not found. Install a LaTeX distribution and add "
            "latexmk to PATH.",
            file=sys.stderr,
        )
        return 1

    result = subprocess.run(
        [
            latexmk,
            "-pdf",
            "-interaction=nonstopmode",
            "-file-line-error",
            source.name,
        ],
        cwd=str(source.parent),
        check=False,
    )
    if result.returncode != 0:
        print(
            "LaTeX compilation failed with exit code {}.".format(
                result.returncode
            ),
            file=sys.stderr,
        )
        return result.returncode

    pdf = source.with_suffix(".pdf")
    if not pdf.is_file():
        print("Error: compilation completed without creating {}.".format(pdf), file=sys.stderr)
        return 1

    print("Created {}".format(pdf))
    return 0


if __name__ == "__main__":
    sys.exit(main())
