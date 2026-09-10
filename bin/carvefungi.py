#!/usr/bin/env python3
"""Run the original reconstruction pipeline in a writable, isolated workspace."""

import argparse
import os
from pathlib import Path
import re
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id", required=True, help="Model identifier (letters, digits, underscores)")
    parser.add_argument("--eggnog", required=True, type=Path, help="Legacy eggNOG-mapper annotation file")
    parser.add_argument("--localization", required=True, type=Path, help="Matching CarveFungi .loc_pred file")
    parser.add_argument("--output", required=True, type=Path, help="Dedicated output directory for this run")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", args.id):
        parser.error("--id must start with a letter or underscore and contain only letters, digits, underscores")
    eggnog, localization = args.eggnog.resolve(), args.localization.resolve()
    for path in (eggnog, localization):
        if not path.is_file():
            parser.error(f"Input file not found: {path}")
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("--output must be a new or empty directory; upstream resume handling is not supported")
    output.mkdir(parents=True, exist_ok=True)

    data = Path(__file__).resolve().parent.parent / "data"
    previous = Path.cwd()
    with tempfile.TemporaryDirectory(prefix="carvefungi-") as directory:
        workspace = Path(directory)
        (workspace / "bin").mkdir()
        (workspace / "data").symlink_to(data, target_is_directory=True)
        (workspace / "results").symlink_to(output, target_is_directory=True)
        # COBRApy creates its cache on import. Keep it writable and per-job,
        # even with a read-only container/home and a shared cluster /tmp.
        previous_cache = os.environ.get("XDG_CACHE_HOME")
        os.environ["XDG_CACHE_HOME"] = str(workspace / "cache")
        try:
            import cobra
            import CreateModelEggNogPool

            # SCIP handles reconstruction; GLPK handles COBRApy's FBA.
            cobra.Configuration().solver = "glpk"
            os.chdir(workspace / "bin")
            CreateModelEggNogPool.carveFungi_pipeline(args.id, str(eggnog), str(localization))
        finally:
            os.chdir(previous)
            if previous_cache is None:
                os.environ.pop("XDG_CACHE_HOME", None)
            else:
                os.environ["XDG_CACHE_HOME"] = previous_cache


if __name__ == "__main__":
    main()
