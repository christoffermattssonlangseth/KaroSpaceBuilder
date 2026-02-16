from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
ENTRYPOINT = ROOT / "packaging" / "entrypoints" / "karospace_builder_gui.py"
LOGO = ROOT / "assets" / "logo_KSB.png"


def _pyinstaller_args(onefile: bool) -> list[str]:
    args = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--clean",
        "--windowed",
        "--name",
        "KaroSpaceBuilder",
    ]
    if onefile:
        args.append("--onefile")

    # AnnData/scientific stacks often need broader collection to avoid missing
    # runtime modules in frozen builds.
    for module in ("karospace", "scanpy", "anndata", "numpy", "pandas", "scipy", "h5py", "PIL"):
        args.extend(["--collect-submodules", module])

    for module in ("karospace", "scanpy", "anndata", "pandas", "scipy", "PIL"):
        args.extend(["--collect-data", module])

    if LOGO.exists():
        sep = ";" if os.name == "nt" else ":"
        args.extend(["--add-data", f"{LOGO}{sep}assets"])

    args.append(str(ENTRYPOINT))
    return args


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build KaroSpaceBuilder desktop app with PyInstaller.")
    parser.add_argument(
        "--onefile",
        action="store_true",
        help="Build a onefile executable (default is onedir/app bundle).",
    )
    parser.add_argument(
        "--keep-build-dirs",
        action="store_true",
        help="Do not delete existing build/dist directories before packaging.",
    )
    args = parser.parse_args(argv)

    if not ENTRYPOINT.exists():
        raise FileNotFoundError(f"Missing PyInstaller entrypoint: {ENTRYPOINT}")

    if not args.keep_build_dirs:
        for path in (ROOT / "build", ROOT / "dist"):
            if path.exists():
                shutil.rmtree(path)

    cmd = _pyinstaller_args(onefile=args.onefile)
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)

    dist_dir = ROOT / "dist"
    if not dist_dir.exists():
        raise RuntimeError("PyInstaller finished but dist/ was not created.")

    print(f"Build complete: {dist_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
