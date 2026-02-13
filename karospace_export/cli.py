from __future__ import annotations

import argparse
from pathlib import Path
import http.server
import socketserver

from .export import ExportConfig, export_h5ad


def _parse_bool(raw: str) -> bool:
    value = raw.strip().lower()
    if value in {"true", "1", "yes", "y"}:
        return True
    if value in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("Expected true/false")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="karospace-export",
        description="Export an AnnData .h5ad into a static KaroSpace-compatible viewer bundle.",
    )
    parser.add_argument("--h5ad", required=True, help="Path to input .h5ad file")
    parser.add_argument("--outdir", required=True, help="Output directory for static export")
    parser.add_argument(
        "--coords",
        choices=["obsm:spatial", "obs:centroid_x_y"],
        default=None,
        help="Coordinate source. Auto-detected if omitted.",
    )
    parser.add_argument(
        "--anno",
        action="append",
        default=None,
        help="Annotation column from obs. Repeatable.",
    )
    parser.add_argument(
        "--genes",
        default="hvgs:500",
        help="Gene mode: hvgs:N, top_mean:N, or list:PATH",
    )
    parser.add_argument("--image", default=None, help="Optional tissue image path")
    parser.add_argument("--downsample", type=int, default=None, help="Downsample to N cells")
    parser.add_argument(
        "--gzip",
        type=_parse_bool,
        default=True,
        help="Compress text/binary assets with gzip (true/false, default: true)",
    )
    parser.add_argument(
        "--max-asset-mb",
        type=float,
        default=16.0,
        help="Split large assets so each file stays below this size",
    )
    parser.add_argument(
        "--serve",
        action="store_true",
        help="Start a tiny local http server in outdir after export",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port for --serve (default: 8000)",
    )
    parser.add_argument(
        "--no-preview",
        action="store_true",
        help="Disable preview.png generation",
    )
    return parser


def _serve_directory(directory: Path, port: int) -> None:
    handler = http.server.SimpleHTTPRequestHandler

    class ReuseTCPServer(socketserver.TCPServer):
        allow_reuse_address = True

    with ReuseTCPServer(("", port), handler) as httpd:
        print(f"Serving {directory} at http://127.0.0.1:{port}")
        print("Press Ctrl+C to stop")
        prev = Path.cwd()
        try:
            import os

            os.chdir(directory)
            httpd.serve_forever()
        finally:
            os.chdir(prev)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    config = ExportConfig(
        h5ad_path=Path(args.h5ad),
        outdir=Path(args.outdir),
        coords=args.coords,
        annotation_columns=args.anno,
        genes_mode=args.genes,
        image_path=Path(args.image) if args.image else None,
        downsample=args.downsample,
        gzip=args.gzip,
        max_asset_mb=args.max_asset_mb,
        preview=not args.no_preview,
    )

    manifest = export_h5ad(config)

    outdir = Path(args.outdir)
    print(f"Export complete: {outdir}")
    print(f"Cells: {manifest.n_cells}")
    print(f"Genes exported: {manifest.n_genes_exported}")
    print(f"Manifest: {outdir / 'manifest.json'}")

    if args.serve:
        _serve_directory(outdir, port=args.port)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
