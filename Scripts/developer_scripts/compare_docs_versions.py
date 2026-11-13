#!/usr/bin/env python3
"""
compare_docs_versions.py

Compare two CGAL documentation trees (n-1 vs n) and list .html pages
that existed in n-1 but are missing in n. You can provide local paths
or ask the script to download the official doc_html archives from GitHub.

Examples:
    # Using local paths
    python compare_docs_versions.py \
      --old-path /path/to/5.6.2 \
      --new-path /path/to/6.0.1

    # Auto-download from GitHub releases and compare
    python compare_docs_versions.py \
      --old-version 5.6.2 \
      --new-version 6.0.1 \
      --download-dir ~/cgal-docs

    # Download new version and compare with local old version
    python compare_docs_versions.py \
      --old-path /path/to/5.6.2 \
      --new-version 6.0.1 \
      --download-dir ~/cgal-docs

Notes:
    - Compares by *relative path* (e.g. Manual/packages.html, Alpha_shapes_3/index.html).
    - Only pages under Manual/ or <PackageName>/ (starting with uppercase) are considered.
Dependencies:
    pip install requests
"""

from __future__ import annotations
import argparse
import os
import sys
import tarfile
import time
from pathlib import Path
from typing import Optional, Set, List
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests


DEFAULT_URL_TEMPLATE = (
    "https://github.com/CGAL/cgal/releases/download/v{version}/CGAL-{version}-doc_html.tar.xz"
)


def info(msg: str) -> None:
    print(f"[INFO] {msg}")


def warn(msg: str) -> None:
    print(f"[WARN] {msg}")


def err(msg: str) -> None:
    print(f"[ERROR] {msg}")


# -------------------------
# Download / extract helpers (requests-only)
# -------------------------

def stream_download(url: str, dest: Path, chunk_size: int = 1 << 20, retries: int = 3, timeout: int = 30) -> None:
    """
    Download 'url' to 'dest' using requests (streaming). Retries on transient failures.
    """
    if dest.exists():
        info(f"Download target already exists, skipping: {dest}")
        return

    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                total = int(r.headers.get("Content-Length", 0))
                read = 0
                started = time.time()
                with open(dest, "wb") as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            read += len(chunk)
                            if total:
                                pct = (read / total) * 100
                                speed = read / max(1e-6, (time.time() - started))
                                print(f"\r[DL] {dest.name}: {read}/{total} bytes ({pct:.1f}%) ~ {int(speed/1024)} KiB/s", end="")
                if total:
                    print()  # newline after progress
            return
        except (requests.RequestException, OSError) as e:
            last_exc = e
            warn(f"Download attempt {attempt}/{retries} failed: {e}")
            time.sleep(1 + 0.5 * attempt)
    raise RuntimeError(f"Failed to download {url}: {last_exc}")


def safe_extract_xz(archive_path: Path, extract_to: Path, strip_components: int = 1) -> None:
    """
    Extract .tar.xz archive to extract_to, stripping 'strip_components' leading path components.
    Prevents path traversal.
    """
    if not archive_path.is_file():
        raise FileNotFoundError(str(archive_path))
    extract_to.mkdir(parents=True, exist_ok=True)
    info(f"Extracting {archive_path.name} to {extract_to} ...")

    with tarfile.open(archive_path, mode="r:xz") as tar:
        for m in tar.getmembers():
            parts = Path(m.name).parts
            if len(parts) <= strip_components:
                continue
            new_rel = Path(*parts[strip_components:])
            if not new_rel:
                continue
            final_path = (extract_to / new_rel).resolve()
            if not str(final_path).startswith(str(extract_to.resolve())):
                raise RuntimeError(f"Blocked suspicious path: {final_path}")
            m.name = str(new_rel)
            tar.extract(m, path=extract_to)


def ensure_doc_tree_from_version(version: str, download_dir: Path, url_template: str = DEFAULT_URL_TEMPLATE) -> Path:
    """
    Ensure doc_html for 'version' is available under download_dir/<version>.
    Downloads & extracts if missing. Returns the root path (contains Manual/).
    """
    root = (download_dir / version).resolve()
    manual_index = root / "Manual" / "index.html"
    if manual_index.is_file():
        info(f"Doc tree for {version} already present at {root}")
        return root

    download_dir.mkdir(parents=True, exist_ok=True)
    url = url_template.format(version=version)
    archive_path = download_dir / f"CGAL-{version}-doc_html.tar.xz"

    info(f"Fetching doc_html for CGAL {version}")
    info(f"URL: {url}")
    stream_download(url, archive_path)

    # Most CGAL doc_html tarballs have one top-level directory -> strip it
    safe_extract_xz(archive_path, root, strip_components=1)

    if not manual_index.is_file():
        warn("Manual/index.html not found after extraction; tree layout might differ.")
    return root


# -------------------------
# Filesystem comparison
# -------------------------

def is_valid_doc_relpath(rel_path: Path) -> bool:
    """
    Keep only paths like:
      Manual/...
      <PackageName>/...  (PackageName starts with uppercase)
    """
    if len(rel_path.parts) < 2:
        return False
    top = rel_path.parts[0]
    if top == "Manual":
        return True
    return top[:1].isupper()


def collect_html(root: Path) -> Set[str]:
    """
    Recursively collect all .html files under root that satisfy is_valid_doc_relpath.
    Returns a set of relative POSIX paths.
    """
    rels: Set[str] = set()
    total = 0
    for p in root.rglob("*.html"):
        total += 1
        rel = p.relative_to(root)
        if is_valid_doc_relpath(rel):
            rels.add(rel.as_posix())
    info(f"Scanned {total} .html files under {root}")
    info(f"Kept {len(rels)} files (Manual/ or <PackageName>/)")
    return rels


def exists_in_new(new_root: Path, rel: str) -> bool:
    return (new_root / Path(rel)).is_file()


# -------------------------
# CLI & main
# -------------------------

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Compare two CGAL documentation versions (local paths or auto-downloaded from GitHub using requests)."
    )

    g_old = ap.add_mutually_exclusive_group(required=True)
    g_old.add_argument("--old-path", type=Path, help="Path to old (n-1) doc root (contains Manual/)")
    g_old.add_argument("--old-version", type=str, help="Old version number to download, e.g. 5.6.2")

    g_new = ap.add_mutually_exclusive_group(required=True)
    g_new.add_argument("--new-path", type=Path, help="Path to new (n) doc root (contains Manual/)")
    g_new.add_argument("--new-version", type=str, help="New version number to download, e.g. 6.0.1")

    ap.add_argument(
        "--download-dir",
        type=Path,
        default=Path.home() / "cgal-docs",
        help="Directory to store downloaded/extracted doc_html trees (default: ~/cgal-docs)",
    )
    ap.add_argument(
        "--url-template",
        type=str,
        default=DEFAULT_URL_TEMPLATE,
        help="Release asset URL template. Must contain '{version}'.",
    )
    ap.add_argument(
        "--report",
        type=Path,
        default=Path("missing_pages_report.txt"),
        help="Output report path (default: missing_pages_report.txt)",
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=max(os.cpu_count() or 8, 8),
        help="Parallelism for filesystem checks (default: number of CPUs)",
    )
    return ap.parse_args()


def resolve_doc_root(path_arg: Optional[Path], version_arg: Optional[str], download_dir: Path, url_template: str) -> Path:
    if path_arg:
        root = path_arg.resolve()
        if not root.is_dir():
            err(f"Not a directory: {root}")
            sys.exit(2)
        return root
    assert version_arg
    return ensure_doc_tree_from_version(version_arg, download_dir, url_template=url_template)


def main() -> None:
    args = parse_args()

    old_root = resolve_doc_root(args.old_path, args.old_version, args.download_dir, args.url_template)
    new_root = resolve_doc_root(args.new_path, args.new_version, args.download_dir, args.url_template)

    info("Comparing doc trees")
    print(f"  OLD: {old_root}")
    print(f"  NEW: {new_root}")

    old_rels = collect_html(old_root)

    missing: List[str] = []
    checked = 0

    # Parallel existence check (pure filesystem)
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = {ex.submit(exists_in_new, new_root, rel): rel for rel in old_rels}
        for i, fut in enumerate(as_completed(futs), 1):
            ok = fut.result()
            rel = futs[fut]
            checked += 1
            if not ok:
                missing.append(rel)
            if i % 2000 == 0:
                print(f"[CHECK] {i} files processed...")

    print(f"\n[SUMMARY] Checked {checked} files. Missing: {len(missing)}")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", encoding="utf-8") as f:
        f.write("Missing HTML pages in NEW (relative paths):\n")
        for rel in sorted(missing):
            f.write(rel + "\n")

    info(f"Report saved to {args.report}")

if __name__ == "__main__":
    main()
