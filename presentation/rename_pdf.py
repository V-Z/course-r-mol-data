#!/usr/bin/env python3
"""Rename a rendered book PDF.

Quarto names the merged book PDF after the book title (e.g.
docs/Molecular-data-in-R.pdf) and this cannot be overridden from YML,
so the beamer/A4 profiles rename it in post-render.

Usage: python3 rename_pdf.py SRC DST   (paths relative to the project root)
"""
import os
import sys


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: python3 rename_pdf.py SRC DST")
    src, dst = sys.argv[1], sys.argv[2]
    if os.path.exists(src):
        os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
        os.replace(src, dst)
        print("rename_pdf.py: %s -> %s" % (src, dst))
    else:
        print("rename_pdf.py: %s not found, nothing to rename" % src)


if __name__ == "__main__":
    main()
