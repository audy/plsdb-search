#!/usr/bin/env python3
"""
Filter samtools coverage output by minimum breadth and mean depth.

samtools coverage columns:
  #rname  startpos  endpos  numreads  covbases  coverage  meandepth  meanbaseq  meanmapq
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",            required=True)
    parser.add_argument("--min_coverage_pct", type=float, default=50.0)
    parser.add_argument("--min_mean_depth",   type=float, default=1.0)
    parser.add_argument("--output",           required=True)
    args = parser.parse_args()

    kept = 0
    total = 0

    with open(args.input) as fh, open(args.output, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            total += 1
            coverage_pct = float(parts[5])
            mean_depth   = float(parts[6])
            if coverage_pct >= args.min_coverage_pct and mean_depth >= args.min_mean_depth:
                out.write(line)
                kept += 1

    print(f"Kept {kept}/{total} plasmids "
          f"(coverage >= {args.min_coverage_pct}%, depth >= {args.min_mean_depth}x)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
