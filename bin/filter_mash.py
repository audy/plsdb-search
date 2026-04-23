#!/usr/bin/env python3
"""
Filter mash screen output by identity and p-value thresholds.

Mash screen output columns (tab-separated):
  identity  shared-hashes  median-multiplicity  p-value  query-ID  query-comment
"""

import argparse
import sys


def parse_accession(query_id):
    """Extract accession from mash query-ID (may include version like NZ_CP123.1)."""
    return query_id.strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",        required=True)
    parser.add_argument("--min_identity", type=float, default=0.90)
    parser.add_argument("--max_pvalue",   type=float, default=1e-5)
    parser.add_argument("--output",       required=True)
    parser.add_argument("--accessions",   required=True)
    args = parser.parse_args()

    header = ["identity", "shared_hashes", "median_multiplicity", "pvalue",
              "query_id", "query_comment"]

    kept = []
    with open(args.input) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            identity   = float(parts[0])
            pvalue     = float(parts[3])
            query_id   = parts[4]
            comment    = parts[5] if len(parts) > 5 else ""

            if identity >= args.min_identity and pvalue <= args.max_pvalue:
                kept.append({
                    "identity":            parts[0],
                    "shared_hashes":       parts[1],
                    "median_multiplicity": parts[2],
                    "pvalue":              parts[3],
                    "query_id":            query_id,
                    "query_comment":       comment,
                })

    # Sort by identity descending
    kept.sort(key=lambda r: float(r["identity"]), reverse=True)

    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in kept:
            out.write("\t".join(row[h] for h in header) + "\n")

    with open(args.accessions, "w") as out:
        for row in kept:
            out.write(parse_accession(row["query_id"]) + "\n")

    print(f"Kept {len(kept)} hits (identity >= {args.min_identity}, "
          f"p-value <= {args.max_pvalue})", file=sys.stderr)


if __name__ == "__main__":
    main()
