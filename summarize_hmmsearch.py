
#!/usr/bin/env python3
"""Summarize HMMER hmmsearch --tblout results."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize DNA-repair HMM search results."
    )
    parser.add_argument(
        "--input", required=True, type=Path,
        help="Directory containing hmmsearch *.tblout files."
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="Output directory for summary tables."
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="Maximum sequence E-value (default: 1e-5)."
    )
    return parser.parse_args()


def parse_tblout(path):
    """Yield parsed hits from an HMMER --tblout file."""
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split(maxsplit=18)

            # HMMER --tblout has 18 fixed fields followed by
            # the free-text description (which may be absent).
            if len(fields) < 18:
                print(f"Warning: skipping malformed line in {path}")
                continue

            try:
                yield {
                    "target": fields[0],
                    "target_accession": fields[1],
                    "query": fields[2],
                    "query_accession": fields[3],
                    "evalue": float(fields[4]),
                    "score": float(fields[5]),
                    "description": fields[18] if len(fields) > 18 else "",
                }
            except ValueError:
                print(f"Warning: invalid numeric field in {path}")


def main():
    args = parse_args()

    if not args.input.is_dir():
        raise SystemExit(f"Input directory does not exist: {args.input}")

    if args.evalue <= 0:
        raise SystemExit("--evalue must be greater than zero.")

    files = sorted(args.input.rglob("*.tblout"))

    if not files:
        raise SystemExit(f"No .tblout files found in {args.input}")

    args.output.mkdir(parents=True, exist_ok=True)

    hits_path = args.output / "hits.tsv"
    no_hits_path = args.output / "no_hits.tsv"

    hit_columns = [
        "sample", "hmm_id", "protein_id",
        "evalue", "bitscore", "description"
    ]

    no_hit_columns = ["sample", "hmm_id"]

    total_hits = 0
    total_files = 0

    with hits_path.open("w", newline="") as hits_handle, \
         no_hits_path.open("w", newline="") as no_hits_handle:

        hit_writer = csv.DictWriter(
            hits_handle, fieldnames=hit_columns, delimiter="\t"
        )
        no_hit_writer = csv.DictWriter(
            no_hits_handle, fieldnames=no_hit_columns, delimiter="\t"
        )

        hit_writer.writeheader()
        no_hit_writer.writeheader()

        for path in files:
            # Expected layout:
            # input/SAMPLE/HMM_ID.tblout
            # Also supports flat directories using filenames
            # SAMPLE-HMM_ID.tblout as a fallback.
            if path.parent != args.input:
                sample = path.parent.name
                hmm_id = path.name.removesuffix(".tblout")
            else:
                stem = path.name.removesuffix(".tblout")
                if "-" not in stem:
                    print(f"Warning: cannot identify sample/HMM: {path}")
                    continue
                sample, hmm_id = stem.rsplit("-", 1)

            total_files += 1
            accepted = 0

            for hit in parse_tblout(path):
                if hit["evalue"] >= args.evalue:
                    continue

                hit_writer.writerow({
                    "sample": sample,
                    "hmm_id": hmm_id,
                    "protein_id": hit["target"],
                    "evalue": hit["evalue"],
                    "bitscore": hit["score"],
                    "description": hit["description"],
                })
                accepted += 1
                total_hits += 1

            if accepted == 0:
                no_hit_writer.writerow({
                    "sample": sample,
                    "hmm_id": hmm_id,
                })

    print(f"Processed HMM result files: {total_files}")
    print(f"Significant protein-HMM matches: {total_hits}")
    print(f"Hits written to: {hits_path}")
    print(f"No-hit combinations written to: {no_hits_path}")


if __name__ == "__main__":
    main()
