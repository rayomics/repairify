
#!/usr/bin/env python3
"""Summarize HMMER hmmsearch results and calculate alignment coverage."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize HMMER searches, calculate HMM/protein coverage, "
            "and report unmatched hits separately."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Root directory containing paired .tblout/.domtblout files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Directory where summary TSV files will be written.",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="Maximum sequence E-value (default: 1e-5).",
    )
    return parser.parse_args()


def parse_tblout(path):
    """Read sequence-level hits from HMMER --tblout."""
    hits = {}

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split(maxsplit=18)

            # Standard tblout has 18 fixed fields and an optional description.
            if len(fields) < 18:
                print(
                    f"Warning: malformed tblout line "
                    f"{line_number} in {path}"
                )
                continue

            try:
                protein_id = fields[0]

                hits[protein_id] = {
                    "evalue": float(fields[4]),
                    "bitscore": float(fields[5]),
                    "description": (
                        fields[18] if len(fields) > 18 else ""
                    ),
                }

            except ValueError:
                print(
                    f"Warning: invalid numeric field on line "
                    f"{line_number} in {path}"
                )

    return hits


def parse_domtblout(path, hmm_id):
    """Read domain coordinates, grouped by protein and HMM filename ID."""
    domains = defaultdict(list)
    lengths = {}

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split(maxsplit=22)

            # Standard domtblout has 22 fixed fields and an optional description.
            if len(fields) < 22:
                print(
                    f"Warning: malformed domtblout line "
                    f"{line_number} in {path}"
                )
                continue

            try:
                protein_id = fields[0]
                protein_length = int(fields[2])
                hmm_length = int(fields[5])

                # HMM alignment coordinates: fields 15-16 (1-based indexes).
                hmm_start = int(fields[15])
                hmm_end = int(fields[16])

                # Protein alignment coordinates: fields 17-18 (1-based indexes).
                protein_start = int(fields[17])
                protein_end = int(fields[18])

                key = (protein_id, hmm_id)

                lengths[key] = (protein_length, hmm_length)

                domains[key].append(
                    (
                        hmm_start,
                        hmm_end,
                        protein_start,
                        protein_end,
                    )
                )

            except (ValueError, IndexError):
                print(
                    f"Warning: invalid domain coordinates on line "
                    f"{line_number} in {path}"
                )

    return domains, lengths


def union_length(intervals):
    """Return the total length of merged, inclusive coordinate intervals."""
    if not intervals:
        return 0

    intervals = sorted(intervals)
    start, end = intervals[0]
    total = 0

    for next_start, next_end in intervals[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            total += end - start + 1
            start, end = next_start, next_end

    return total + end - start + 1


def identify_sample_and_hmm(path, input_dir):
    """
    Identify sample and HMM IDs.

    Nested layout:
        input/SAMPLE/HMM_ID.tblout

    Flat layout:
        input/SAMPLE-HMM_ID.tblout

    For flat filenames, the final hyphen-separated component is
    treated as the HMM ID. Nested sample folders are recommended
    when sample names and HMM IDs contain hyphens.
    """
    if path.parent != input_dir:
        return path.parent.name, path.name.removesuffix(".tblout")

    stem = path.name.removesuffix(".tblout")

    if "-" not in stem:
        return None, None

    return stem.rsplit("-", 1)


def main():
    args = parse_args()

    if not args.input.is_dir():
        raise SystemExit(
            f"Input directory does not exist: {args.input}"
        )

    if args.evalue <= 0:
        raise SystemExit("--evalue must be greater than zero.")

    files = sorted(args.input.rglob("*.tblout"))

    if not files:
        raise SystemExit(
            f"No .tblout files found in {args.input}"
        )

    # Avoid accidentally processing output files as inputs if the
    # output directory is inside the input directory.
    output_resolved = args.output.resolve()
    files = [
        path for path in files
        if output_resolved not in path.resolve().parents
        and path.resolve() != output_resolved
    ]

    if not files:
        raise SystemExit("No input .tblout files remain after excluding output.")

    args.output.mkdir(parents=True, exist_ok=True)

    hits_path = args.output / "hits.tsv"
    no_hits_path = args.output / "no_hits.tsv"
    unmatched_path = args.output / "unmatched_hits.tsv"

    hit_columns = [
        "sample",
        "hmm_id",
        "protein_id",
        "evalue",
        "bitscore",
        "hmm_length",
        "hmm_covered_positions",
        "hmm_coverage_pct",
        "protein_length",
        "protein_covered_positions",
        "protein_coverage_pct",
        "description",
    ]

    no_hit_columns = [
        "sample",
        "hmm_id",
        "evalue_threshold",
        "status",
    ]

    unmatched_columns = [
        "sample",
        "hmm_id",
        "protein_id",
        "evalue",
        "bitscore",
        "reason",
        "tblout_file",
        "domtblout_file",
    ]

    total_files = 0
    total_sequence_hits = 0
    total_significant_hits = 0
    total_unmatched = 0
    missing_domtblout = 0
    malformed_or_unidentified = 0

    with (
        hits_path.open("w", newline="") as hits_handle,
        no_hits_path.open("w", newline="") as no_hits_handle,
        unmatched_path.open("w", newline="") as unmatched_handle,
    ):
        hit_writer = csv.DictWriter(
            hits_handle,
            fieldnames=hit_columns,
            delimiter="\t",
        )
        no_hit_writer = csv.DictWriter(
            no_hits_handle,
            fieldnames=no_hit_columns,
            delimiter="\t",
        )
        unmatched_writer = csv.DictWriter(
            unmatched_handle,
            fieldnames=unmatched_columns,
            delimiter="\t",
        )

        hit_writer.writeheader()
        no_hit_writer.writeheader()
        unmatched_writer.writeheader()

        for path in files:
            sample, hmm_id = identify_sample_and_hmm(
                path, args.input
            )

            if sample is None or hmm_id is None:
                print(
                    f"Warning: cannot identify sample/HMM from filename: {path}"
                )
                malformed_or_unidentified += 1
                continue

            dom_path = path.with_suffix(".domtblout")

            sequence_hits = parse_tblout(path)
            total_files += 1
            total_sequence_hits += len(sequence_hits)

            # Only sequence-level hits passing the E-value threshold
            # are considered significant for this summary.
            eligible_hits = {
                protein: hit
                for protein, hit in sequence_hits.items()
                if hit["evalue"] < args.evalue
            }

            total_significant_hits += len(eligible_hits)

            # A genuine no-hit is recorded only if there are no
            # sequence-level hits passing the threshold.
            if not eligible_hits:
                no_hit_writer.writerow(
                    {
                        "sample": sample,
                        "hmm_id": hmm_id,
                        "evalue_threshold": args.evalue,
                        "status": "no_sequence_hits_passing_evalue",
                    }
                )
                continue

            # If domain output is missing, do not call the search a no-hit.
            if not dom_path.is_file():
                print(f"Warning: missing paired file: {dom_path}")
                missing_domtblout += 1

                for protein, sequence in eligible_hits.items():
                    unmatched_writer.writerow(
                        {
                            "sample": sample,
                            "hmm_id": hmm_id,
                            "protein_id": protein,
                            "evalue": sequence["evalue"],
                            "bitscore": sequence["bitscore"],
                            "reason": "missing_domtblout_file",
                            "tblout_file": str(path),
                            "domtblout_file": str(dom_path),
                        }
                    )
                    total_unmatched += 1

                continue

            domains, lengths = parse_domtblout(dom_path, hmm_id)

            for protein, sequence in eligible_hits.items():
                key = (protein, hmm_id)

                if key not in domains or key not in lengths:
                    print(
                        f"Warning: no domain coordinates for protein "
                        f"{protein} in {dom_path}"
                    )

                    unmatched_writer.writerow(
                        {
                            "sample": sample,
                            "hmm_id": hmm_id,
                            "protein_id": protein,
                            "evalue": sequence["evalue"],
                            "bitscore": sequence["bitscore"],
                            "reason": "no_matching_domain_record",
                            "tblout_file": str(path),
                            "domtblout_file": str(dom_path),
                        }
                    )
                    total_unmatched += 1
                    continue

                protein_length, hmm_length = lengths[key]
                intervals = domains[key]

                hmm_positions = union_length(
                    [(domain[0], domain[1]) for domain in intervals]
                )

                protein_positions = union_length(
                    [(domain[2], domain[3]) for domain in intervals]
                )

                hmm_coverage = (
                    100.0 * hmm_positions / hmm_length
                    if hmm_length > 0 else 0.0
                )

                protein_coverage = (
                    100.0 * protein_positions / protein_length
                    if protein_length > 0 else 0.0
                )

                hit_writer.writerow(
                    {
                        "sample": sample,
                        "hmm_id": hmm_id,
                        "protein_id": protein,
                        "evalue": sequence["evalue"],
                        "bitscore": sequence["bitscore"],
                        "hmm_length": hmm_length,
                        "hmm_covered_positions": hmm_positions,
                        "hmm_coverage_pct": round(hmm_coverage, 2),
                        "protein_length": protein_length,
                        "protein_covered_positions": protein_positions,
                        "protein_coverage_pct": round(
                            protein_coverage, 2
                        ),
                        "description": sequence["description"],
                    }
                )

    print("\nHMMER summary complete.")
    print(f"Processed result files: {total_files}")
    print(f"Sequence hits parsed: {total_sequence_hits}")
    print(
        f"Sequence hits passing E-value < {args.evalue}: "
        f"{total_significant_hits}"
    )
    print(f"Hits with coverage calculated: {total_significant_hits - total_unmatched}")
    print(f"Hits without matched domain data: {total_unmatched}")
    print(f"Missing .domtblout files: {missing_domtblout}")
    print(f"Files with unidentified sample/HMM: {malformed_or_unidentified}")
    print(f"\nHits: {hits_path}")
    print(f"No sequence hits passing threshold: {no_hits_path}")
    print(f"Unmatched hits requiring review: {unmatched_path}")


if __name__ == "__main__":
    main()

