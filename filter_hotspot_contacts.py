#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Filter RFD3 designs based on hotspot-contact summary CSV."
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Path to hotspot_contacts_summary.csv"
    )

    parser.add_argument(
        "--cif-dir",
        required=True,
        help="Directory containing original .cif.gz files"
    )

    parser.add_argument(
        "--out-dir",
        required=True,
        help="Directory where passing .cif.gz files will be copied"
    )

    parser.add_argument(
        "--min-contacts",
        type=int,
        default=1,
        help="Minimum number of hotspot atom contacts required. Default: 1"
    )

    parser.add_argument(
        "--max-distance",
        type=float,
        default=4.0,
        help="Maximum minimum hotspot distance in Angstroms. Default: 4.0"
    )

    parser.add_argument(
        "--copy",
        action="store_true",
        help="Copy passing files instead of symlinking them."
    )

    args = parser.parse_args()

    summary_path = Path(args.summary)
    cif_dir = Path(args.cif_dir)
    out_dir = Path(args.out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(summary_path)

    print(f"Loaded {len(df)} designs from {summary_path}")

    # Make sure numeric columns are numeric
    df["n_contacts"] = pd.to_numeric(df["n_contacts"], errors="coerce").fillna(0)
    df["min_distance_A"] = pd.to_numeric(df["min_distance_A"], errors="coerce")

    passing = df[
        (df["n_contacts"] >= args.min_contacts) &
        (df["min_distance_A"] <= args.max_distance)
    ].copy()

    passing = passing.sort_values(
        ["n_contacts", "min_distance_A"],
        ascending=[False, True]
    )

    passing_csv = out_dir / "passing_hotspot_contacts_summary.csv"
    passing.to_csv(passing_csv, index=False)

    print(f"Passing designs: {len(passing)}")
    print(f"Wrote passing summary to: {passing_csv}")

    copied = 0
    missing = 0

    for _, row in passing.iterrows():
        filename = row["file"]
        src = cif_dir / filename
        dst = out_dir / filename

        if not src.exists():
            print(f"Missing file: {src}")
            missing += 1
            continue

        if dst.exists():
            continue

        if args.copy:
            shutil.copy2(src, dst)
        else:
            dst.symlink_to(src.resolve())

        copied += 1

    print(f"Linked/copied passing structures: {copied}")
    print(f"Missing structures: {missing}")

    print("\nTop passing designs:")
    print(
        passing[
            ["file", "binder_chain", "n_contacts", "min_distance_A"]
        ].head(20).to_string(index=False)
    )


if __name__ == "__main__":
    main()