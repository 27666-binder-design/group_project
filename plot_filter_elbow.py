#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


# Folder containing passing_hotspot_contacts_1, ..., passing_hotspot_contacts_6
base_dir = Path("work/exp_04")

filter_nums = range(1, 8)

counts = []

for num in filter_nums:
    summary_csv = (
        base_dir
        / f"passing_hotspot_contacts_{num}"
        / "passing_hotspot_contacts_summary.csv"
    )

    if not summary_csv.exists():
        print(f"Missing: {summary_csv}")
        counts.append(0)
        continue

    # Count number of rows excluding header
    df = pd.read_csv(summary_csv)
    n_designs = len(df)

    counts.append(n_designs)
    print(f"Filter {num}: {n_designs} passing designs")


# Make plot
plt.figure(figsize=(6, 4))
plt.plot(list(filter_nums), counts, marker="o", linewidth=2)

plt.xlabel("Number of hotspot atom contacts")
plt.ylabel("Number of passing designs")
plt.title("Passing designs after hotspot filtering")
plt.xticks(list(filter_nums))
plt.grid(True, alpha=0.3)

plt.tight_layout()

out_png = base_dir / "hotspot_filter_elbow.png"
plt.savefig(out_png, dpi=300)

print(f"\nSaved plot to: {out_png}")