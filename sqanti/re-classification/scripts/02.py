#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path

def merge_classification_results(base_dir: Path, n_rounds: int, filename: str) -> pd.DataFrame:
    dfs = []
    for i in range(1, n_rounds + 1):
        f = base_dir / str(i) / filename
        df = pd.read_csv(f, sep="\t")

        # old → i-1, new → i にリネーム
        df = df.rename(columns={
            "structural_category_old": f"structural_category_{i-1}",
            "subcategory_old": f"subcategory_{i-1}",
            "structural_category_new": f"structural_category_{i}",
            "subcategory_new": f"subcategory_{i}",
        })
        dfs.append(df)

    # isoform をキーに順に outer merge
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, how="outer")

    return merged


def main():
    parser = argparse.ArgumentParser(description="Merge classification results across rounds.")
    parser.add_argument("--base_dir", type=Path, required=True,
                        help="Parent directory containing round subdirs (1/,2/,3/...)")
    parser.add_argument("--n_rounds", type=int, required=True)
    parser.add_argument("--filename", type=str, default="new_old_classification_comparison.txt")
    parser.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()

    merged = merge_classification_results(args.base_dir, args.n_rounds, args.filename)
    merged.to_csv(args.output, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
