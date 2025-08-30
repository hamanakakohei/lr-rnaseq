#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--round1_class", required=True, help="Path to old isoform_classification.txt (round 0)")
    parser.add_argument("--round2_class", required=True, help="Path to new isoform_classification.txt (round 1)")
    parser.add_argument("--overlap_anno", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # 旧ラウンド isoform_classification.txt
    old_df = pd.read_csv(args.round1_class, sep="\t", usecols=["isoform", "associated_gene"])
    old_df = old_df.rename(columns={"isoform": "transcript_id", "associated_gene": "gene_id_old"})

    # 新ラウンド isoform_classification.txt
    new_df = pd.read_csv(args.round2_class, sep="\t", usecols=["isoform", "associated_gene"])
    new_df = new_df.rename(columns={"isoform": "transcript_id", "associated_gene": "gene_id"})

    # overlap_txs.tsv
    overlap_df = pd.read_csv(args.overlap_anno, sep="\t")
    overlap_df = overlap_df.rename(columns={"group_id": "gene_id"})

    # 新ラウンドと overlap を縦に結合
    # transcript_idごとに最初(new)→最後(overlap)の順序を保持
    # transcript_id ごとに gene_id を決定（overlap優先）
    combined = pd.concat([new_df, overlap_df], ignore_index=True)
    combined = combined.drop_duplicates(subset=["transcript_id"], keep="last")
    combined = combined.rename(columns={"gene_id": "gene_id_new"})

    # old の gene_id と結合
    result = pd.merge(old_df, combined, on="transcript_id", how="outer")

    # 出力
    result = result[["transcript_id", "gene_id_old", "gene_id_new"]]
    result.to_csv(args.out, sep="\t", index=False, na_rep="NA")

if __name__ == "__main__":
    main()
