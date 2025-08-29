#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
from gtfparse import read_gtf

novel_loci_categories = {"antisense", "genic_intron", "intergenic"}


def filter_gtf(gtf_df, transcript_ids):
    return gtf_df[gtf_df["transcript_id"].isin(transcript_ids)].copy()


def write_gtf(df, output_path):
    # 属性フィールドを整形
    df['frame'] = "."
    df['score'] = "."
    df["attribute"] = (
        'transcript_id "' + df["transcript_id"].astype(str) + '"; ' +
        'gene_id "' + df["gene_id"].astype(str) + '";'
    )

    df = df[["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]]
    df.to_csv(output_path, sep="\t", header=False, index=False, quoting=3)


def main():
    parser = argparse.ArgumentParser(description="Filter isoforms and save based on category.")
    parser.add_argument("--previous_sqanti_res", type=Path, required=True, help="sqantiによる分類結果")
    parser.add_argument("--previous_ref", type=Path, required=True, help="前のsqantiラウンドでrefとしたトランスクリプトgtf")
    parser.add_argument("--all_tx_gtf", type=Path, required=True, help="Original GTF file")
    parser.add_argument("--next_target", type=Path, required=True, help="次のsqantiラウンドでtargetとなるトランスクリプトgtf")
    parser.add_argument("--next_ref", type=Path, required=True, help="次のsqantiラウンドでrefとなるトランスクリプトgtf")
    args = parser.parse_args()

    df = pd.read_csv(args.previous_sqanti_res, sep="\t")
    gtf_df = read_gtf(args.all_tx_gtf).to_pandas()


    # このsqantiラウンドでtargetとなるトランスクリプト（= 前のラウンドのnew loci tx）
    target_df = df[df["structural_category"].isin(novel_loci_categories)]
    target_ids = set(target_df["isoform"])
    target_gtf = filter_gtf(gtf_df, target_ids)
    write_gtf(target_gtf, args.next_target)


    # このsqantiラウンドで新たにrefとなるトランスクリプト（= 前のラウンドのnew isoform tx）
    ref_df = df[
        ( df["structural_category"] != "full-splice_match") &
        (~df["structural_category"].isin(novel_loci_categories)) &
        (~df["subcategory"].str.contains("intron_retention", na=False)) &
        ( df["structural_category"] != "fusion") &
        (
            (df["structural_category"] != "genic") |
            (df["subcategory"] != "mono-exon")
        )
    ]
    ref_ids = set(ref_df["isoform"])


    # 前のラウンドのref txはref側のベースとして追加（重複は自動で除かれる）
    previous_ref_gtf_df = read_gtf(args.previous_ref).to_pandas()
    previous_ref_ids = set(previous_ref_gtf_df["transcript_id"].dropna().unique())
    ref_ids.update(previous_ref_ids)

    ref_gtf = filter_gtf(gtf_df, ref_ids)
    write_gtf(ref_gtf, args.next_ref)


if __name__ == "__main__":
    main()
