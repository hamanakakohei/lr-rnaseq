#!/usr/bin/env python3
import argparse
import pandas as pd
from gtfparse import read_gtf


def main():
    parser = argparse.ArgumentParser( description="新旧のsqanti分類表を比べる" )
    parser.add_argument("--old", required=True, help="")
    parser.add_argument("--new", required=True, help="")
    parser.add_argument("--gtf", required=True, help="txの位置情報を付けるためのGTF")
    parser.add_argument("--out", required=True, help="")
    parser.add_argument("--out_bed", required=True, help="")
    args = parser.parse_args()

    # old/new の表を読み込み結合する
    usecols = ["isoform", "structural_category", "subcategory"]
    old = pd.read_csv(args.old, sep="\t", usecols=usecols)
    new = pd.read_csv(args.new, sep="\t", usecols=usecols)

    old = old.rename(
        columns={
            "structural_category": "structural_category_old",
            "subcategory": "subcategory_old",
        }
    )
    new = new.rename(
        columns={
            "structural_category": "structural_category_new",
            "subcategory": "subcategory_new",
        }
    )

    new_old_class = pd.merge(old, new, on="isoform", how="inner")
    new_old_class.to_csv(args.out, sep="\t", index=False)

    # GTF の読み込み
    gtf = read_gtf(args.gtf).to_pandas()
    tx_pos = gtf.query("feature == 'transcript'")[["transcript_id", "seqname", "start", "end"]]

    # 分類表にgtfで位置情報を加える
    final = pd.merge(
        new_old_class,
        tx_pos[["transcript_id", "seqname", "start", "end"]],
        left_on="isoform",
        right_on="transcript_id",
        how="left"
    ).drop(columns=["transcript_id"])

    # 変化した行だけ残す
    changed = final[
        (final["structural_category_old"] != final["structural_category_new"]) |
        (final["subcategory_old"] != final["subcategory_new"])
    ].copy()

    # BEDとして保存する
    changed["name"] = (
        changed["isoform"] + ";" +
        changed["structural_category_old"] + ";" +
        changed["subcategory_old"] + ";" +
        changed["structural_category_new"] + ";" +
        changed["subcategory_new"]
    )
    changed["score"] = "."
    changed["strand"] = "."
    changed[["seqname", "start", "end", "name", "score", "strand"]].to_csv(
        args.out_bed, sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    main()
