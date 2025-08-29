#!/usr/bin/env python3

import argparse
from gtfparse import read_gtf
import pandas as pd
import polars as pl
import pybedtools
from pybedtools import BedTool
import networkx as nx
from pathlib import Path


parser = argparse.ArgumentParser(description="Detect overlapping transcript groups from a GTF file.")
parser.add_argument("--new_loci_gtf", type=Path, required=True, help="new loci txs GTF")
parser.add_argument("--all_tx_gtf", type=Path, required=True, help="位置情報を付けるための全てのtxのGTF")
parser.add_argument("--out_overlap_txs", type=Path, required=True, help="")
parser.add_argument("--out_overlap_txs_bed", type=Path, required=True, help="")
parser.add_argument("--out_nonoverlap_txs", type=Path, required=True, help="")
parser.add_argument("--out_nonoverlap_txs_bed", type=Path, required=True, help="")
args = parser.parse_args()


# 1. gtf を df として読み込み
df = read_gtf( args.new_loci_gtf )

exons = df.\
  filter(df["feature"] == "exon")\
  [["seqname", "start", "end", "strand", "transcript_id"]].\
  to_pandas()

df_all_txs = read_gtf( args.all_tx_gtf ).\
  to_pandas().\
  query("feature == 'transcript'")\
  [["transcript_id", "seqname", "start", "end"]]


# 2. BEDに変換
exons["chrom"] = exons["seqname"]
exons["name"] = exons["transcript_id"]
exons["score"] = 0
exons["start"] = exons["start"] - 1
bed_df = exons[["chrom", "start", "end", "name", "score", "strand"]]


# 3. bedtoolsでエクソンがかぶるtranscriptを見つける（-sでstrand一致、-wa -wbでペア取得）
bed = BedTool.from_dataframe(bed_df)
intersected = bed.intersect(bed, s=True, wa=True, wb=True)
overlap_df = intersected.to_dataframe(names=[
  "chrom_a", "start_a", "end_a", "name_a", "score_a", "strand_a",
  "chrom_b", "start_b", "end_b", "name_b", "score_b", "strand_b"
])


# 4. 自分自身とのかぶりを除く
filtered = overlap_df[
    (overlap_df["name_a"] != overlap_df["name_b"])
][["name_a", "name_b"]].drop_duplicates()


# 6. グループ化
G = nx.Graph()
G.add_edges_from(filtered.values)
groups = list(nx.connected_components(G))


# 7. グループをファイルに保存
grouped_list = []

for i, group in enumerate(groups, 1):
  for tx in sorted(group):
    grouped_list.append({"group_id": f"Group{i}", "transcript_id": tx})

df_groups = pd.DataFrame(grouped_list)
df_groups.to_csv( args.out_overlap_txs, sep="\t", index=False )

df_groups = pd.merge(df_groups, df_all_txs, on="transcript_id", how="left")
df_groups["name"] = (
    df_groups["group_id"] + ";" +
    df_groups["transcript_id"]
)
df_groups["score"] = "."
df_groups["strand"] = "."
df_groups[["seqname", "start", "end", "name", "score", "strand"]].to_csv(
    args.out_overlap_txs_bed, sep="\t", index=False, header=False
)


# 8. グループに含まれていない transcript_id を取得して保存
all_txs = set(exons["transcript_id"])
overlapping_txs = set(d["transcript_id"] for d in grouped_list)
non_overlap_txs = all_txs - overlapping_txs

df_noOverlap = pd.DataFrame({"transcript_id": sorted( non_overlap_txs )})
df_noOverlap.to_csv( args.out_nonoverlap_txs, sep="\t", index=False )

df_noOverlap = pd.merge(df_noOverlap, df_all_txs, on="transcript_id", how="left")
df_noOverlap["name"] = df_noOverlap["transcript_id"]
df_noOverlap["score"] = "."
df_noOverlap["strand"] = "."
df_noOverlap[["seqname", "start", "end", "name", "score", "strand"]].to_csv(
    args.out_nonoverlap_txs_bed, sep="\t", index=False, header=False
)
