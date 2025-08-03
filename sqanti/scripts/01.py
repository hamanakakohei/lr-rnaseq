#!/usr/bin/env python3
from gtfparse import read_gtf
import pandas as pd
import polars as pl
import pybedtools
from pybedtools import BedTool
import networkx as nx


gtf_path = '/path/to/sqanti/SQANTI3/data/reference/gencode.v38.basic_chr22.gtf'
out0 =     "/path/to/sqanti/practice/results/01/all_gene_tx.tsv"
out1 =     "/path/to/sqanti/practice/results/01/overlapping_transcript_groups.tsv"
out2 =     "/path/to/sqanti/practice/results/01/non-overlap_transcripts.tsv"


# 1. gtf を df として読み込み
df = read_gtf( gtf_path )
df\
  [['gene_id','transcript_id']].\
  filter(pl.col("transcript_id") != "").\
  unique().\
  write_csv( out0, separator="\t" )

exons = df.\
  filter(df["feature"] == "exon")\
  [["seqname", "start", "end", "strand", "transcript_id"]].\
  to_pandas()


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

pd.DataFrame(grouped_list).\
  to_csv( out1, sep="\t", index=False )


# 8. グループに含まれていない transcript_id を取得して保存
all_txs = set(exons["transcript_id"])
overlapping_txs = set(d["transcript_id"] for d in grouped_list)
non_overlap_txs = all_txs - overlapping_txs

pd.DataFrame({"transcript_id": sorted( non_overlap_txs )}).\
  to_csv( out2, sep="\t", index=False )


