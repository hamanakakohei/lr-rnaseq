#!/usr/bin/env python3
import pandas as pd
from itertools import permutations
from itertools import combinations
from gtfparse import read_gtf
import polars as pl
import os


# 1. tx1 - tx2 ペアのリストを作って保存する
group_transcript_path = '/path/to/sqanti/practice/results/01/overlapping_transcript_groups.tsv'
out_pairs =             '/path/to/sqanti/practice/results/02/pairs.txt'
out_trios =             '/path/to/sqanti/practice/results/02/trios.txt'

group_transcript_df = pd.read_table( group_transcript_path )

pairs = []
for group, group_df in group_transcript_df.groupby('group_id'):
  txs = group_df['transcript_id'].tolist()
  for tx1, tx2 in permutations(txs, 2):
    if tx1 != tx2:
      pairs.append({'group_id': group, 'tx1': tx1, 'tx2': tx2})

pl.DataFrame(pairs, schema=["group_id", "tx1", "tx2"]).\
  write_csv(out_pairs, separator="\t")


# 2. tx1 vs tx2/3 トリオ（tx2/3 がrefでtx1が対象）のリストを保存する
seen = set()
trios = []

for group, group_df in group_transcript_df.groupby('group_id'):
  txs = group_df['transcript_id'].tolist()
  for tx1, tx2, tx3 in permutations(txs, 3):
    if tx1 != tx2 and tx1 != tx3 and tx2 != tx3:
      key = (group, tx1, tuple(sorted([tx2, tx3])))  # tx2, tx3の順序は無視
      if key not in seen:
        seen.add(key)
        trios.append({
          'group_id': group,
          'tx1': tx1,
          'tx2': key[2][0],
          'tx3': key[2][1]
        })

pl.DataFrame(trios, schema=["group_id", "tx1", "tx2", "tx3"]).\
  write_csv(out_trios, separator="\t")


# 3. gtf を読み込んで、各transcript_idごとのgtfにして保存する
gtf_path =     '/path/to/sqanti/SQANTI3/data/reference/gencode.v38.basic_chr22.gtf'
out_gtfs_dir = '/path/to/sqanti/practice/results/02/gtfs/'

def save_transcript_gtfs( gtf_df: pl.DataFrame, transcript_ids: list[str], out_dir: str):
  os.makedirs(out_dir, exist_ok=True)
  for tx in transcript_ids:
    tx_df = gtf_df.filter(pl.col("transcript_id") == tx)
    tx_path = os.path.join(out_dir, f"{tx}.gtf")
    tx_df\
      .with_columns([
          pl.format('gene_id "{}"; transcript_id "{}";', pl.col("gene_id"), pl.col("transcript_id")).alias("attribute")
      ])\
      .select("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")\
      .to_pandas()\
      .to_csv(tx_path, sep="\t", header=False, index=False, quoting=3)

gtf_df = read_gtf( gtf_path )
save_transcript_gtfs( gtf_df, group_transcript_df['transcript_id'], out_gtfs_dir )


# 4. そのgtfを二つずつ足しつつ、gene_idを全て別にする（のでいいのか？）
