#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import itertools
import networkx as nx


# 入力ファイル（tx1, tx2ペアが入っている）
PAIR_FILE =                 '/path/to/sqanti/practice/results/02/pairs.txt'
non_overlap_tx_path =       '/path/to/sqanti/practice/results/01/non-overlap_transcripts.tsv'
truth_set_of_gene_tx_path = '/path/to/sqanti/practice/results/01/all_gene_tx.tsv'


# 出力結果をまとめるリスト
# 読み込む列番号（0-based index）
results = []
col_indices = [0, 5, 6, 14]

with open(PAIR_FILE) as f:
  _ = next(f)  # ヘッダー行はスキップ
  #for line in f:
  for line in itertools.islice(f, 8000):
    _, tx1, tx2 = line.strip().split()
    #
    path = Path(f"/path/to/sqanti/practice/results/03/{tx1}-{tx2}/tx1_tx2_classification.txt")
    try:
      with open(path) as f2:
        _ = next(f2)  # ヘッダー行はスキップ
        second_line = next(f2).strip().split('\t')
      #
      selected = [second_line[i] for i in col_indices] # 必要な列だけ取り出す
      results.append([tx1, tx2] + selected)
    #
    except Exception as e:
      print(f"Error reading {path}: {e}")


# まとめてDataFrameにする
columns = ["tx1", "tx2", "isoform", "structural_category", "associated_gene", "subcategory"]
df = pd.DataFrame(results, columns=columns)



# 1. tx1, tx2 をノードとして、同じグループに属するようなペアを集める（グラフ的アプローチ）
# 条件を満たす行のみでエッジを作成
# 対象となる structural_category のセット
same_gene_categories = {'full-splice_match', 'novel_in_catalog', 'incomplete-splice_match', 'novel_not_in_catalog'} #, 'genic'}

G = nx.Graph()
for _, row in df.iterrows():
  if row['structural_category'] in same_gene_categories:
    G.add_edge(row['tx1'], row['tx2'])


# 各ノードにグループIDを割り当てる
group_map = {}
for i, component in enumerate(nx.connected_components(G)):
  for tx in component:
    group_map[tx] = i

df['group'] = df['tx1'].map(group_map).fillna(-1).astype(int)


# グループ分けの正しさを評価する
## non_overlap_txを加えつつ、group==-1も含めてgroup idを振りなおす
#overlapping_tx = df.rename({'isoform': 'transcript_id'}, axis=1)[['transcript_id']].drop_duplicates()
non_overlap_tx = pd.read_table( non_overlap_tx_path )
truth = pd.read_table( truth_set_of_gene_tx_path )

our = (
  pd.concat([
    df.rename({'isoform': 'transcript_id'}, axis=1),
    non_overlap_tx
  ])
  [['transcript_id', 'group']]
  .fillna(-1)
  .astype({'group': int})
  .drop_duplicates()
  .reset_index(drop=True)
)

# -1 の行に対して、新しい group ID を付け直す
max_group = our['group'].max()
unassigned_mask = our['group'] == -1
n_unassigned = unassigned_mask.sum()

# -1 のところに、max_group+1 から順番に振っていく
our.loc[unassigned_mask, 'group'] = range(max_group + 1, max_group + 1 + n_unassigned)

merged_df = pd.merge(our, truth, on="transcript_id")
merged_df['group_size'] = merged_df.groupby('group')['group'].transform('count')
merged_df['gene_id_size'] = merged_df.groupby('gene_id')['gene_id'].transform('count')


# B列グループの集合（Aのリスト）
# C列グループの集合（Aのリスト）
truth_groups = merged_df.groupby("gene_id")["transcript_id"].apply(lambda x: frozenset(x)).to_dict()
our_groups = merged_df.groupby("group"  )["transcript_id"].apply(lambda x: frozenset(x)).to_dict()

# グループ数
num_truth_groups = len(truth_groups)
num_our_groups = len(our_groups)

# 完全一致しているグループ数（frozensetで集合の中身が完全一致するものを数える）
intersection = set(truth_groups.values()) & set(our_groups.values())
num_common_groups = len(intersection)

# 出力
print(f"B列のグループ数: {num_b_groups}")
print(f"C列のグループ数: {num_c_groups}")
print(f"BとCで完全一致するグループ数: {num_common_groups}")






# ラベルを割り当てる
def label_tx1(subdf):
  # 優先順位①: intron_retention が一度でもある
  if (subdf['subcategory'] == 'intron_retention').any():
    return 'intron_retention'
  #
  # 優先順位②: target categories が一度も登場せず、かつ genic + mono-exon の行が存在する
  target_categories = {'novel_in_catalog', 'incomplete-splice_match', 'novel_not_in_catalog'}
  if (
    (~subdf['structural_category'].isin(target_categories)).all() and
    ((subdf['structural_category'] == 'genic') & (subdf['subcategory'] == 'mono-exon')).any()
  ):
    return 'genic_mono-exon'
  #
  # それ以外
  return 'primary'


tx1_labels = df.groupby('tx1').apply(label_tx1).reset_index()
tx1_labels.columns = ['tx1', 'tx1_label']

df = df.merge(tx1_labels, on='tx1', how='left')


