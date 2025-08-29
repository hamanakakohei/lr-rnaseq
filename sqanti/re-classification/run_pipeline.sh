#!/bin/bash

set -euo pipefail
eval "$(conda shell.bash hook)"


SQANTI_SCRIPT="./sqanti3_qc.py"
ALL_TX_GTF="data/all_txs.gtf"
FASTA="data/GRCh38.primary_assembly.genome.fa"
rounds=(1 2)


# 0：最初のsqanti結果を、おかわり0回目として準備する
mkdir -p results/01/0/
ln -s data/isoforms_classification.txt results/01/0/


# 1：sqantiを何度かおかわりする
for i in ${rounds[@]}; do
  PREVIOUS_RES=results/01/$((i-1))/isoforms_classification.txt
  PREVIOUS_REF=results/01/$((i-1))/isoforms.ref.gtf
  OUT_DIR=results/01/${i}/
  NEXT_TARGET=results/01/${i}/isoforms.target.gtf
  NEXT_REF=results/01/${i}/isoforms.ref.gtf

  mkdir -p "$OUT_DIR"

  # 1-1 前回のsqanti分類結果を使って、new isoformを取り込んだnext ref gtfとnew lociのみのnext targetをつくる
  conda activate misc_20250301
  scripts/01.py \
    --previous_sqanti_res $PREVIOUS_RES \
    --previous_ref $PREVIOUS_REF \
    --all_tx_gtf $ALL_TX_GTF \
    --next_target $NEXT_TARGET \
    --next_ref $NEXT_REF


  # 1-2 それをもとにsqantiする
  conda activate sqanti3
  "$SQANTI_SCRIPT" \
    --isoforms "$NEXT_TARGET" \
    --refGTF   "$NEXT_REF" \
    --refFasta "$FASTA" \
    -d "$OUT_DIR" \
    --report skip \
    -l DEBUG \
    --skipORF \
    --cpus 1


  # 1-3 前回と今回の分類の対応表を作る
  scripts/01-3.py \
    --old $PREVIOUS_RES \
    --new ${OUT_DIR}isoforms_classification.txt \
    --gtf $ALL_TX_GTF \
    --out ${OUT_DIR}new_old_classification_comparison.txt \
    --out_bed ${OUT_DIR}new_old_different_txs.bed
done


# 2：全ラウンドの分類の対応を合わせて一つにする
mkdir -p results/02
final_round="${rounds[@]: -1}"

scripts/02.py \
  --base_dir results/01/ \
  --n_rounds $final_round \
  --filename new_old_classification_comparison.txt \
  --output results/02/new_old_classification_comparison.merged.txt


# 3：エクソン領域がかぶるnew loci txをグループ化する（同じ遺伝子と見なす）
conda activate misc_20250301
mkdir -p results/03

scripts/03.py \
  --new_loci_gtf results/01/$final_round/isoforms.target.gtf \
  --all_tx_gtf $ALL_TX_GTF \
  --out_overlap_txs results/03/overlap_txs.tsv \
  --out_overlap_txs_bed results/03/overlap_txs.bed \
  --out_nonoverlap_txs results/03/nonoverlap_txs.tsv \
  --out_nonoverlap_txs_bed results/03/nonoverlap_txs.bed
