#!/bin/bash
#SBATCH --job-name=sqanti_array
#SBATCH --output=/path/to/sqanti/practice/logs/03/%A_%a.out
#SBATCH --error=/path/to/sqanti/practice/logs/03/%A_%a.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1


# 使い方：
# sbatch --array=1-7356:1%200 scripts/run_sqanti_array.sh


# 環境設定
source =/miniconda3/etc/profile.d/conda.sh
conda activate sqanti3


# 入力情報
PAIR_FILE="/path/to/sqanti/practice/results/02/pairs.txt"
GTF_DIR="/path/to/sqanti/practice/results/02/gtfs/"
REF_FASTA="/path/to/sqanti/SQANTI3/data/reference/GRCh38.p13_chr22.fasta"
SQANTI_SCRIPT="/path/to/sqanti/SQANTI3/sqanti3_qc.py"


# タスクIDに対応する行を取得
read -r GROUP_ID TX1 TX2 < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $PAIR_FILE)

echo "tx1: $TX1"
echo "tx2: $TX2"
echo "group: $GROUP_ID"

OUT_DIR="/path/to/sqanti/practice/results/03/${TX1}-${TX2}"
mkdir -p "$OUT_DIR"

: "${GTF_DIR:?Error: GTF_DIR is not set or empty}"
python "$SQANTI_SCRIPT" \
  --isoforms "$GTF_DIR/${TX1}.gtf" \
  --refGTF   "$GTF_DIR/${TX2}.gtf" \
  --refFasta "$REF_FASTA" \
  -o tx1_tx2 \
  -d "$OUT_DIR" \
  --report skip \
  -l DEBUG \
  --skipORF \
  --cpus 1
