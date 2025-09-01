### Download raw reads from ENA by bioproject and species name

```
curl -sG 'https://www.ebi.ac.uk/ena/portal/api/search' \
  --data-urlencode "result=read_run" \
  --data-urlencode 'query=study_accession=PRJEB5065 AND scientific_name="Serratia marcescens"' \
  --data-urlencode "fields=run_accession,scientific_name,fastq_ftp" \
  --data-urlencode "format=tsv" \
  > runs.tsv

cut -f3 runs.tsv | tail -n +2 | tr ';' '\n' | while read url; do
  wget "ftp://$url"
done
```

### Download from file with accessions
```
# accessions.txt contains one run per line (e.g., ERR11895879)
parallel -j 4 '/home/britto/tools/enaBrowserTools-1.7.1/python3/enaDataGet -f fastq -d fastq {}' :::: /home/britto/data/Sam/Serratia/all_raw/accessions.txt```

### Get metadata
# 1) Convert Excel to CSV (keeps your header)
xlsx2csv -d tab Map1.xlsx Map1.tsv

# 2) Pull unique BioSamples
awk -F'\t' 'NR>1{print $2}' Map1.tsv | sort -u > biosamples.txt
# (Assumes header is: run_accession  sample_accession  experiment_accession  study_accession)

# show what you’ll loop over
nl -ba biosamples.txt
# normalize line endings & trim whitespace just in case
tr -d '\r' < biosamples.txt | sed 's/[[:space:]]\+$//' | awk 'NF' > biosamples.clean.txt
mv biosamples.clean.txt biosamples.txt

printf "biosample_accession\tscientific_name\tcollection_date\tcountry\thost\tisolation_source\tfirst_public\n" > biosample_meta.tsv
while read BS; do
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${BS}&result=sample&format=tsv&fields=sample_accession,scientific_name,collection_date,country,host,isolation_source,first_public" \
  | sed '1d' >> biosample_meta.tsv
done < biosamples.txt``
```

### run fastqc and kraken2
```
#### 0) activate your tools
conda activate kraken2

#### 1) set paths & threads
READS_DIR="/home/britto/data/Sam/Serratia/all_raw/trimmed"
DB="/home/britto/reference_database/refseq_kraken"
THREADS=8
OUT="work"

#### 2) create output folders
mkdir -p "$OUT/fastqc" "$OUT/kraken" "$OUT/multiqc"

#### 3) loop over all R1 files (supports .fq.gz and .fastq.gz)
shopt -s nullglob
for r1 in "$READS_DIR"/*_1_val_1.fq.gz "$READS_DIR"/*_1_val_1.fastq.gz; do
  [ -e "$r1" ] || continue  # skip if no matches

  # figure out extension and sample name
  if [[ "$r1" == *.fastq.gz ]]; then
    ext=".fastq.gz"
  else
    ext=".fq.gz"
  fi
  sample=$(basename "${r1%_1_val_1$ext}")
  r2="$READS_DIR/${sample}_2_val_2$ext"

  if [[ ! -f "$r2" ]]; then
    echo "!! Missing mate for $sample: $r2" >&2
    continue
  fi

  echo "==> $sample"
  # FastQC
  fastqc --quiet -t "$THREADS" -o "$OUT/fastqc" "$r1" "$r2"

  # Kraken2 (report only)
  kraken2 --db "$DB" --threads "$THREADS" \
          --paired "$r1" "$r2" \
          --report "$OUT/kraken/${sample}.report"
done

#### 4) MultiQC summary
multiqc "$OUT" -o "$OUT/multiqc"

echo "Done. FastQC: $OUT/fastqc ; Kraken reports: $OUT/kraken ; MultiQC: $OUT/multiqc"
```

### Trimming
```
conda activate trimgalore
for r1 in \
    *_R1.fastq.gz *_R1_001.fastq.gz \
    *_R1.fq.gz    *_1.fastq.gz \
    *_1.fq.gz
do
    # Skip if the glob expanded to nothing
    [[ -f "$r1" ]] || continue

    # Derive the R2 mate based on the detected pattern
    if   [[ "$r1" == *_R1.fastq.gz ]];      then r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    elif [[ "$r1" == *_R1_001.fastq.gz ]];  then r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    elif [[ "$r1" == *_R1.fq.gz ]];         then r2="${r1/_R1.fq.gz/_R2.fq.gz}"
    elif [[ "$r1" == *_1.fastq.gz ]];       then r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    elif [[ "$r1" == *_1.fq.gz ]];          then r2="${r1/_1.fq.gz/_2.fq.gz}"
    else
        echo "Skipping unrecognized file: $r1"
        continue
    fi

    if [[ -f "$r2" ]]; then
        # Sample label (strip the R1 marker and extension variations)
        sample="$(basename "$r1" | sed -E 's/(_R1(_001)?|_1)\.f(ast)?q\.gz$//')"
        echo "Trimming $sample"
        trim_galore --paired --fastqc -o trimmed "$r1" "$r2"
    else
        echo "Warning: missing R2 mate for: $r1"
    fi
done

conda deactivate
```
### Bracken
```
cd /home/britto/data/Sam/Serratia/all_raw/work1/work/kraken

for file in *.report; do
    sample=$(basename "$file" .report)
    bracken -d /home/britto/reference_database/refseq_kraken \
            -i "$file" \
            -o "${sample}.bracken" \
            -r 150 \
            -l S
done
```

#### combine bracken reports in R
```
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

# --- set your folder here ---
bracken_dir <- "/home/britto/data/Sam/Serratia/all_raw/work1/work/kraken"

# list all .bracken files in that folder (set recursive = TRUE if they’re in subfolders)
files <- list.files(path = bracken_dir,
                    pattern = "\\.bracken$",
                    full.names = TRUE,
                    recursive = FALSE)

stopifnot(length(files) > 0)  # fail early if none found

# read and reshape: keep taxon name + new_est_reads, name column by sample (file stem)
dfs <- map(files, ~ read_tsv(.x, show_col_types = FALSE) %>%
             transmute(
               name,
               !!str_remove(basename(.x), "\\.bracken$") := new_est_reads
             ))

# full join across all samples; NA -> 0
merged <- reduce(dfs, full_join, by = "name") %>%
  mutate(across(-name, ~ replace_na(., 0)))

# write output next to the inputs (or change the path if you prefer)
out_path <- file.path(bracken_dir, "bracken_merged.tsv")
write_tsv(merged, out_path)

message("Wrote: ", out_path)
```

### Move selected trimmed reads
```
SRC="/home/britto/data/Sam/Serratia/all_raw/all"
DEST="/home/britto/data/Sam/Serratia/all_raw/include"
LIST="accessions.txt"

mkdir -p "$DEST"

while read acc; do
  [ -z "$acc" ] && continue   # skip empty lines
  echo "==> Moving $acc"
  for mate in 1 2; do
    file="$SRC/${acc}_${mate}_val_${mate}.fq.gz"
    if [[ -f "$file" ]]; then
      mv "$file" "$DEST/"
    else
      echo "   Missing: $file"
    fi
  done
done < "$LIST"
```

## Snippy
```
conda activate snippy
snippy-multi /home/britto/data/Sam/Serratia/all_raw/samples.tsv --ref /home/britto/data/Sam/Serratia/all_raw/db11/db11.cleaned.fna --cpus 24 > runme.sh
less runme.sh   # check the script makes sense
sh ./runme.sh   # leave it running over lunch
snippy-clean_full_aln core.full.aln > clean.full.aln
conda deactivate
```

### qc with snippy bam_files
```
# ====== CONFIG ======
BASE="/home/britto/data/Sam/Serratia/all_raw/snippy"   # parent containing <STRAIN>/snps.bam
OUT="snippy_coverage_qc.tsv"                            # output TSV in current working dir

# QC thresholds for optional flagging
TH_MEAN_PASS=20       # PASS if mean_depth >= 20, REVIEW if 10–20, FAIL if <10
TH_MEAN_REVIEW=10

TH_BREADTH10_PASS=90  # PASS if breadth_ge10x_pct >= 90, REVIEW if 80–90, FAIL if <80
TH_BREADTH10_REVIEW=80

TH_ZERO_PASS=5        # PASS if zero_cov_pct <= 5, REVIEW if 5–10, FAIL if >10
TH_ZERO_REVIEW=10

TH_MAXGAP_PASS=10000  # PASS if max_zero_gap_bp < 10k, REVIEW if 10–50k, FAIL if >=50k
SPLIT_REVIEW_GAP=50000

TH_NGAPS_PASS=0       # PASS if n_gaps_ge10kb == 0, REVIEW if 1–2, FAIL if >=3
TH_NGAPS_REVIEW=2
# ====================

# Header (always overwrite)
echo -e "sample\tgenome_bp\tpositions_reported\ttotal_depth_sum\tmean_depth\tmin_depth\tmax_depth\tbreadth_ge1x_pct\tbreadth_ge4x_pct\tbreadth_ge10x_pct\tzero_cov_pct\tmax_zero_gap_bp\tn_gaps_ge10kb\tqc_flag\tqc_reasons" > "$OUT"

# Find exactly one snps.bam per strain dir
mapfile -d '' BAMS < <(find "$BASE" -mindepth 2 -maxdepth 2 -type f -name 'snps.bam' -print0 | sort -z)

if (( ${#BAMS[@]} == 0 )); then
  echo "No snps.bam files found under: $BASE" >&2
  exit 1
fi

# Process each BAM
for BAM in "${BAMS[@]}"; do
  SAMPLE="$(basename "$(dirname "$BAM")")"

  # Genome length from BAM header
  GENOME_BP=$(samtools view -H "$BAM" | awk 'BEGIN{sum=0} /^@SQ/ { if (match($0,/LN:([0-9]+)/,m)) sum+=m[1] } END{print sum+0}')
  if [ "${GENOME_BP:-0}" -le 0 ]; then
    echo "WARN: $SAMPLE -> could not read genome length from $BAM header; skipping." >&2
    continue
  fi

  # One pass over per-base depth
  # - stats: mean/min/max
  # - breadth >=1/4/10x
  # - zero coverage percent
  # - longest 0x run & number of 0x runs >=10 kb
  # positions_reported = NR (rows emitted by samtools depth -a)
  read -r POSITIONS SUM MEAN MIN MAX B1 B4 B10 ZERO_PCT MAXZERO NGAPS <<EOF
$(samtools depth -a "$BAM" | awk -v g="$GENOME_BP" -v t1=1 -v t4=4 -v t10=10 '
  BEGIN{
    sum=0; min=-1; max=-1; c1=0; c4=0; c10=0;
    prev_chrom=""; zero_run=0; max_zero=0; gaps_ge10k=0;
  }
  {
    chrom=$1; pos=$2; d=$3+0;

    if (NR==1) { prev_chrom=chrom }

    if (chrom!=prev_chrom) {
      if (zero_run>0) {
        if (zero_run>max_zero) max_zero=zero_run;
        if (zero_run>=10000) gaps_ge10k++;
      }
      zero_run=0;
      prev_chrom=chrom;
    }

    if (d==0) {
      zero_run++;
    } else if (zero_run>0) {
      if (zero_run>max_zero) max_zero=zero_run;
      if (zero_run>=10000) gaps_ge10k++;
      zero_run=0;
    }

    sum+=d;
    if (min<0 || d<min) min=d;
    if (max<0 || d>max) max=d;
    if (d>=t1) c1++;
    if (d>=t4) c4++;
    if (d>=t10) c10++;
  }
  END{
    if (zero_run>0) {
      if (zero_run>max_zero) max_zero=zero_run;
      if (zero_run>=10000) gaps_ge10k++;
    }
    mean = (NR? sum/NR:0.0);
    breadth1 = (g? (c1/g*100.0):0.0);
    breadth4 = (g? (c4/g*100.0):0.0);
    breadth10 = (g? (c10/g*100.0):0.0);
    zero_pct = 100.0 - breadth1;
    printf "%d %.0f %.6f %.6f %.6f %.2f %.2f %.2f %.2f %d %d\n",
           NR, sum, mean, min, max, breadth1, breadth4, breadth10, zero_pct, max_zero, gaps_ge10k;
  }')
EOF

  # QC flag & reasons
  FLAG="PASS"; REASONS=""

  # mean depth
  awk -v m="$MEAN" -v p="$TH_MEAN_PASS" -v r="$TH_MEAN_REVIEW" '
    BEGIN{
      if (m < r)      print "FAIL\tlow_mean_depth<" r;
      else if (m < p) print "REVIEW\tmean_depth_" r "-" p;
      else            print "PASS\t";
    }' | {
      read flag reason; [ -n "$reason" ] && REASONS+="$reason;"
      [ "$flag" = "FAIL" ] && FLAG="FAIL"
      [ "$flag" = "REVIEW" ] && [ "$FLAG" = "PASS" ] && FLAG="REVIEW"
    }

  # breadth at 10x
  awk -v b10="$B10" -v p="$TH_BREADTH10_PASS" -v r="$TH_BREADTH10_REVIEW" '
    BEGIN{
      if (b10 < r)      print "FAIL\tbreadth_ge10x<" r;
      else if (b10 < p) print "REVIEW\tbreadth_ge10x_" r "-" p;
      else              print "PASS\t";
    }' | {
      read flag reason; [ -n "$reason" ] && REASONS+="$reason;"
      [ "$flag" = "FAIL" ] && FLAG="FAIL"
      [ "$flag" = "REVIEW" ] && [ "$FLAG" = "PASS" ] && FLAG="REVIEW"
    }

  # zero coverage pct
  awk -v z="$ZERO_PCT" -v pass="$TH_ZERO_PASS" -v rev="$TH_ZERO_REVIEW" '
    BEGIN{
      if (z > rev)      print "FAIL\tzero_cov_pct>" rev;
      else if (z > pass) print "REVIEW\tzero_cov_pct_" pass "-" rev;
      else               print "PASS\t";
    }' | {
      read flag reason; [ -n "$reason" ] && REASONS+="$reason;"
      [ "$flag" = "FAIL" ] && FLAG="FAIL"
      [ "$flag" = "REVIEW" ] && [ "$FLAG" = "PASS" ] && FLAG="REVIEW"
    }

  # max zero gap
  if (( MAXZERO >= SPLIT_REVIEW_GAP )); then
    FLAG="FAIL"; REASONS+="max_zero_gap>="SPLIT_REVIEW_GAP";"
  elif (( MAXZERO >= TH_MAXGAP_PASS )); then
    [ "$FLAG" = "PASS" ] && FLAG="REVIEW"
    REASONS+="max_zero_gap_"TH_MAXGAP_PASS"-"SPLIT_REVIEW_GAP";"
  fi

  # number of big gaps
  if (( NGAPS > TH_NGAPS_REVIEW )); then
    FLAG="FAIL"; REASONS+="n_gaps_ge10kb>="(TH_NGAPS_REVIEW+1)";"
  elif (( NGAPS > TH_NGAPS_PASS )); then
    [ "$FLAG" = "PASS" ] && FLAG="REVIEW"
    REASONS+="n_gaps_ge10kb>=1;"
  fi

  printf "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%s\t%s\n" \
    "$SAMPLE" "$GENOME_BP" "$POSITIONS" "$SUM" "$MEAN" "$MIN" "$MAX" "$B1" "$B4" "$B10" "$ZERO_PCT" "$MAXZERO" "$NGAPS" "$FLAG" "$REASONS" \
    >> "$OUT"
done

echo "Wrote $(realpath "$OUT")"

```







