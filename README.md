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
```
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




