mkdir -p 01.nfcore_smrnaseq/trimmed_fastq_10-20nt

find work_smrnaseq \( -name "*.trim.gz" \) \
  -printf "%T@ %p\n" \
  | sort -n \
  | awk '{
      path=$2
      n=split(path,a,"/")
      file=a[n]
      latest[file]=path
    }
    END {for (file in latest) print latest[file]}' \
  | while read -r f; do
      sample=$(basename "$f" .trim.gz)
      ln -sf "../../$f" "01.nfcore_smrnaseq/trimmed_fastq_10-20nt/${sample}.trim.gz"
      echo "$sample -> $f"
    done