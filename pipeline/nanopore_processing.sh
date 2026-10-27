#!/usr/bin/env bash
set -euo pipefail

: "${DATA_DIR:?Set DATA_DIR to the directory containing the sequencing data}"

THREADS="${THREADS:-8}"
DORADO_BIN="${DORADO_BIN:-dorado}"
DORADO_DEVICE="${DORADO_DEVICE:-cuda:0,1}"
STAGE="${1:-}"
GENOME="${DATA_DIR}/genome/hvo.fasta"

map_reads() {
  local reference="$1" reads="$2" output="$3"
  mkdir -p "$(dirname "$output")"

  minimap2 -y -ax map-ont --MD -t "$THREADS" "$reference" "$reads" \
    | samtools view -b -F 2048 - \
    | samtools sort -@ "$THREADS" -o "$output" -

  samtools index "$output"
  samtools flagstat "$output" > "${output%.bam}.flagstat.txt"
}

case "$STAGE" in
  basecall)
    : "${POD5_DIR:?Set POD5_DIR}"
    mkdir -p "$DATA_DIR/rebasecalling"

    "$DORADO_BIN" basecaller sup@v5.0.0 "$POD5_DIR" \
      --device "$DORADO_DEVICE" --kit-name SQK-PCB114-24 --no-trim \
      > "$DATA_DIR/rebasecalling/rebasecalling_pod5.bam"

    "$DORADO_BIN" summary "$DATA_DIR/rebasecalling/rebasecalling_pod5.bam" \
      > "$DATA_DIR/rebasecalling/summary.tsv"
    ;;

  demux)
    mkdir -p "$DATA_DIR/demux" "$DATA_DIR/fastq_pass"

    "$DORADO_BIN" demux --kit-name SQK-PCB114-24 --no-trim --emit-fastq \
      --output-dir "$DATA_DIR/demux" \
      "$DATA_DIR/rebasecalling/rebasecalling_pod5.bam"

    for file in "$DATA_DIR"/demux/*.fastq; do
      [[ -f "$file" ]] || continue
      barcode=$(basename "$file" | grep -oE 'barcode[0-9]+' | head -n 1) || continue
      gzip -c "$file" > "$DATA_DIR/fastq_pass/${barcode}.fastq.gz"
    done
    ;;

  trim)
    : "${PRIMER_FASTA:?Set PRIMER_FASTA to the custom pychopper primer sequences}"
    mkdir -p "$DATA_DIR/pychopper" "$DATA_DIR/fastq_pass_fl" "$DATA_DIR/cutadapt_VNP"

    for file in "$DATA_DIR"/fastq_pass/barcode*.fastq.gz; do
      [[ -f "$file" ]] || continue
      barcode=$(basename "$file" .fastq.gz)
      out="$DATA_DIR/pychopper/$barcode"
      mkdir -p "$out"

      pychopper -m edlib -k PCB114 -b "$PRIMER_FASTA" -U -y -t "$THREADS" \
        -r "$out/${barcode}_report.pdf" \
        -u "$out/${barcode}_unclassified.fastq" \
        -w "$out/${barcode}_rescued.fastq" \
        -S "$out/${barcode}_stats.txt" \
        "$file" "$DATA_DIR/fastq_pass_fl/${barcode}_full_length.fastq"

      cutadapt -a CTGTAGGCACCATCAAT -g ^GGG -j "$THREADS" \
        -o "$DATA_DIR/cutadapt_VNP/${barcode}.fl.cutVNP.fastq.gz" \
        "$DATA_DIR/fastq_pass_fl/${barcode}_full_length.fastq"
    done
    ;;

  map)
    for file in "$DATA_DIR"/cutadapt_VNP/barcode*.fl.cutVNP.fastq.gz; do
      [[ -f "$file" ]] || continue
      barcode=$(basename "$file" | cut -d. -f1)
      map_reads "$GENOME" "$file" \
        "$DATA_DIR/mapped_fl_only/$barcode/$barcode.sorted.bam"
    done
    ;;

  dedup)
    for bam in "$DATA_DIR"/mapped_fl_only/barcode*/barcode*.sorted.bam; do
      [[ -f "$bam" ]] || continue
      barcode=$(basename "$bam" | cut -d. -f1)
      out="$DATA_DIR/mapped_fl_dedup/$barcode"
      mkdir -p "$out"

      samtools view -h "$bam" \
        | awk 'BEGIN {OFS="\t"} /^@/ {print; next} {
            for (i = 12; i <= NF; i++) {
              if ($i ~ /^RX:Z:/) {
                umi = substr($i, 6)
                umi = substr(umi, 1, 28)
                while (length(umi) < 28) umi = umi "N"
                $i = "RX:Z:" umi
              }
            }
            print
          }' \
        | samtools view -b - \
        | samtools sort -@ "$THREADS" -o "$out/${barcode}_umi.sorted.bam" -

      samtools index "$out/${barcode}_umi.sorted.bam"
      umi_tools dedup --extract-umi-method=tag --umi-tag=RX \
        --stdin="$out/${barcode}_umi.sorted.bam" \
        --stdout="$out/${barcode}_umi_dedup.sorted.bam"
      samtools index "$out/${barcode}_umi_dedup.sorted.bam"

      samtools fastq -@ "$THREADS" "$out/${barcode}_umi_dedup.sorted.bam" \
        | cutadapt -e 0.1 -j "$THREADS" \
            -g '^TTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG' \
            -o "$out/${barcode}_deduplicated_reads_umi_removed.fastq.gz" -
    done
    ;;

  remap)
    for file in "$DATA_DIR"/mapped_fl_dedup/barcode*/*_deduplicated_reads_umi_removed.fastq.gz; do
      [[ -f "$file" ]] || continue
      barcode=$(basename "$file" | cut -d_ -f1)
      map_reads "$GENOME" "$file" \
        "$DATA_DIR/mapped_fl_only_umi_removed/$barcode/$barcode.sorted.bam"
    done
    ;;

  circular)
    : "${SAMPLE_INFO:?Set SAMPLE_INFO to sample_info.txt}"
    for file in "$DATA_DIR"/mapped_fl_dedup/barcode*/*_deduplicated_reads_umi_removed.fastq.gz; do
      [[ -f "$file" ]] || continue
      barcode=$(basename "$file" | cut -d_ -f1)

      for entry in '4 16S' '5 23S'; do
        read -r column rrna <<< "$entry"
        reference=$(awk -F, -v bc="$barcode" -v col="$column" '$1 == bc {print $col; exit}' "$SAMPLE_INFO")
        [[ -n "$reference" && "$reference" != NA ]] || continue

        map_reads "$DATA_DIR/genome/${reference}.fasta" "$file" \
          "$DATA_DIR/mapped_fl_only_umi_removed_circ${rrna}/$barcode/$barcode.sorted.bam"
      done
    done
    ;;

  coverage)
    for bam in "$DATA_DIR"/mapped_fl_only_umi_removed/barcode*/barcode*.sorted.bam; do
      [[ -f "$bam" ]] || continue
      barcode=$(basename "$bam" | cut -d. -f1)
      out="$DATA_DIR/bamCoverage/$barcode"
      mkdir -p "$out"

      for entry in 'pos reverse' 'neg forward'; do
        read -r label strand <<< "$entry"
        for end in end3 end5; do
          offset=-1
          [[ "$end" == end5 ]] && offset=1

          bamCoverage -b "$bam" \
            -o "$out/${barcode}.${end}.sorted.${label}.cpm.bedgraph" \
            --numberOfProcessors "$THREADS" -of bedgraph \
            --filterRNAstrand "$strand" --normalizeUsing CPM \
            --Offset "$offset" --binSize 1
        done
      done
    done
    ;;

  peaks)
    for condition in wt RNaseR; do
      sample1=barcode20
      sample2=barcode21
      if [[ "$condition" == RNaseR ]]; then
        sample1=barcode22
        sample2=barcode23
      fi

      for strand in pos neg; do
        for end in end3 end5; do
          first="$DATA_DIR/bamCoverage/$sample1/${sample1}.${end}.sorted.${strand}.cpm.bedgraph"
          second="$DATA_DIR/bamCoverage/$sample2/${sample2}.${end}.sorted.${strand}.cpm.bedgraph"
          out="$DATA_DIR/$end"
          prefix="$out/t${condition}.${strand}.peaks"
          mkdir -p "$out"

          termseq_peaks "$first" "$second" \
            --peaks "$prefix" --idr-threshold 0.05 --strand "$strand"

          bedtools intersect -wao -a "${prefix}.oracle.narrowPeak" \
            -b "$first" "$second" > "${prefix}.oracle.narrowPeak.counts"
        done
      done
    done
    ;;

  *)
    echo "Usage: DATA_DIR=/path/to/data bash pipeline/nanopore_processing.sh {basecall|demux|trim|map|dedup|remap|circular|coverage|peaks}" >&2
    exit 2
    ;;
esac
