---
layout: post
title: "Nanopore seqeuncing"
date: 2024-12-04
categories: [Sequence analysis]
tags: [sequencing, nanopore]
---



# 1. 작업 환경 설정
# Miniconda 환경 활성화
conda activate nanoplot_env

# 작업 디렉토리로 이동
cd /mnt/c/Users/MARS/soyoon/nanopore/240229

# 2. 데이터 품질 시각화
# NanoPlot으로 품질 시각화
for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/*.fastq; do
    NanoPlot --fastq "$file" --outdir /mnt/c/Users/MARS/soyoon/nanopore/240229/nanoploted_fastq --title "$(basename "$file")"
done

# 3. 품질이 낮거나 짧은 리드 제거
# Filtlong 설치
sudo apt-get update && sudo apt-get install filtlong

# Filtlong으로 필터링
for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/*.fastq; do
    filtlong --min_length 1000 --min_mean_q 10 "$file" > "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/$(basename "$file")"
done

# 4. 필터링 후 품질 재확인
# 필터링된 FASTQ 파일 시각화
NanoPlot --fastq /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/*.fastq --outdir /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/nanoploted_fastq

# 5. 레퍼런스 서열에 리드 정렬
# Minimap2와 SAMtools를 사용한 정렬
for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/*.fastq; do
    minimap2 -ax map-ont /mnt/c/Users/MARS/soyoon/nanopore/240229/refseq.fasta "$file" > "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/$(basename "$file" .fastq).sam"
    samtools view -Sb "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/$(basename "$file" .fastq).sam" | samtools sort -o "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/$(basename "$file" .fastq).sorted.bam"
    samtools index "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/$(basename "$file" .fastq).sorted.bam"
done

# 6. 정렬 품질 평가
# SAMtools로 매핑 품질 평가
samtools flagstat /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/*.sorted.bam > /mnt/c/Users/MARS/soyoon/nanopore/240229/mapping_stats.txt

# 7. 변이 탐지
# Longshot으로 변이 탐지
longshot -F /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/*.sorted.bam -o /mnt/c/Users/MARS/soyoon/nanopore/240229/variants.vcf

# 8. 커버리지 분석
# SAMtools로 커버리지 계산
samtools depth /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/*.sorted.bam > /mnt/c/Users/MARS/soyoon/nanopore/240229/coverage.txt

# 커버리지 평균 계산
awk '{sum += $3; count++} END {print "Average coverage:", sum / count}' /mnt/c/Users/MARS/soyoon/nanopore/240229/coverage.txt

# Python으로 커버리지 시각화
python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt

coverage_file = "/mnt/c/Users/MARS/soyoon/nanopore/240229/coverage.txt"
data = pd.read_csv(coverage_file, sep="\t", header=None, names=["Chromosome", "Position", "Depth"])
plt.figure(figsize=(15, 5))
plt.plot(data["Position"], data["Depth"], color="blue")
plt.xlabel("Genomic Position")
plt.ylabel("Read Depth")
plt.title("Coverage Plot")
plt.grid(True)
plt.show()
EOF

# 9. VCF에서 최종 서열 생성
# VCF 파일을 사용해 참조 서열 수정
bcftools consensus -f /mnt/c/Users/MARS/soyoon/nanopore/240229/refseq.fasta -o /mnt/c/Users/MARS/soyoon/nanopore/240229/fasta/consensus.fasta /mnt/c/Users/MARS/soyoon/nanopore/240229/variants.vcf

# 10. FASTQ에서 FASTA로 변환
mkdir -p /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/fastq_fasta
for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/*.fastq; do
    seqtk seq -a "$file" > "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/fastq_fasta/$(basename "$file" .fastq).fasta"
done

# 11. 태그 기반 리드 처리
python3 <<EOF
import os
from Bio import SeqIO

input_dir = "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/fastq_fasta"
output_dir = "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/tagged_reads"
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    if filename.endswith(".fasta"):
        filepath = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, f"tagged_{filename}")
        with open(filepath, "r") as fasta_file, open(output_path, "w") as tagged_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                record.id = f"Tagged_{record.id}"
                SeqIO.write(record, tagged_file, "fasta")
EOF

# 12. 최종 보고서 생성
echo "File Name, Total Reads, Mapped Reads, Percentage Mapped" > /mnt/c/Users/MARS/soyoon/nanopore/240229/mapping_summary.txt
for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/minimaped_bam/*.sorted.bam; do
    result=$(samtools flagstat "$file")
    total=$(echo "$result" | grep "in total" | awk '{print $1}')
    mapped=$(echo "$result" | grep "mapped (" | awk '{print $1}')
    percentage=$(echo "$result" | grep "mapped (" | awk -F'[()]' '{print $2}')
    echo "$(basename "$file"), $total, $mapped, $percentage" >> /mnt/c/Users/MARS/soyoon/nanopore/240229/mapping_summary.txt
done
