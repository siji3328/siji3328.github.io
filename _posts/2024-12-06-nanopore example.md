---
layout: post
title: "Nanopore seqeuncing example"
date: 2024-12-06
categories: [Bioinformatics, Sequence analysis]
tags: [sequencing, nanopore]
---

# NGS 분석 워크플로우

## 1. 작업환경 확인

### WSL에서 Miniconda 설치
WSL에서 터미널을 열고 `nanoplot_env` 환경을 활성화합니다.

conda activate nanoplot_env

작업할 디렉토리로 이동합니다.

cd /mnt/c/Users/MARS/soyoon/nanopore/240229

## 2. 데이터 품질 시각화

### NanoPlot: FASTQ 파일의 품질을 시각화
NanoPlot을 사용하여 FASTQ 파일의 품질, 리드 길이, 품질 점수(Q Score) 분포 등을 시각화합니다. 이 단계에서 리드의 전반적인 품질을 확인할 수 있습니다.

for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/*.fastq; do
    NanoPlot --fastq "$file" --outdir /mnt/c/Users/MARS/soyoon/nanopore/240229/nanoploted_fastq --title "$(basename "$file")"
done

**결과**:
- 품질 점수(Q Score) 히스토그램.
- 리드 길이 분포.
- 품질과 길이 간의 상관관계.

## 3. 데이터 처리: Filtlong으로 품질이 낮거나 짧은 리드 제거

### WSL에서 Filtlong 설치
sudo apt-get update
sudo apt-get install filtlong

### 필터링 실행
for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/*.fastq; do
    filtlong --min_length 1000 --min_mean_q 10 "$file" > "/mnt/c/Users/MARS/soyoon/nanopore/240229/qced_fastq/$(basename "$file")"
done

**결과**:
- 품질 점수(Q Score)와 리드 길이 기준에 따라 필터링된 FASTQ 파일.

## 4. 다시 데이터 품질 시각화하여 품질 향상 여부 확인

NanoPlot --fastq /mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/*.fastq --outdir /mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/nanoploted_fastq


## 5. 필터링된 FASTQ 파일 정렬

### Minimap2: 필터링된 모든 FASTQ 파일을 레퍼런스 서열에 정렬합니다.

for file in /mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/*.fastq; do
    # SAM 파일 생성
    minimap2 -ax map-ont /mnt/c/Users/MARS/soyoon/nanopore/240229/refseq.fasta "$file" > "/mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/$(basename "$file" .fastq).sam"
    
    # BAM 파일 생성 및 정렬
    samtools view -Sb "/mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/$(basename "$file" .fastq).sam" | samtools sort -o "/mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/$(basename "$file" .fastq).sorted.bam"
    
    # BAM 파일 인덱싱
    samtools index "/mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/$(basename "$file" .fastq).sorted.bam"
done

**정렬 과정**:
- SAM 파일은 임시로 생성되고, 이후 BAM 파일로 변환됩니다.
- 정렬된 BAM 파일은 `minimaped_bam` 디렉토리에 저장됩니다.

**BAM 파일 구조**:
- `*.sam`: 정렬된 SAM 파일.
- `*.sorted.bam`: 정렬된 BAM 파일.
- `*.sorted.bam.bai`: BAM 파일의 인덱스.

### 출력 파일 위치
- **입력**: `C:\Users\MARS\soyoon\nanopore\240229\filtlonged_fastq`의 FASTQ 파일.
- **출력**: 정렬된 BAM 파일과 인덱스는 `C:\Users\MARS\soyoon\nanopore\240229\filtlonged_fastq\minimaped_bam`에 저장됩니다.

`C:\Users\MARS\soyoon\nanopore\240229\filtlonged_fastq\minimaped_bam` 폴더에 저장된 BAM 파일들은 레퍼런스 서열에 정렬된 결과입니다. 이 데이터를 기반으로 변이 분석(Variant Calling), 커버리지 확인(Coverage Analysis), 또는 시각화(Visualization) 등의 후속 작업을 진행할 수 있습니다.

## 6. 후속 분석

### 정렬 품질 확인
BAM 파일의 매핑 품질(Mapping Quality)을 확인합니다.

samtools flagstat /mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/*.sorted.bam

**결과 예시**:
- 총 리드 수
- 매핑된 리드 수
- 중복된 리드 수

### 매핑률 확인
매핑률(%)을 통해 데이터가 얼마나 잘 매핑되었는지 평가할 수 있습니다. 이 값이 낮다면 레퍼런스 서열과 데이터가 적합하지 않을 수 있습니다.


## 6. 변이 분석 (Variant Calling)

### BAM 파일을 기반으로 레퍼런스 서열과 비교하여 변이를 탐지합니다.

### Longshot (SNP 및 Indel 탐지)
Longshot은 나노포어 데이터를 처리할 때 유용한 변이 탐지 도구입니다. 다음 명령어로 VCF 파일을 생성합니다:

longshot -F /mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/*.sorted.bam -o /mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/variants.vcf

### VCF 파일 확인
생성된 VCF 파일에서 변이 정보를 확인합니다:

less /mnt/c/Users/MARS/soyoon/nanopore/240229/filtlonged_fastq/minimaped_bam/variants.vcf





















