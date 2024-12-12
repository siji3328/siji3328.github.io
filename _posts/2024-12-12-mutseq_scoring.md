---
layout: post
title: "MutSequence Scoring"
date: 2024-11-07
categories: [Bioinformatics, Tool]
tags: [Bioinformatics, Tool]
---




# RNA 접근성, GC 함량, CAI 계산 및 점수화

이 코드는 DNA 서열에서 GC 함량, 코돈 적응성 지수(CAI), RNA 접근성을 계산하고 점수화하여 최종 점수를 산출하는 Python 코드입니다. Google Colab 환경에서 실행되며, ViennaRNA를 사용해 RNA 접근성을 계산합니다.

---

## **1. 라이브러리 설치 및 데이터 업로드**
필요한 라이브러리를 설치하고 데이터를 업로드합니다.

```python
!pip install viennarna
!pip install pandas openpyxl

import RNA
from google.colab import files
import pandas as pd
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex

# 엑셀 파일 업로드
uploaded = files.upload()
file_name = list(uploaded.keys())[0]
df = pd.read_excel(file_name)
2. 주요 계산 함수
GC 함량 계산
서열의 GC 함량을 계산합니다.

python
코드 복사
def calculate_gc_content(seq):
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2)
코돈 적응성 지수 (CAI) 계산
E. coli의 코돈 사용 빈도표를 사용해 CAI를 계산합니다.

python
코드 복사
def calculate_cai(seq):
    cai_calc = CodonAdaptationIndex()
    cai_calc.set_cai_index(ecoli_codon_usage)
    trimmed_seq = seq[:len(seq) // 3 * 3]  # 3의 배수로 맞춤
    return round(cai_calc.cai_for_gene(trimmed_seq), 3)
RNA 접근성 계산
서열의 최소 자유 에너지(MFE)를 기반으로 RNA 접근성을 계산합니다.

python
코드 복사
def calculate_accessibility(seq):
    seq = seq.replace('T', 'U')[:50]  # DNA → RNA 변환 후 50bp로 제한
    if len(seq) < 3:
        return 0.5  # 기본 접근성 값 반환
    fc = RNA.fold_compound(seq)
    mfe_structure, mfe = fc.mfe()  # 최소 자유 에너지 계산
    mfe = float(mfe)  # 문자열 → 실수 변환
    return round(1 - (abs(mfe) / 100), 3)  # 접근성 계산
3. 점수화 함수
값을 기준으로 점수를 부여하는 일반화된 함수입니다.

python
코드 복사
def score(value, thresholds, scores):
    for threshold, score in zip(thresholds, scores):
        if value >= threshold:
            return score
    return scores[-1]
4. 데이터 처리
서열별로 GC, CAI, RNA 접근성 값을 계산하고 점수를 부여합니다.

python
코드 복사
# 점수 기준 설정
thresholds_cai = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5]
thresholds_gc = [50, 45, 40, 35, 30, 25, 20, 15, 10, 5]
thresholds_accessibility = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
scores = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

cai_values, gc_values, accessibility_values = [], [], []
cai_scores, gc_scores, accessibility_scores = [], [], []

for seq in df['seq']:
    gc_content = calculate_gc_content(seq)
    cai = calculate_cai(seq)
    accessibility = calculate_accessibility(seq)

    # 실제 값 저장
    gc_values.append(gc_content)
    cai_values.append(cai)
    accessibility_values.append(accessibility)

    # 점수 계산
    cai_scores.append(score(cai, thresholds_cai, scores))
    gc_scores.append(score(gc_content, thresholds_gc, scores))
    accessibility_scores.append(score(accessibility, thresholds_accessibility, scores))
5. 최종 점수 계산
가중치를 적용해 최종 점수를 계산합니다.

python
코드 복사
# 점수 데이터프레임에 추가
df['GC_Value'] = gc_values
df['CAI_Value'] = cai_values
df['Accessibility_Value'] = accessibility_values

df['GC_Score'] = gc_scores
df['CAI_Score'] = cai_scores
df['Accessibility_Score'] = accessibility_scores

# 가중치 적용 후 최종 점수 계산
CAI_WEIGHT, ACCESSIBILITY_WEIGHT, GC_WEIGHT = 2, 2, 1
df['Final_Score'] = (
    df['CAI_Score'] * CAI_WEIGHT +
    df['Accessibility_Score'] * ACCESSIBILITY_WEIGHT +
    df['GC_Score'] * GC_WEIGHT
)

# 소수점 두 자리로 포맷팅
df['Final_Score'] = df['Final_Score'].apply(lambda x: round(x, 2))
6. 결과 저장
최종 점수를 포함한 데이터를 Excel 파일로 저장합니다.

python
코드 복사
output_file = "final_scores_with_values.xlsx"
df.to_excel(output_file, index=False)
files.download(output_file)
출력 결과
GC_Value: GC 함량.
CAI_Value: 코돈 적응성 지수.
Accessibility_Value: RNA 접근성.
GC_Score, CAI_Score, Accessibility_Score: 각 값에 대한 점수.
Final_Score: 최종 점수.
이 코드는 RNA 서열의 생물학적 특성을 정량화하고 점수화하여 Excel 파일로 저장할 수 있도록 설계되었습니다. 😊