---
layout: post
title: "Oligo designer"
date: 2024-11-08
categories: [Bioinformatics, Tool]
tags: [Oliognucleotide, Oiligo designer]
---

# Oligo designer (PCA를 위한) 코드 사용법 


- 중첩되는 서열의 Tm 값을 만족하는 조각을 만들어 PCA(Polymerase Chain Assembly) 프라이머를 설계하는 Python 코드
- Oligo Design Tool을 직접 실행해 보세요. 클릭 ->
[![Colab에서 실행하기](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/siji3328/OligoDesign/blob/master/oligodesigner.ipynb)


### 1. Tm(녹는 온도) 계산 (Wallace 공식)

`calculate_tm` 함수는 DNA 조각의 Tm을 Wallace 공식을 사용해 계산합니다.

\[
Tm = 2 \times (A + T) + 4 \times (G + C)
\]

각 염기서열의 A, T, G, C 개수를 센 후 공식을 적용하여 Tm을 계산합니다.

### 2. 프라이머 설계 함수: `design_pca_primers`

`design_pca_primers` 함수는 PCA를 위한 DNA 조각들을 생성하고, 각 조각의 Tm 값을 확인하여 설정된 Tm 범위에 맞는지 검사합니다.

#### 매개변수

- **`sequence`**: 전체 DNA 서열
- **`min_fragment_length`**, **`max_fragment_length`**: 조각의 최소 및 최대 길이
- **`min_tm`**, **`max_tm`**: 중첩 서열의 Tm 범위

#### 작동 방식

1. `while` 루프를 통해 서열을 조각 단위로 나눕니다.
2. 각 조각의 Tm 값이 `min_tm`과 `max_tm` 범위 내에 있으면 조각을 선택합니다.
3. 선택된 각 조각에 대해 Forward와 Reverse 프라이머를 설계하고, 결과는 `primers` 리스트에 저장됩니다.

### 3. 사용 예시

1. 전체 DNA 서열(`sequence`)을 입력하고, 프라이머를 설계하기 위해 `design_pca_primers` 함수를 호출합니다.
2. 각 조각에 대해 Fragment 번호, 서열, 길이, Tm, Forward 및 Reverse 프라이머를 출력합니다.




## Code

```python
def calculate_tm(sequence):
    # Wallace 공식: Tm = 2(A+T) + 4(G+C)
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    tm = 2 * (a_count + t_count) + 4 * (g_count + c_count)
    return tm

def design_pca_primers(sequence, min_fragment_length=75, max_fragment_length=95, min_tm=50, max_tm=55):
    primers = []
    fragments = []

    i = 0
    while i < len(sequence):
        selected_fragment = None
        overlap_tm = 0
        for fragment_length in range(max_fragment_length, min_fragment_length - 1, -1):
            if i + fragment_length > len(sequence):
                fragment_length = len(sequence) - i

            fragment = sequence[i:i + fragment_length]

            if i + fragment_length < len(sequence):
                overlap_sequence = sequence[i + fragment_length - (fragment_length - min_fragment_length):i + fragment_length]
                overlap_tm = calculate_tm(overlap_sequence)

                if min_tm <= overlap_tm <= max_tm:
                    selected_fragment = fragment
                    i += fragment_length - (fragment_length - min_fragment_length)
                    break
            else:
                selected_fragment = fragment
                i += fragment_length

        if selected_fragment:
            fragments.append((selected_fragment, overlap_tm))

    # 짝수 개의 조각을 유지하기 위해 조정
    if len(fragments) % 2 != 0:
        fragments.pop()

    # 프라이머 설계 (Forward 및 Reverse)
    for i, (fragment, overlap_tm) in enumerate(fragments):
        forward_primer = fragment
        if i < len(fragments) - 1:
            reverse_primer = fragments[i + 1][0][:min_fragment_length]  # 다음 조각의 오버랩 서열을 사용
        else:
            reverse_primer = ""  # 마지막 조각의 Reverse Primer는 필요하지 않음

        primers.append({
            "fragment_number": i + 1,
            "fragment_sequence": fragment,
            "fragment_length": len(fragment),
            "overlap_tm": overlap_tm,
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer
        })

    return primers

# 제공된 전체 DNA 서열
sequence = "ATGAGCAAAGGTGAAGAACTGTTTACCGGCGTTGTGCCGATTCTGGTGGAACTGGATGGCGATGTGAACGGTCACAAATTCAGCGTGCGTGGTGAAGGTGAAGGCGATGCCACGATTGGCAAACTGACGCTGAAATTTATCTGCACCACCGGCAAACTGCCGGTGCCGTGGCCGACGCTGGTGACCACCCTGACCTATGGCGTTCAGTGTTTTAGTCGCTATCCGGATCACATGAAACGTCACGATTTCTTTAAATCTGCAATGCCGGAAGGCTATGTGCAGGAACGTACGATTAGCTTTAAAGATGATGGCAAATATAAAACGCGCGCCGTTGTGAAATTTGAAGGCGATACCCTGGTGAACCGCATTGAACTGAAAGGCACGGATTTTAAAGAAGATGGCAATATCCTGGGCCATAAACTGGAATACAACTTTAATAGCCATAATGTTTATATTACGGCGGATAAACAGAAAAATGGCATCAAAGCGAATTTTACCGTTCGCCATAACGTTGAAGATGGCAGTGTGCAGCTGGCAGATCATTATCAGCAGAATACCCCGATTGGTGATGGTCCGGTGCTGCTGCCGGATAATCATTATCTGAGCACGCAGACCGTTCTGTCTAAAGATCCGAACGAAAAACGGGACCACATGGTTCTGCACGAATATGTGAATGCGGCAGGTATTACGTGGAGCCATCCGCAGTTCGAAAAATAA"

# 프라이머 설계
primers = design_pca_primers(sequence, min_fragment_length=75, max_fragment_length=95, min_tm=50, max_tm=54)

# 결과 출력
for primer in primers:
    print(f"Fragment {primer['fragment_number']}:")
    print(f"  Sequence: {primer['fragment_sequence']}")
    print(f"  Fragment Length: {primer['fragment_length']} bp")
    print(f"  Overlap Tm: {primer['overlap_tm']} °C")
    print(f"  Forward Primer: {primer['forward_primer']}")
    print(f"  Reverse Primer: {primer['reverse_primer']}")
    print("\n")
```
