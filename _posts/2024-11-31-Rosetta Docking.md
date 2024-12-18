---
layout: post
title: "RosettaDock-5.0"
date: 2024-11-31
categories: [Bioinformatics, Docking]
tags: [Bioinformatics, docking]
---


---


# **RosettaDock-5.0**

- RosettaDock은 몬테카를로(Monte Carlo) 기반 다단계 도킹 알고리즘을 사용

---

## **1. RosettaDock-5.0 Quick Start**

### **입력 파일 준비**
- **단백질 파트너의 초기 위치를 추정**하여 서로의 표면이 마주 보도록 설정.
- **PyMOL**을 사용하여 분자를 번역 및 회전 가능:
  - **PyMOL 팁**:
    - 오른쪽 패널에서 "편집 모드" 활성화.
    - `Shift + 좌클릭 + 중간 버튼`으로 이동 및 회전.
    - `File → Export Molecule`로 PDB 저장.

### **기본 실행 명령**

- **`-s`**: 입력 구조 파일(PDB)을 지정합니다.
- **`-out:path:all`**: 출력 파일 경로를 지정합니다.
- **`-nstruct`**: 생성할 구조의 수를 지정합니다.

  ```bash
./bin/docking_protocol.default.linuxgccrelease \
  -s [파일 경로 및 입력파일] \
  -docking: partners [도킹 파트너] \
  -out:path:all [결과 파일 경로] \
  -nstruct [실행 횟수] \
  ```

### **예시**

  ```bash
./bin/docking_protocol.default.linuxgccrelease \
    -s input/chemi_A_B_rosedock.pdb \
    -docking:partners A_B \
    -out:path:all /mnt/c/Users/MARS/rosettadock/main/source/output/ \
    -nstruct 10
  ```
---

## **2. 결과 파일**

### **구조 파일 (`.pdb`)**
- **PDB 파일**: 각 도킹 후보 모델에 대해 생성.

### **점수 파일 (`.sc`)**
- **점수 파일(scorefile)**: 모든 생성된 모델의 에너지 및 메트릭 요약.

---

## **3. 결과 파일 해석**

- 성공적 도킹 기준: 상위 5개의 모델 중 **3개 이상이 "허용" 이상**일 경우 성공적인 도킹 (`N5 >= 3`).

| 지표       | 설명                                                                           |
|------------|--------------------------------------------------------------------------------|
| **Total**  | 복합체의 전체 에너지.                                                         |
| **I_sc**   | 인터페이스 점수 (결합 인터페이스의 에너지). 일반적으로 -5 ~ -10 사이의 값 추천. |

### **CAPRI 평가 지표**
도킹 모델의 정확성을 다음 기준으로 평가:
- **0**: 부정확 모델 (검정)
- **1**: 허용 수준 모델 (노랑)
- **2**: 중간 품질 모델 (빨강)
- **3**: 고품질 모델 (초록)

---

## **후처리**
1. **점수 파일 정렬**:
   - `Total Score` 기준으로 점수 파일을 정렬.
   - `I_sc` 값 확인: 좋은 데코이는 보통 -5 ~ -10 범위에 속함.

2. **클러스터링**:
   - 상위 **200개**의 데코이를 **RMSD** 기준으로 클러스터링.
   - 클러스터 크기와 점수를 비교하여 최적의 구조를 선택.

3. **데코이 수**:
   - **글로벌 도킹**: 최소 10,000개의 데코이 생성 필요. 이상적으로는 100,000개 생성 권장.
   - **로컬 도킹**: 최소 1,000개의 데코이 생성.

4. **추가 점수 계산**:
   - 인터페이스 점수가 표시되지 않는 경우 `-score:docking_interface_score 1` 플래그 사용.

5. **스크립트 사용**:
   - 클러스터링 및 후처리를 자동화하려면 Rosetta++ 튜토리얼에서 제공하는 스크립트 활용.

---

## **4. 옵션 설명**

### **General Docking**
- **`-docking:partners`**: 도킹할 체인을 지정합니다. 
  ```bash
  -docking:partners A_B
  ```
  - `A`는 리간드 체인, `B`는 수용체 체인입니다.

- **`-docking:randomize1`**, **`-docking:randomize2`**: 각각의 체인을 초기 위치에서 무작위화합니다.
  ```bash
  -docking:randomize1 true
  -docking:randomize2 true
  ```

- **`-docking:low_res_protocol_only`**: 저해상도 도킹 프로토콜만 실행합니다.
  ```bash
  -docking:low_res_protocol_only true
  ```
---


---

### **리간드 도킹 관련 옵션**

- **`-docking:ligand`**: 리간드 도킹 작업을 활성화합니다.
  ```bash
  -docking:ligand true
  ```

- **`-docking:ligand:protocol`**: 리간드 도킹 프로토콜을 지정합니다.
  ```bash
  -docking:ligand:protocol full
  ```

- **`-docking:ligand:soft_rep`**: 부드러운 반발력을 활성화합니다.
  ```bash
  -docking:ligand:soft_rep true
  ```

#### **Grid-Based Ligand Docking**
- **`-docking:ligand:grid`**: 그리드 기반 계산을 활성화합니다.
  ```bash
  -docking:ligand:grid true
  ```

- **`-docking:ligand:grid_kin`**: 생성된 그리드를 `.kin` 파일 형식으로 저장합니다.
  ```bash
  -docking:ligand:grid_kin grid_output.kin
  ```

- **`-docking:ligand:grid_map`**: 그리드를 BRIX 맵 형식으로 저장합니다.
  ```bash
  -docking:ligand:grid_map grid_output.map
  ```

---

### **예제**

```bash
./bin/docking_protocol.default.linuxgccrelease \
  -s receptor_ligand.pdb \
  -docking:ligand true \
  -docking:ligand:protocol full \
  -docking:ligand:soft_rep true \
  -nstruct 20 \
  -out:path:all /path/to/output/
```
---

### **Symmetric Docking**
대칭성을 고려한 도킹 시 사용되는 옵션입니다.

- **`-docking:symmetry`**: 대칭성을 활성화합니다.
  ```bash
  -docking:symmetry true
  ```

- **`-docking:symmetry:minimize_backbone`**: 단백질 백본의 최소화를 허용합니다.
  ```bash
  -docking:symmetry:minimize_backbone true
  ```

---

### **예제**

```bash
./bin/docking_protocol.default.linuxgccrelease \
  -s symmetric_complex.pdb \
  -docking:symmetry true \
  -docking:symmetry:minimize_backbone true \
  -nstruct 10 \
  -out:path:all /path/to/output/
```

---

## 6. Links

1. **Github:** <https://github.com/ccsb-scripps/AutoDock-GPU>
2. **RosettaCommons:** <https://rosettacommons.org/>
3. **RosettaCommons github:**<https://github.com/RosettaCommons>
4. **Documentation:** <https://docs.rosettacommons.org/docs/latest/application_documentation/docking/docking-protocol>

