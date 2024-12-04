---
layout: post
title: "Rosetta Docking"
date: 2024-11-31
categories: [Bioinformatics, Docking]
tags: [Bioinformatics, docking]
---


---


# **Rosetta Docking Options Guide**

## **1. Introduction**
Rosetta 도킹은 단백질-단백질 및 단백질-리간드 상호작용을 예측하기 위한 강력한 도구입니다. 이 가이드에서는 도킹 프로토콜 실행 시 사용할 수 있는 주요 옵션과 그 사용법을 정리합니다.

---

## **2. General Options**

### **Input/Output**
- **`-s`**: 입력 구조 파일(PDB)을 지정합니다.
  ```bash
  -s input.pdb
  ```
- **`-out:path:all`**: 출력 파일 경로를 지정합니다.
  ```bash
  -out:path:all /path/to/output/
  ```
- **`-nstruct`**: 생성할 구조의 수를 지정합니다.
  ```bash
  -nstruct 10
  ```

### **Score Functions**
- **`-score:weights`**: 사용할 점수 함수(weight set)를 지정합니다.
  ```bash
  -score:weights ref15
  ```

---

## **3. Docking Options**

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

### **Ligand Docking**
리간드 도킹 관련 옵션입니다.

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

## **4. Running Docking Protocol**

### **Basic Command**
아래는 Rosetta 도킹 프로토콜의 기본 실행 명령입니다:

```bash
./bin/docking_protocol.default.linuxgccrelease \
  -s input.pdb \
  -docking:partners A_B \
  -nstruct 10 \
  -out:path:all /path/to/output/
```

### **Ligand Docking Example**
리간드 도킹을 실행하는 예제:

```bash
./bin/docking_protocol.default.linuxgccrelease \
  -s receptor_ligand.pdb \
  -docking:ligand true \
  -docking:ligand:protocol full \
  -docking:ligand:soft_rep true \
  -nstruct 20 \
  -out:path:all /path/to/output/
```

### **Symmetric Docking Example**
대칭성을 고려한 도킹 예제:

```bash
./bin/docking_protocol.default.linuxgccrelease \
  -s symmetric_complex.pdb \
  -docking:symmetry true \
  -docking:symmetry:minimize_backbone true \
  -nstruct 10 \
  -out:path:all /path/to/output/
```

---

## **5. Tips for Using Docking Options**
1. **최적화된 파라미터 사용:** 기본값으로 실행해 보고, 필요 시 특정 파라미터를 조정하세요.
2. **출력 디렉토리 설정:** 출력 경로를 명시적으로 지정하여 결과를 관리하기 쉽게 만드세요.
3. **로그 파일 확인:** 실행 오류 시 Rosetta 로그(`ROSETTA_CRASH.log`)를 확인하여 원인을 파악하세요.

---

## **6. References**
- [Rosetta Documentation](https://www.rosettacommons.org/docs/latest/home)
- `doc/full-options-list.md`: Rosetta 소스 디렉토리에 포함된 옵션 목록 파일


---

## 결과에 영향을 미칠만한 Ligand Docking 옵션

### `-docking:ligand:protocol`
- **영향**:
  - `abbreviated`는 간소화된 도킹 프로토콜로 빠르게 결과를 제공하지만, 정밀도가 낮을 수 있습니다.
  - `full`은 고해상도 도킹을 수행하며, 더 정확한 결과를 제공합니다.
- **추천 상황**:
  - 초기 스크리닝에서는 `abbreviated`, 최종 구조 최적화에서는 `full`을 사용하세요.

---

### `-docking:ligand:soft_rep`
- **영향**:
  - 부드러운 반발력을 활성화하면 리간드와 단백질 간의 겹침을 더 관대하게 처리합니다.
  - 초기 구조의 품질이 낮거나 샘플링 범위를 넓히고 싶을 때 유용합니다.
- **추천 상황**:
  - 구조 예측의 초기 단계나 결합 포켓의 정확한 정의가 불확실한 경우에 사용하세요.

---

### `-docking:ligand:grid`
- **영향**:
  - 그리드 기반 계산은 도킹 속도를 크게 향상시킬 수 있습니다.
  - 복잡한 리간드와 단백질 상호작용을 다룰 때 효율적입니다.
- **추천 상황**:
  - 대규모 도킹 시뮬레이션이나 계산 자원이 제한적인 경우에 활성화하세요.

---

## 결과에 영향을 미칠만한 Symmetric Docking 옵션

### `-docking:symmetry`
- **영향**:
  - 대칭성을 고려하지 않으면 구조가 물리적, 화학적 현실성을 잃을 수 있습니다.
  - 대칭성을 활성화하면 계산 효율성과 정확도가 모두 향상됩니다.
- **추천 상황**:
  - 대칭성을 가지는 복합체(예: 이량체, 사량체)를 다룰 때 반드시 사용하세요.

---

### `-docking:symmetry:minimize_backbone`
- **영향**:
  - 백본 구조의 유연성을 추가하면 상호작용 면적을 최적화할 가능성이 높아집니다.
  - 그러나 계산 시간이 증가할 수 있습니다.
- **추천 상황**:
  - 복합체의 유연성을 포함하여 도킹 결과를 최적화해야 할 때 사용하세요.

---

### `-docking:symmetry:use_symmetry_definition`
- **영향**:
  - 정확한 대칭 정의를 사용하면 계산이 더 정밀해집니다.
  - 잘못된 대칭 정의 파일을 사용하면 오류가 발생하거나 비현실적인 결과를 초래할 수 있습니다.
- **추천 상황**:
  - 대칭 복합체가 명확하게 정의된 경우, 반드시 대칭 정의 파일을 활용하세요.
