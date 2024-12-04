---
layout: post
title: "Autodock-GPU"
date: 2024-11-07
categories: [Bioinformatics, Docking]
tags: [Bioinformatics, docking]
---



---
# AutoDock-GPU 사용법 및 옵션 정리

## 1. AutoDock-GPU 실행 방법

AutoDock-GPU를 사용하여 도킹 작업을 수행하려면 다음 단계를 따릅니다:

### **필요한 파일**
- **Receptor Grid Map 파일** (`.fld`): `--ffile`로 지정.
- **Ligand 파일** (`.pdbqt`): `--lfile`로 지정.
- 기타 옵션에 따라 필요에 따라 추가 파일을 준비합니다.

### **실행 명령어 형식**

```bash
./bin/autodock_gpu_64wi \
--ffile [Receptor Grid Map 경로] \
--lfile [Ligand 파일 경로] \
--nrun [실행 횟수] \
--resnam [결과 파일 경로] \
> [로그 파일 경로]
```

### **예시**

```bash
./bin/autodock_gpu_64wi \
--ffile /mnt/c/Users/MARS/AutoDock-GPU/input/soyoon/1ac8/1ac8_protein.maps.fld \
--lfile /mnt/c/Users/MARS/AutoDock-GPU/input/soyoon/1ac8/1ac8_ligand.pdbqt \
--nrun 10 \
--resnam /mnt/c/Users/MARS/AutoDock-GPU/output/soyoon/1ac8_result.dlg \
> /mnt/c/Users/MARS/AutoDock-GPU/output/soyoon/1ac8_log.txt
```

---

## 2. 주요 옵션 설명

### **입력 파일 관련**
  ```
| 옵션               | 설명                                            | 기본값        |
|--------------------|------------------------------------------------|--------------|
| `--lfile` / `-L`  | Ligand 파일 경로 (`.pdbqt` 파일)                 | 필수         |
| `--ffile` / `-M`  | Receptor Grid Map 파일 경로 (`.fld` 파일)        | 필수         |
| `--flexres` / `-F`| Flexible residue 파일 경로 (`.pdbqt`)            | 없음         |
| `--filelist` / `-B`| Batch 모드로 처리할 파일 리스트 지정             | 없음         |
  ```
### **출력 관련**
  ```
| 옵션                | 설명                                            | 기본값        |
|---------------------|------------------------------------------------|--------------|
| `--resnam` / `-N`   | 출력 로그 파일 이름                             | Ligand 이름  |
| `--contact_analysis`| 도킹 후 거리를 기반으로 한 접촉 분석 수행         | `0` (비활성) |
| `--dlgoutput`       | `.dlg` 형식 결과 파일 생성 여부                  | `1` (활성)   |
| `--xmloutput`       | `.xml` 형식 결과 파일 생성 여부                  | `1` (활성)   |
| `--output-cluster-poses` | 클러스터링 결과에서 출력할 자세 수            | `0` (모두)   |
  ```
### **탐색 알고리즘 설정**
  ```
| 옵션                | 설명                                            | 기본값        |
|---------------------|------------------------------------------------|--------------|
| `--nrun` / `-n`     | LGA (Local Search Genetic Algorithm) 실행 횟수  | `20`         |
| `--nev` / `-e`      | LGA 실행당 최대 평가 횟수                       | `2500000`    |
| `--ngen` / `-g`     | 세대 수 (Generations)                           | `42000`      |
| `--psize` / `-p`    | LGA 인구 크기                                   | `150`        |
| `--lsrat`           | Local Search 비율 (%)                           | `100`        |
| `--crat`            | Crossover 비율 (%)                              | `80`         |
| `--mrat`            | Mutation 비율 (%)                               | `2`          |
| `--dmov`            | 최대 LGA 이동 변위 (`Å`)                        | `6.0`        |
| `--dang`            | 최대 LGA 회전 변위 (`°`)                        | `90.0`       |
  ```
### **스코어링 및 에너지 설정**
  ```
| 옵션                | 설명                                            | 기본값        |
|---------------------|------------------------------------------------|--------------|
| `--smooth`          | van der Waals 상호작용의 스무딩 파라미터 (`Å`)   | `0.5`        |
| `--elecmindist`     | 최소 전기적 상호작용 거리 (`Å`)                 | `0.01`       |
  ```
### **기타 설정**
  ```
| 옵션                | 설명                                            | 기본값        |
|---------------------|------------------------------------------------|--------------|
| `--devnum` / `-D`   | CUDA/OpenCL 디바이스 번호                        | `1`          |
| `--seed` / `-s`     | 랜덤 시드 (세 개의 정수 지정 가능)               | `시간, PID`  |
| `--autostop`        | 수렴 기준에 따라 자동 정지 여부                   | `1` (활성)   |
| `--asfreq`          | AutoStop 테스트 빈도 (세대 수)                   | `5`          |
  ```
---

## 3. 결과 파일 설명

### **결과 파일 (`.dlg`)**
- **경로:** `--resnam` 옵션에서 지정한 위치.
- **주요 정보:**
  - **Estimated Free Energy of Binding:** 도킹된 리간드의 예상 결합 자유 에너지.
  - **RMSD (Root Mean Square Deviation):** 도킹 자세의 변동성을 보여주는 값.
  - **Clusters:** 각 클러스터의 평균 에너지 및 발생 빈도.

### **로그 파일 (`.txt`)**
- **경로:** 명령 끝에 `>`로 지정한 파일.
- **내용:**
  - 명령 실행 과정.
  - CUDA 디바이스 및 메모리 사용량.
  - 오류 및 디버그 메시지.

---

## 4. 결과 파일 해석 예시

### **결과 파일 주요 섹션**

#### **Free Energy of Binding**
- 파일의 `Estimated Free Energy of Binding` 항목에서 리간드의 결합 에너지를 확인합니다.
- 값이 **더 낮을수록** 결합이 안정적입니다.
  ```
  DOCKED: USER    Estimated Free Energy of Binding    =  -8.25 kcal/mol
  ```

#### **RMSD Table**
- RMSD 값을 통해 도킹 자세의 변동성을 평가합니다.
- RMSD가 **2.0 Å 이하**인 클러스터를 주로 사용합니다.

#### **클러스터 분석**
- Clustering Histogram 섹션을 확인하여 가장 빈도가 높은 클러스터를 찾습니다.
  ```
  Clus | Lowest    | Run | Mean      | Num | Histogram
  Rank | Binding   |     | Binding   | Clus|    5    10   15   20
  _____|___________|_____|___________|_____|____:____|____:____
     1 |     -8.31 |  10 |     -8.05 |  10 |##########
  ```

---

## 5. 옵션 활용 예시

### **`--nrun` 조정**
- 실행 횟수를 늘리면 더 많은 자세를 탐색하여 결합 안정성을 확인할 수 있습니다.
  ```bash
  ./bin/autodock_gpu_64wi \
  --ffile [Receptor Grid Map 경로] \
  --lfile [Ligand 파일 경로] \
  --nrun 50 \
  --resnam [결과 파일 경로]
  ```

### **`--smooth` 사용**
- van der Waals 스무딩 파라미터를 변경하여 더 부드러운 상호작용을 모델링할 수 있습니다.
  ```bash
  ./bin/autodock_gpu_64wi \
  --ffile [Receptor Grid Map 경로] \
  --lfile [Ligand 파일 경로] \
  --smooth 1.0 \
  --resnam [결과 파일 경로]
  ```

### **Flexible Residues 추가**
- 리간드와 결합부위를 유연하게 처리할 때 사용합니다.
  ```bash
  ./bin/autodock_gpu_64wi \
  --ffile [Receptor Grid Map 경로] \
  --lfile [Ligand 파일 경로] \
  --flexres [Flexible Residue 파일 경로] \
  --resnam [결과 파일 경로]
  ```

---

## 6. 주의사항
1. **경로 확인:** 입력 파일과 출력 경로가 올바른지 반드시 확인하세요.
2. **출력 디렉토리:** 결과가 저장될 폴더가 존재하지 않으면 명령어가 실패합니다.
3. **GPU 사용량:** GPU 메모리 부족으로 인해 도킹이 중단되지 않도록 주의하세요.