---
layout: post
title: "RFAA (RoseTTAFold-All-Atom)"
date: 2024-11-07
categories: [Bioinformatics, Model]
tags: [RFAA, RoseTTAFold, Bioinformatics, Protein Structure]
---

## RFAA (RoseTTAFold-All-Atom) 

### 1. Docker 프로그램 실행
   - 바탕화면에 있는 Docker 아이콘을 더블 클릭하여 Docker 프로그램을 실행합니다.

### 2. 터미널 열기
   - 바탕화면의 MARS 폴더로 이동한 후, RFAA 폴더를 찾습니다.
   - RFAA 폴더에서 마우스 오른쪽 버튼을 클릭하여 **'터미널에서 열기'**를 선택합니다.

### 3. Docker 실행 및 GPU 테스트
   - 아래 명령어를 차례로 입력하여 Docker와 GPU 인식을 확인합니다:

     ```bash
     docker run hello-world  # Docker 실행 테스트
     docker run --gpus all nvidia/cuda:11.8.0-runtime-ubuntu22.04 nvidia-smi  # GPU 인식 테스트
     ```

### 4. RFAA Docker 이미지 빌드 및 실행
   - 아래 명령어로 Docker 이미지를 빌드하고, GPU 지원을 활성화한 상태로 컨테이너를 실행합니다:

     ```bash
     docker build -t rfaa .  # Docker 이미지 빌드 (약 11분 소요)
     docker run --gpus all -it rfaa /bin/bash  # GPU 지원과 함께 컨테이너 실행
     ```

### 5. 환경 활성화 및 TensorFlow GPU 인식 테스트
   - 컨테이너 내에서 아래 명령어를 실행하여 RFAA 환경을 활성화하고, TensorFlow가 GPU를 인식하는지 확인합니다:

     ```bash
     source activate RFAA  # 환경 활성화
     python -c "import tensorflow as tf; print('Num GPUs Available:', len(tf.config.list_physical_devices('GPU')))"
     ```

   - **Num GPUs Available: 1**라고 표시되면 다음 단계로 진행합니다.  
   - **Num GPUs Available: 0**이라면 아래 명령어를 추가로 입력하여 TensorFlow를 업그레이드하고 환경 변수를 설정합니다:

     ```bash
     mamba install tensorflow-gpu=2.11  # TensorFlow 업그레이드
     docker run --gpus all -e NVIDIA_VISIBLE_DEVICES=all -e NVIDIA_DRIVER_CAPABILITIES=compute,utility -it rfaa /bin/bash  # 환경변수 설정
     ```

   - **Num GPUs Available: 1**이 표시되면 성공적으로 GPU가 인식된 것입니다!

### 6. RFAA 모델 예측 테스트
   - RFAA 모델 예측을 위한 파일을 다운로드하고, 예시로 단백질 예측을 수행해봅니다:

     ```bash
     git clone https://github.com/baker-laboratory/RoseTTAFold-All-Atom.git
     cd RoseTTAFold-All-Atom
     python -m rf2aa.run_inference --config-name protein  # 단백질 예측 테스트
     ```

### 7. 내 서열로 테스트하기
   - 위 단계를 완료한 후, 내 서열로 예측을 진행할 수 있습니다.
