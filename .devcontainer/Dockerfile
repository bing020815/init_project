# ==========
# Bioinfo Python Base Dockerfile
# Python 3.11 + JupyterLab + Notebook + 基本分析繪圖套件
# 無 R，適合未來做多系統擴充（機器學習、資料科學、API、可視化、Dash/Streamlit App）
# 可部署在 Linux / macOS / Windows 共通 Docker 環境
# ==========

FROM python:3.11-slim

# ==========
# OS 基本套件
# ==========
RUN apt-get update && apt-get install -y \
    gcc \
    python3-dev \
    curl \
    wget \
    git \
    vim \
    sudo \
    build-essential \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ==========
# Python 環境變數
# ==========
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# ==========
# 設定工作目錄
# ==========
WORKDIR /app

# ==========
# 安裝基本 Python 開發套件（支援 Notebook + 分析 + 繪圖）
# ==========
RUN pip install --no-cache-dir \
    ipykernel \
    jupyterlab \
    notebook \
    matplotlib \
    seaborn \
    pandas \
    numpy \
    scikit-learn

# ==========
# 複製 requirements.txt（如有）
# ==========
COPY ../requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt || true

# ==========
# 建立 Jupyter Notebook 設定（無密碼）
# ==========
RUN jupyter notebook --generate-config && \
    echo "c.NotebookApp.allow_root = True" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.token = ''" >> /root/.jupyter/jupyter_notebook_config.py

# ==========
# 預設啟動指令
# ==========
CMD ["bash"]