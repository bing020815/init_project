# 使用官方的 Python 映像 3.9.2
FROM python:3.9-slim 

# 設置工作目錄
WORKDIR /app

# 複製當前目錄的內容到容器中
COPY . /app

# 啟用所有terminal功能
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    vim && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# 安裝依賴
RUN pip install --no-cache-dir -r requirements.txt

# 設置啟動時進入交互式 shell
CMD ["sh"]