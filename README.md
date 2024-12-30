# Container Usage
# Use Vscode with docker image for dev
在docker的環境下，運行terminal \ interactive jupyter
1. 啟動 Docker Desktop
2. `docker build -t [image namespace] .`建置Docker映像檔案 (通常是一次性)
3. 啟動保持容器使用狀態 
   1. `docker run -it --name [container namespace] -v $(pwd):/app [image namespace]`: 建立 container 容器, 數據持有化
   2. `docker start [containerID]`: 啟動已建立 container 容器
4. `Ctrl + Shift + P`: 輸入`>Dev Containers: Attach to Running Container`
5. 存擋後docker虛擬端與local端會資料同步


## instruction reference for working with docker container terminal
# docker (old way)
運用Dockerfile建立鏡像和容器
1. `docker build -t qpcr_app .`: 構建 Docker 映像 (通常是一次性)
2. `docker run -it --name qpcr_container -v $(pwd):/app qpcr_app`:運行容器並進入交互式終端, 持有化數據
    +	-it：啟動交互模式，允許輸入。
    +   --name qpcr_container：將容器命名為 qpcr_container。
    +   -v $(pwd):/app : 持有化數據
	+	qpcr_app：指定基於上一步構建的映像。
   1. `docker start qpcr_container`: 重啟容器
   2. `docker exec -it qpcr_container sh`: 開啟容器並進入交互終端
3. `exit`: 輸入指令退出
4. `docker start -ai qpcr_container`: 重新啟動容器
    +   start：啟動一個已停止的容器。
	+	-a（--attach）：將容器的標準輸出（STDOUT）和標準錯誤（STDERR）重新附加到當前終端，讓您可以看到容器內的輸出。
	+	-i（--interactive）：允許交互式輸入，讓您可以在容器內執行指令。
	+	qpcr_container：容器的名稱或 ID。
5. `docker start <CONTAINER_NAME or CONTAINER_ID>` 啟動停止的容器
6. `dokcer ps -a`: 查詢所有容器清單
7. `docker images`: 查詢所有印象檔案清單
8. `docker stop [containerID]`: 暫停容器
9.  `docker rm [containerID]`: 刪除容器

# docker-compose (new way)
運用yml檔案建立鏡像和容器
1. `docker-compose up --build` : 運行命令, 當相依設定更動時
2. `docker-compose up -d`: 啟動容器
3. `docker-compose up` : 運行命令
4. `docker-compose exec python_docker sh`： 進入容器的交互式終端
5. `docker-compose logs` : 檢查 Docker Compose 日誌
6. `docker ps -a` : 檢查容器狀態
7. `docker build -t myapp .` : 手動構建映像
8. `docker run -p 8000:8000 myapp` : 手動運行容器
9.  `docker pause my_container` : 暫停容器
10. `docker unpause my_container` : 暫停容器
11. `docker-compose pause` : 暫停所有服務
12. `docker-compose unpause` : 暫停所有服務
13. `docker stop my_container` : 停止容器
14. `docker-compose down` : 停止所有服務


