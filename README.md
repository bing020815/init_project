# Container Usage
## Clean cache
會徹底清掉所有 builder cache，如果 rebuild 很頻繁，每個禮拜一次清。
`docker builder prune --all` 
# Use Vscode Terminal with docker image for dev:
在docker的環境下，運行terminal \ interactive jupyter
1. 啟動 Docker Desktop
2. `docker build -t [image namespace] .`建置Docker映像檔案 (通常是一次性)
3. 啟動保持容器使用狀態 
   1. `docker run -it --name [container namespace] -v $(pwd):/app [image namespace]`: 建立 container 容器, 數據持有化
   (通常是一次性)
   2. `docker start [containerID]`: 啟動已建立 container 容器
4. `Ctrl + Shift + P`: 輸入`>Dev Containers: Attach to Running Container`
5. 存擋後docker虛擬端與local端會資料同步

# Use Vscode Ctrl + Shift + P Cmd with docker image for dev: VScode Command Palette 開發者(推薦)
在 Ctrl + Shift + P 的指令下，運行 Dev Containers command
1. 啟動 Docker Desktop
2. `Ctrl + Shift + P`: 輸入 `>Ctrl + Shift + P → Dev Containers: Open Folder in Container...`
    1. 之後只要要進去開發: 用 `>Rebuild and Reopen in Container`
    2. 若要停用: 用 `Dev Containers: Close Remote Connection`
    3. 若要啟用: 用 `Dev Containers: Open Folder in Container`

# docker-compose 相關指令
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


