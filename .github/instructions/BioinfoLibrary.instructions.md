---
applyTo: '**/*.py'
---
The instruction of developling Bioinfo Library tools.

# Bioinfo Library 工具箱開發規範 v0.1

適用範圍：使用者自建的 `Bioinfo Library` 工具模組，如 `plot.py`, `preprocess.py`, `stats.py`, `ml.py` 等

---

## 1. 命名慣例與共通變數

### 1.1 主要資料表命名
| 變數名稱           | 類型           | 說明                                           |
|--------------------|----------------|------------------------------------------------|
| `a_df`             | `pd.DataFrame` | feature 資料（如菌相 abundance table）         |
| `meta_df`          | `pd.DataFrame` | metadata 資料，與 `samples_info_df` 等價       |
| `samples_info_df`  | `pd.DataFrame` | 與 `meta_df` 相同，為樣本 metadata 表格        |

### 1.2 其他共通參數
| 變數名稱          | 類型   | 預設值 | 說明                                           |
|-------------------|--------|--------|------------------------------------------------|
| `value_col`       | `str`  | 無     | 欲分析的欄位名稱                               |
| `group_col`       | `str`  | 無     | 分組欄位                                       |
| `time_col`        | `str`  | 無     | 時間欄位                                       |
| `target_col`      | `str`  | 無     | 作為目標端點使用的欄位                         |
| `output_dir`      | `str`  | None   | 圖表或分析結果的輸出資料夾                    |
| `filename_suffix` | `str`  | ""     | 輸出檔案的後綴字串                             |
| `palette`         | `str`  | "Set2" | Seaborn 繪圖配色                               |
| `dropna_group`    | `bool` | True   | 是否排除 group 為 NA 的資料                   |
| `show_plot`       | `bool` | True   | 是否執行 plt.show()                            |
| `save_data`       | `bool` | False  | 是否儲存繪圖用的資料表為 `.csv` 備份          |

---

## 2. 裝飾器使用原則

| 裝飾器            | 說明                                                   |
|-------------------|--------------------------------------------------------|
| `@classmethod`    | 與 class 狀態相關，使用 class method 撰寫              |
| `@staticmethod`   | 與 class 狀態無關的純工具函式                           |

原則：
- 若函式需統一呼叫風格並可能擴充，使用 `@classmethod`
- 若為純工具性處理函式，使用 `@staticmethod`

---

## 3. 輸出資料與備份規範

### 3.1 資料夾處理
- 所有輸出資料夾需事先建立，避免錯誤
- 若 `output_dir` 為 None，需設預設輸出資料夾：
```python
if output_dir is None:
    output_dir = './data/stats/default'
os.makedirs(output_dir, exist_ok=True)
```

### 3.2 圖表與資料備份
- 所有圖表需輸出 `.png` 和 `.pdf`
- 若 `save_data=True`，則將繪圖所用 dataframe 輸出為 `.csv`
```python
if save_data:
    data_df.to_csv(os.path.join(output_dir, f"{filename_suffix}_data.csv"))
```

### 3.3 標準輸出 print 訊息
建議每個主流程函式加入輸出路徑提示：
```python
print(f"##### Running: {filename_suffix} at {output_dir}")
```

### 3.4 函數輸入與輸出資料型態需要明確定義
- 每個函式的參數都要有明確的資料型態定義
- 每個class底下的函式都要有預設的檔案格式，像是pd.DataFrame, str, int等。
```python
def addition(variable1: int, variable2: int) -> int:
    return variable1 + variable2
```

### 3.5 圖表函數輸入會將資料與meta拆開
- 菌的資料會有自己的dataframe，用`OG_sample`作為index。
- 路徑的資料會有自己的dataframe，用`OG_sample`作為index。
- metadat也有自己的dataframe，含有組別、血檢資料以及其他補充資料，用`OG_sample`作為index。

---

## 4. 統計分析使用習慣

### 4.1 Spearman correlation 檢定
- 統一使用既有共用函式：`StatsBioData.calculate_spearman()`
- 禁止於主流程中重複撰寫計算邏輯
- 建議結果包含：rho, p-value, 有效樣本數（n）

---

## 5. 視覺化繪圖規則

### 5.1 圖表類型
- 折線圖：`plot_linechart_by_group()`
- Boxplot：`plot_boxplot_by_group()`

### 5.2 標準設定
```python
sns.despine()
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
```

---

## 6. 資料預處理規則（`preprocess.py`）

### 函式命名與功能
| 函式名稱              | 說明                             |
|-----------------------|----------------------------------|
| `tabular_to_long()`   | wide 轉 long 格式               |
| `long_to_tabular()`   | long 轉 wide 格式               |
| `merge_df_by_index()` | 依 index 合併多個 dataframe     |

原則：
- 所有資料以 `OG_SampleID` 為 index，確保可對齊樣本的ID與metadata結合使用。
- 所有資料需包含 `OG_SampleID` (SampleID在不同run) 欄位，以便於後續分析與比對。
- `SampleID` 必須為唯一值，避免重複。可對應最後選用的樣本 ID(OG_SampleID)。


---

## 7. 機器學習模組（`ml.py`）規則

- 支援 XGBoost、LightGBM、SHAP 等分類器
- 所有模型輸出需標示版本（使用 `note` 或 timestamp）
- 輸出結果包括：模型檔、重要性排序圖、預測結果表格等

---

## 8. 專案資料夾結構建議

```
project_root/
├── data/
│   ├── raw/
│   ├── preprocess/
│   ├── stats/
│   ├── ml/
│   └── plots/
├── app/
│   └── subproject1/
│   └── subproject2/
└── docs
```

---

## 9. 建議擴充項目（未來版本）
- 待補齊功能