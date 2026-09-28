# ECG/EDA 自动预处理与 QC

这是从现有 SC26CONC 代码中独立出的预处理工具。ECG 与 EDA 分开，预处理与 QC 也分开，总入口负责按顺序运行。

## 文件结构

| 文件 | 内容 |
|---|---|
| `project_config.py` | 每次实验主要修改这里：数据路径、实验名、通道和阈值 |
| `main.py` | 总入口 |
| `ecg_preprocess.py` | ECG 读取、滤波、R 峰检测及阶段图 |
| `ecg_qc.py` | ECG 窗口自动 QC、连续区间选择及 QC 图 |
| `ecg_finalize.py` | 从源 FIF 裁取 QC 区间，重新处理并计算最终 HR/HRV |
| `eda_preprocess.py` | EDA 重采样、清理、分解、SCR 检测及阶段图 |
| `eda_qc.py` | EDA 窗口自动 QC、连续区间选择及 QC 图 |
| `eda_finalize.py` | 从源 FIF 裁取 QC 区间，重新处理并计算最终 EDA 特征 |
| `common.py` | 共用的文件与绘图辅助函数 |

## 使用

先在 `project_config.py` 设置：

```python
EXPERIMENT_NAME = "SC26CONC"
DATA_DIR = Path(r"E:\experiment_data\SC26CONC\physiology")
FILE_PATTERN = "A*_music*_physio_raw.fif"
ECG_CHANNEL = "ECG"
EDA_CHANNEL = "GSR"
```

再运行：

```powershell
cd .\physio_autoqc
python -m pip install -r .\requirements.txt
python .\main.py
```

只处理一种信号时，把 `RUN_ECG` 或 `RUN_EDA` 设为 `False`。

## ECG 流程

```text
源FIF → 0.5–40 Hz滤波和50 Hz陷波
→ R峰检测和修正 → 滤波波形/R峰红点/心率图
→ 30秒窗口自动QC → QC时间线
→ PASS/REVIEW/TRIMMED/FAIL → 区间检查图
→ 从源FIF重新裁剪、滤波和检测 → 最终HR/HRV和最终图
```

ECG QC 以 RR 合理比例、局部 RR 异常比例和坏段长度为主。最多连续 3 个且不超过 5 秒的短异常允许由峰修正流程处理；更长坏段不插值，而是分割后选择最长连续干净区间。信号平坦度、峰数和 QRS 模板相关作为保护或诊断信息。频域 HRV 资格单独要求至少 300 秒连续数据。

## EDA 流程

```text
源FIF → 源采样率下记录缺失并插值 → 重采样到100 Hz
→ NeuroKit清理 → tonic/phasic分解 → SCR检测和绝对幅度筛选
→ 一张EDA预处理三联图
→ 30秒窗口自动QC → QC时间线
→ PASS/REVIEW/TRIMMED/FAIL → 区间检查图
→ 从源FIF重新裁剪、清理、分解和检测 → 最终EDA特征和最终图
```

现有 FIF 未保存设备到 µS 的校准，因此 EDA 幅度保持为 `source units`。

## 输出

```text
outputs/<实验名称>/
├─ config_snapshot.json
├─ ecg/
│  ├─ preprocessed/*.npz
│  ├─ ecg_preprocessing_report.csv
│  ├─ ecg_window_qc.csv
│  ├─ ecg_rr_bad_runs.csv
│  ├─ ecg_qc_report.csv
│  ├─ final/*.npz
│  ├─ final_ecg_hrv_features.csv
│  ├─ final_ecg_processing_report.csv
│  └─ plots/01_preprocess、02_qc、03_final/
└─ eda/
   ├─ preprocessed/*.npz
   ├─ eda_preprocessing_report.csv
   ├─ eda_window_qc.csv
   ├─ eda_qc_report.csv
   ├─ final/*.npz
   ├─ final_eda_features.csv
   ├─ final_eda_processing_report.csv
   └─ plots/01_preprocess、02_qc、03_final/
```

`config_snapshot.json` 保存本次运行设置。当前自动阈值是初始建议值，正式分析前仍需通过盲审样本校准并冻结版本。

每个 QC 指标的计算公式、阈值、判定顺序和已知限制见 `QC_METHODS.md`。
