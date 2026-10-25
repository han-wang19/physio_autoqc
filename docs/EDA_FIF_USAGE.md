# EDA 输出与审核终版 FIF

EDA 与 ECG 使用一致的目录布局：

- `results/`：后续分析入口，包含 `final_eda_features.csv`；完成人工审阅后生成 `reviewed_fif/<记录编号>_reviewed_raw.fif`。
- `internal/`：预处理和最终 NPZ，程序内部使用。
- `audit/`：自动 QC、有效 QC、最终处理等报告。
- `plots/`：各阶段图片。

审核终版 FIF 保留完整时间轴，只有一个 `EDA` 通道，数值来自预处理 clean 波形。网页还可能显示按当前有效区间插值并重新分解的预览曲线；该预览不等同于此 FIF 的逐样本值。为保持现有数据数值和单位，使用 MNE 的 misc 通道类型；实际单位由 `raw.info['description']` 中 JSON 的 `unit` 字段记录（当前配置为微西门子），不能将数值按伏特解释。

所有不可用区间统一使用 `BAD_segment` annotation（坏段、未解决的可疑段、保留区间外和整条排除）。重叠及首尾相接的区间合并，不再重叠显示。保留区间起止时间存于 `json.loads(raw.info['description'])['analysis_retained']`，不再绘制区间标记。原始判定原因保留在审计表中。EDA 不生成 R 峰 annotation。

```python
import json
import mne

raw = mne.io.read_raw_fif('A1_music1_reviewed_raw.fif', preload=True)
eda = raw.get_data(picks=['EDA'])[0]
unit = json.loads(raw.info['description'])['unit']
bad_segments = [a for a in raw.annotations if a['description'].startswith('BAD_')]
```

尚未完成人工审阅的记录不会标为审核终版。现有特征是否通过最终计算检查，以 `audit/final_eda_processing_report.csv` 为准。
