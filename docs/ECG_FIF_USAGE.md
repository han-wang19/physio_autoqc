# 审核后 ECG：后续分析入口

完成每条 ECG 审阅后自动生成 `outputs/<实验名称>/ecg/results/reviewed_fif/<记录编号>_reviewed_raw.fif`。结束本次审核并更新汇总也会同步导出。文件保存完整记录时间轴，不拼接或删除坏段。

- `ECG`：审核页面使用的预处理心电波形，单位 V。
- 博睿康原始 `uV/mm` 电压数据按厂商兼容规则换算至 V；换算来源记录在 FIF 的 `info['description']`。旧缓存不会被重复改写，导出会识别其校准标记。依据见 [单位核查](ECG_UNIT_CHECK.md)。
- 文件只有 `ECG` 一个通道。`R_PEAK` annotation：最终人工确认的 R 峰，onset 为峰时间，duration 为 0。包含坏段内的峰，提取特征时必须结合坏段筛选。
- `raw.annotations`：所有不可用区间统一为 `BAD_segment`，包括坏段、未解决的可疑段、分析区间外和整条排除。重叠及首尾相接的区间取并集，因此坏段互不重叠。原始原因保留在 `audit/` 报告和实验级人工决定表中。
- 保留区间不再绘制 annotation；起止秒数存于 `json.loads(raw.info['description'])['analysis_retained']`，不适用时为 `None`。

这是审核结果交换文件，不代表记录已通过后续 HRV 最终质量检查。是否被现有 HRV 流程纳入，以 `audit/final_ecg_processing_report.csv` 为准。

```python
import mne
import numpy as np

raw = mne.io.read_raw_fif('A1_music1_reviewed_raw.fif', preload=True)
sfreq = raw.info['sfreq']
ecg = raw.get_data(picks=['ECG'])[0]
# 将零时长 R_PEAK annotation 转回采样点索引。
peak_mask = raw.annotations.description == 'R_PEAK'
peaks = np.rint(raw.annotations.onset[peak_mask] * sfreq).astype(int)
valid = np.ones(raw.n_times, dtype=bool)
for onset, duration, label in zip(raw.annotations.onset, raw.annotations.duration,
                                  raw.annotations.description):
    if label.startswith('BAD_'):
        left = max(0, int(round(onset * sfreq)))
        right = min(raw.n_times, int(round((onset + duration) * sfreq)))
        valid[left:right] = False
valid_peaks = peaks[valid[peaks]]
# RR 必须来自原始相邻峰，且两峰之间没有任何坏样本；不要跨坏段计算 RR。
rr_ms = [(b-a) / sfreq * 1000 for a, b in zip(peaks[:-1], peaks[1:])
         if valid[a:b+1].all()]
```

其他目录中的 NPZ 是程序内部缓存及历史计算产物。使用上述 FIF 进行新分析无需手动读取它们；现有审核和 HRV 程序仍依赖这些文件，因此暂不删除。
