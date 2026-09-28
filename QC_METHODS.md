# ECG 与 EDA 自动质量控制方法

本文档说明 `physio_autoqc` 当前实际实现的 QC 算法。所有阈值集中在 `project_config.py`；修改阈值后，`config_snapshot.json` 会保存该次运行的设置。

当前阈值是根据 SC26CONC 既有处理经验提出的初始项目规则，尚不是经过独立标注集验证的通用标准。正式分析前应完成盲审校准并冻结协议版本。

## 1. 共同设计

### 1.1 两级 QC

QC 分为两级：

1. **初始处理后的窗口 QC**：在完整曲目上定位信号问题，决定 PASS、TRIMMED、REVIEW 或 FAIL；
2. **最终区间复检**：对 PASS/TRIMMED 片段从源 FIF 重新裁剪和预处理，再检查一次关键质量指标。最终复检失败的片段不会进入特征表。

初始预处理产生的滤波信号、R 峰、tonic/phasic 或 SCR 只用于完整片段检查和区间选择。最终特征不通过截取这些中间结果得到。

### 1.2 窗口设置

ECG 和 EDA 默认均使用：

| 参数 | 默认值 |
|---|---:|
| 窗口长度 | 30 秒 |
| 窗口步长 | 5 秒 |
| 重叠长度 | 25 秒 |

设第一个窗口为 `[0, 30)` 秒，之后依次为 `[5, 35)`、`[10, 40)`。不足 30 秒的整个记录只产生一个短窗口。窗口的目的不是计算最终特征，而是定位局部掉线、伪迹或检测失败。

### 1.3 连续区间选择

只有状态为 PASS 的窗口参与自动区间选择。相邻 PASS 窗口合并为连续候选区间，选择持续时间最长的一个；若长度相同，保留时间更早的区间。

不会执行以下操作：

- 不把两个被坏段隔开的区间拼接；
- 不根据哪个区间的 HRV 或 EDA 结果更“好看”来选择；
- 不在 QC 后手工删除少量 RR 再假装信号连续；
- 不让 REVIEW 窗口自动进入最终区间。

## 2. ECG 初始预处理

### 2.1 信号处理

对完整 FIF 片段执行：

1. 选择配置中的 ECG 通道；
2. 保留源采样率；
3. MNE IIR 0.5–40 Hz 带通；
4. MNE IIR 50 Hz 陷波；
5. 主检测器使用 NeuroKit2 `promac`；
6. 使用 `signal_fixpeaks(iterative=True)` 修正主检测器的异常峰；
7. 在每个修正峰 ±20 ms 内搜索相对局部中位数绝对偏差最大的采样点；
8. 根据最终 R 峰计算 RR 间隔和瞬时心率。

局部细化使用绝对极值而不是单纯正向最大值，因此比旧项目更能兼容反相 ECG。

## 3. ECG 窗口 QC 指标

### 3.1 ECG 标准差

对窗口内滤波 ECG 计算：

```text
ECG_SD = std(filtered_ECG)
```

满足以下任一条件，窗口直接 FAIL：

```text
ECG_SD 不是有限数
ECG_SD <= 1e-12
```

该门槛用于识别近零、平坦或无效信号。它依赖当前数据量纲；迁移到其他设备时必须重新检查。

### 3.2 最少 R 峰数量

窗口内最少峰数根据 QC 心率范围下限计算：

```text
minimum_peaks = floor(窗口时长 × 40 bpm / 60) - 1
```

同时保证至少为 3 个峰。30 秒完整窗口默认至少需要 19 个峰。低于该数量时窗口直接 FAIL，理由为 `TOO_FEW_RPEAKS`。这项检查用于防止峰数不足时 RR 比例失去意义。

### 3.3 QRS 模板相关性

对每个主检测器 R 峰截取：

```text
R峰前 200 ms 至 R峰后 400 ms
```

只保留没有越过窗口边界的心搏。每个心搏先减去自身中位数，再除以自身标准差。随后：

1. 对全部标准化心搏逐采样点取中位数，形成中位数 QRS 模板；
2. 计算每个心搏与模板的 Pearson 相关系数；
3. 取所有有限相关系数的中位数作为窗口指标。

少于 3 个完整心搏时，该指标记为 0。

这一指标用于辅助发现噪声峰、形态不一致和错峰，但**当前不直接决定 PASS/FAIL**。原因是形态变化并不一定意味着 RR 序列无效，而且真实多形态心搏或反相波形可能使相关性降低。该值保留在窗口表和 QC 图中，供 REVIEW 复核。

### 3.4 RR 合理比例

对窗口内相邻主 R 峰计算：

```text
RR_seconds = diff(R_peak_samples) / sampling_rate
HR_bpm     = 60 / RR_seconds
```

默认宽松范围为 40–180 bpm：

```text
RR_Plausible_Fraction = 处于40–180 bpm的RR数 / 全部RR数
```

少于 3 个峰时，合理比例记为 0。

| RR 合理比例 | 解释 |
|---:|---|
| ≥0.98 | PASS 范围 |
| 0.90–<0.98 | REVIEW 范围 |
| <0.90 | FAIL 范围 |

40–180 bpm 是信号检测保护范围，不是临床正常心率定义。

### 3.5 RR 局部异常比例

当前实现还记录：

```text
局部RR中位数 = 当前RR及其前后各2个RR的中位数
异常RR = abs(RR - 局部RR中位数) > 20% × 局部RR中位数
RR_Outlier_Fraction = 异常RR数 / 全部RR数
```

当前阈值为：

| RR 异常比例 | 解释 |
|---:|---|
| ≤0.05 | PASS 范围 |
| >0.05–0.20 | REVIEW 范围 |
| >0.20 | FAIL 范围 |

局部中位数比整窗单一中位数更能容许缓慢心率趋势，但 20% 规则仍可能标记真实的快速心率变化，因此边界结果应结合 QRS 图复核。

### 3.6 RR 坏段长度与短异常修正

在 `signal_fixpeaks` 修正之前，使用初始检测峰识别不合理 RR，并把时间上连续的异常 RR 合并为坏段。每个坏段记录异常 RR 数量、开始时间、结束时间和持续时间。

当前项目规则为：

```text
异常RR数 <= 3 且坏段持续时间 <= 5秒
    → 可修正短段

异常RR数 > 3 或坏段持续时间 > 5秒
    → 不可插值长段，相应窗口FAIL

坏段持续时间 > 30秒
    → 额外标记为严重长坏段
```

短段仍由 NeuroKit `signal_fixpeaks` 完成漏峰/多峰修正，代码不会对 ECG 波形插值。最终 HRV 使用修正后的 R 峰；初始峰和修正峰均保存以便追踪。短异常占全部 RR 的比例超过 5% 时，窗口转为 REVIEW，最终区间则拒绝纳入。

输出中的 `short_bad_rr_count` 和 `short_bad_rr_fraction` 表示修正前检测到的短异常候选，并不声称每一个都被成功恢复成真实心搏。因此还需结合 `detected_rpeak_count`、`corrected_rpeak_count` 和最终波形红点检查。

## 4. ECG 窗口状态

判定顺序固定如下：

```text
若 ECG_SD 无效或 <= 1e-12
    → FAIL / FLAT_SIGNAL

否则，若 R峰数 < 根据40 bpm计算的最低峰数
    → FAIL / TOO_FEW_RPEAKS

否则，若存在超过3个连续异常RR或持续>5秒的坏段
    → FAIL / UNCORRECTABLE_RR_RUN

否则，若短异常RR比例 > 0.05
    → REVIEW / EXCESSIVE_RR_CORRECTION

否则，若 RR_plausible_fraction < 0.90
      或 RR_outlier_fraction > 0.20
    → FAIL / IMPLAUSIBLE_RR

否则，若 RR_plausible_fraction < 0.98
      或 RR_outlier_fraction > 0.05
    → REVIEW / BORDERLINE_RR

否则
    → PASS
```

采用“任一指标触发”的保守逻辑，而不是几个指标取平均。这样不会让一个很高的指标掩盖另一个明显失败的指标。

## 5. ECG 片段状态

### 5.1 PASS

同时满足：

```text
PASS窗口数 / 全部窗口数 = 1.00
完整片段主R峰数 >= 30
完整片段时长 >= 120秒
```

PASS 最终区间暂定为完整片段。随后仍从源 FIF 重新处理，而不是复用初始滤波结果。

### 5.2 TRIMMED

完整片段不满足 PASS，但存在由连续 PASS 窗口形成的最长区间，并且：

```text
最长连续区间 >= 120秒
```

该最长区间进入最终重处理。

由于窗口宽 30 秒、步长 5 秒，候选区间边界来自窗口起点与末窗口终点，不是逐采样点伪迹边界。

### 5.3 REVIEW

无法达到 PASS/TRIMMED，但窗口中同时出现 PASS 和 REVIEW 时，片段标为 REVIEW。REVIEW 不自动生成最终特征，需要后续受限人工复核。

### 5.4 FAIL

其余情况标为 FAIL，主要原因代码为：

```text
NO_SUFFICIENT_CONTINUOUS_CLEAN_INTERVAL
```

FAIL 不进入最终 ECG 处理。

## 6. ECG 最终复检

对初始状态为 PASS 或 TRIMMED 的片段：

1. 从源 FIF 重新读取 ECG；
2. 按 QC 起止时间裁剪；
3. 重新执行带通和陷波；
4. 重新运行两个 R 峰检测器；
5. 重新修正和细化主峰；
6. 对整个最终区间重新计算关键 QC；
7. 通过后才计算 HR/HRV。

最终硬性排除条件为：

| 条件 | 结果 |
|---|---|
| ECG SD 非有限或 ≤1e-12 | 排除 |
| R 峰少于 30 | 排除 |
| RR 合理比例 <0.90 | 排除 |
| RR 局部异常比例 >0.20 | 排除 |
| 仍存在不可插值 RR 坏段 | 排除 |
| 短异常 RR 比例 >0.05 | 排除 |

QRS 模板相关性仍写入最终报告，但不作为硬性排除条件。最终复检使用 RR 的 REVIEW/FAIL 分界作为硬排除门槛，而不是要求全部指标达到 PASS 分界。失败时记录为 `EXCLUDED_AFTER_FINAL_CHECK`。

### 6.1 HRV 指标资格

```text
有效时长 >= 120秒 → TimeDomain_Eligible
有效时长 >= 300秒 → FrequencyDomain_Eligible
```

当前实现通过最终复检后会计算时域 HRV；只有达到 300 秒才计算频域 HRV。后续版本还应使“不足 120 秒时不计算时域 HRV”在代码中成为严格门槛，目前自动 TRIMMED 已要求至少 120 秒，完整 PASS 通常也超过该长度。

## 7. EDA 初始预处理

对完整 FIF 片段执行：

1. 选择配置中的 GSR/EDA 通道；
2. 重采样至 100 Hz；
3. 在源采样率下保存原始有限样本掩码；
4. 在源采样率下对 NaN/Inf 做双向线性插值；
5. 将修复后的信号重采样至 100 Hz，并把有限样本掩码映射到目标时间轴；
6. 使用 NeuroKit `eda_clean(method="neurokit")` 清理；
7. 使用 `eda_phasic(method="highpass")` 分解 tonic 和 phasic；
8. 使用 `eda_peaks(method="neurokit", amplitude_min=0.01)` 检测候选 SCR；
9. 只保留 `SCR_Amplitude >= 0.05 source units` 的候选。

0.01 是检测器内部相对阈值，0.05 是当前项目的绝对 source-unit 阈值。由于没有设备校准，不能把 0.05 写成 0.05 µS。

## 8. EDA 窗口 QC 指标

EDA 窗口指标从“源采样率插值后、重采样到 100 Hz、但尚未经过 NeuroKit 清理”的信号计算；有限样本比例来自同步映射到目标时间轴的原始缺失掩码。

### 8.1 有限样本比例

```text
Finite_Fraction = 窗口内有限样本数 / 全部样本数
```

默认要求：

```text
Finite_Fraction >= 0.999
```

否则窗口 FAIL，理由为 `NONFINITE_OR_DROPOUT`。

### 8.2 信号范围

只使用有限样本计算：

```text
Signal_Range = max(raw_EDA) - min(raw_EDA)
```

若范围 ≤`1e-12`，窗口 FAIL，理由为 `FLAT_OR_NEAR_ZERO`。该阈值同样依赖源信号量纲。

### 8.3 平线比例

先计算相邻样本一阶差分 `d`，再计算稳健尺度：

```text
scale = 1.4826 × median(abs(d - median(d)))
flat_tolerance = max(机器精度, 0.01 × scale)
Flatline_Fraction = mean(abs(d) <= flat_tolerance)
```

默认：

```text
Flatline_Fraction > 0.20 → FAIL / EXCESSIVE_FLATLINE
```

这里检查的是采样点级近似不变，不是生理上的低反应。设备量化步长可能明显影响该指标，因此跨设备迁移时必须校准。

### 8.4 跳变比例

基于同一差分序列：

```text
jump_threshold = max(10 × scale, 机器精度)
异常跳变 = abs(d - median(d)) > jump_threshold
Jump_Fraction = 异常跳变数 / 全部差分数
```

默认：

```text
Jump_Fraction > 0.02 → REVIEW / EXCESSIVE_STEPS
```

跳变目前触发 REVIEW 而不是直接 FAIL，因为真实 SCR 和运动伪迹都可能产生快速变化，需要结合波形和其他指标确认。

## 9. EDA 窗口与片段状态

窗口按以下优先级判定：

```text
Finite_Fraction < 0.999
    → FAIL
否则 Signal_Range <= 1e-12
    → FAIL
否则 Flatline_Fraction > 0.20
    → FAIL
否则 Jump_Fraction > 0.02
    → REVIEW
否则
    → PASS
```

片段状态为：

- **PASS**：全部窗口均为 PASS，统一去掉首尾各 5 秒后仍不少于 120 秒；
- **TRIMMED**：完整片段不满足 PASS，但最长连续 PASS 区间去除边缘限制后仍不少于 120 秒；
- **REVIEW**：没有合格连续区间，但存在 REVIEW 窗口；
- **FAIL**：没有满足最低时长的连续干净区间，也没有可供复核的边界状态。

EDA 的 5 秒首尾裁剪用于减小清理、重采样和分解的边缘效应。

## 10. EDA 最终复检

PASS/TRIMMED 片段从源 FIF 重新执行：

```text
裁剪 → 源采样率下记录缺失并插值 → 重采样
→ EDA清理 → tonic/phasic分解
→ SCR检测 → 0.05 source-unit筛选
→ 最终特征和最终复核图
```

最终硬性排除条件为：

```text
Finite_Fraction < 0.999
或 Signal_Range <= 1e-12
```

最终报告同时保留平线比例、跳变比例和 SCR 数。当前实现没有在最终阶段因跳变比例单独排除，因为区间已经由连续 PASS 窗口选定；如果最终报告仍显示明显异常，应进入 REVIEW，而不是静默接受。

## 11. 输出与可追溯性

主要 QC 输出为：

| 文件 | 内容 |
|---|---|
| `ecg/ecg_window_qc.csv` | 每个 ECG 窗口的全部指标、状态和理由 |
| `ecg/ecg_rr_bad_runs.csv` | 每个连续异常 RR 坏段的边界、长度和是否可短段修正 |
| `ecg/ecg_qc_report.csv` | ECG 片段状态与自动选择区间 |
| `ecg/final_ecg_processing_report.csv` | ECG 最终重处理和复检结果 |
| `eda/eda_window_qc.csv` | 每个 EDA 窗口的全部指标、状态和理由 |
| `eda/eda_qc_report.csv` | EDA 片段状态与自动选择区间 |
| `eda/final_eda_processing_report.csv` | EDA 最终重处理和复检结果 |
| `config_snapshot.json` | 该次运行使用的所有配置 |

所有主要阶段均有独立图片目录。QC 结果应结合窗口表与波形图复核，不能只查看最终标签。

## 12. 当前已知限制

1. 阈值尚未通过双人盲审标注集校准；当前适合方法开发，不宜直接称为已验证标准。
2. ECG 自动 QC 依赖单一 R 峰检测流程；若主检测器系统性错峰，RR 指标仍可能受到影响，因此必须保留波形、红点和 HR 图。
3. QRS 模板相关性仅作为诊断信息，尚未形成经过验证的人工复核阈值。
4. ECG 当前记录的峰修正比例只反映峰数量变化，不反映所有峰位置移动，尚未进入决策。
5. EDA 平线和跳变标准依赖采集分辨率、单位和设备噪声。
6. EDA 没有 µS 校准，绝对 SCR 阈值只能用于相同采集流程。
7. REVIEW 的双人盲审、理由代码和裁决表尚未实现为独立界面。
8. 当前窗口区间边界精度为 5 秒；如需更精细裁剪，应在 REVIEW 阶段使用更小步长复核。

## 13. 正式使用前的校准建议

1. 从数据中抽取约 20%–30%，覆盖明显优质、边界和失败信号；
2. 隐藏被试、实验条件和自动标签，由两名审阅者独立标注；
3. 计算人工一致性，并核对每个自动指标的分布；
4. 只在校准集上调整阈值；
5. 将调整后的设置冻结为新的协议版本；
6. 用剩余数据运行一次，不再根据结果修改阈值；
7. 保存旧人工 QC 与新自动 QC 的状态转换表并进行敏感性分析。
