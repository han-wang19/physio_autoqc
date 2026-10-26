# ECG 与 EDA 自动质量控制方法

本文档说明 `physio_autoqc` 当前实际实现的 QC 算法。所有阈值集中在 `project_config.py`；修改阈值后，`config_snapshot.json` 会保存该次运行的设置。

**实现同步状态（2026-10-25）**：已核对实验配置、预处理、QC、最终处理、人工审核和 FIF 导出。方法阈值为项目规则，不能把工程回归检查视为科学有效性验证。

### 实验与条件适用范围

The default inclusion floor is 120 seconds. Private configuration can set a condition-specific floor, e.g. 60 seconds for a short baseline. Such overrides do not establish validity of every HRV measure; standard 5-minute and frequency-domain eligibility remain separate.

每条规则在本文中标为三类：**指南支持**、**文献启发并按设备适配**、**项目待验证**。后两类不能表述为通用标准，正式分析前仍应完成盲审校准并冻结协议版本。

## 0. Presentation 级方法汇总

### 总体原则

```text
源 FIF → 标准化预处理 → 局部窗口 QC → 记录级处置
       → 人工波形复核 → 从源 FIF 重处理 → 最终复检 → 特征
```

- 自动 QC 用于定位问题和筛查，不代替波形复核；
- 最终特征从源 FIF 重新处理，不直接截取开发阶段中间数组；
- ECG 不跨坏段拼接 RR；EDA 允许为保持分解连续性而插值，但插值样本不进入统计；
- 所有自动阈值、人工裁决、峰位修改和区间标注均保留可追溯输出。

### ECG：主要步骤与依据

| 步骤 | Presentation 中的通俗表述 | 主要依据 |
|---|---|---|
| 波形清理与心搏定位 | 先去除缓慢漂移、高频噪声和 50 Hz 工频干扰，再结合多种算法找出每次心搏。最后把标记统一对齐到每个 QRS 波的同一个位置，避免在正、负尖峰之间跳动（0.5–40 Hz；50 Hz 陷波；局部对齐 ±50 ms）。 | Task Force (1996)；Laborde et al. (2017)；Makowski et al. (2021)；精确滤波参数为项目实现，自动检测必须配合波形复核 |
| 逐段质量检查 | 用滑动时间窗逐段查看信号是否平坦或噪声过大、心搏是否漏检/多检、心率与相邻心搏间隔是否可信。每段标为 PASS、REVIEW 或 FAIL（30 s 窗，每 5 s 向前移动）。 | Zhao & Zhang (2018)；Clifford et al. (2012)；窗长与数值阈值属项目待验证 |
| 选择可用时段 | 如果整首曲子都合格，保留全段；如果只有局部问题，保留最长的连续干净段。不把被坏信号隔开的两段 ECG 硬拼在一起（最低 120 s）。 | HRV 必须使用连续、经校验的 NN 间期；Task Force (1996) |
| 最终 HRV | 对选中时段从原始文件重新处理和检查，通过后才计算 HRV。不足 5 分钟的结果只作探索性时域分析；达到 5 分钟才标为标准短时 HRV，而频域分析还要求心搏修正负担很低（120–299 s：探索性；≥300 s；频域修正比例 ≤2%）。 | Task Force (1996)；Laborde et al. (2017)；APS (2026)；频域 2% 为保守项目规则 |

### EDA：主要步骤与依据

| 步骤 | Presentation 中的通俗表述 | 主要依据 |
|---|---|---|
| 统一单位并分离慢/快变化 | 先把皮肤电导换算为统一的 µS，保留原始缺失位置，再清理波形。随后将信号分为缓慢变化的 tonic 和短暂反应的 phasic（统一为 100 Hz）。 | Boucsein et al. (2012)；单位、采样率、设备与处理方法必须报告 |
| 识别皮电反应 | 在快变化成分中寻找短暂的皮电反应（SCR），只保留幅度足够、不太可能是微小噪声的反应（≥0.05 µS）。 | SPR 建议报告 0.01–0.05 µS 范围内的最小幅度标准 |
| 逐段质量检查 | 每 5 秒检查一次：是否掉线、超出设备范围，或出现过快而不像生理变化的跳变。硬伪迹及其前后 5 秒不纳入；低电导或很平稳的信号只提醒人工复核，不自动删除（硬跳变 >10 µS/s）。 | Taylor et al. (2015)；Kleckner et al. (2018)；低反应不应自动等同于伪迹 |
| 处理局部坏段 | 保留整首曲子的原始时间轴，不因中间一个短伪迹删掉前后所有数据。坏段只为了让滤波和信号分解能连续运行而临时补齐；补齐值不进入特征或 SCR 统计（实际好数据合计至少 120 s）。 | Boucsein et al. (2012)；Kelsey et al. (2023) |
| 最终特征与复检 | 在整段信号上完成清理和分解后，只用标记为有效的真实样本计算 tonic、phasic 和 SCR 特征，并再检查这些样本是否超量程或仍有非生理跳变。 | SPR 可追溯报告原则；Kleckner et al. (2018) |

> Presentation 中应将“指南支持”与“项目阈值”分开表述：文献支持 QC 框架和报告原则，但不意味着 30 s 窗、20% RR 偏差或 120 s 最低时长是通用金标准。

## 1. 共同设计

### 1.1 两级 QC

QC 分为两级：

1. **初始处理后的窗口 QC**：在完整曲目上定位信号问题，决定 PASS、TRIMMED、REVIEW 或 FAIL；
2. **最终区间复检**：PASS/TRIMMED 记录从源 FIF 重新处理。ECG 先裁切选定的连续区间再滤波；EDA 先在完整时间轴上完成修复和分解，再裁切记录级分析区间。两者都再检查一次关键质量指标；失败时不进入特征表。

初始预处理产生的滤波信号、R 峰、tonic/phasic 或 SCR 只用于完整片段检查和区间选择。最终特征不通过截取这些中间结果得到。

### 1.2 窗口设置

两种信号使用不同窗口：

| 信号 | 窗口 | 步长 | 依据 |
|---|---:|---:|---|
| ECG | 30 秒 | 5 秒 | **项目待验证**：兼顾 RR 稳定性与坏段定位 |
| EDA | 5 秒 | 5 秒 | **文献启发**：EDA 专家标注和自动伪迹研究常以 5 秒 epoch 评估 |

窗口只用于定位掉线、伪迹或检测失败，不直接计算最终特征。ECG 使用重叠窗口；EDA 使用不重叠 5 秒窗口，并在硬伪迹前后增加 5 秒过渡缓冲。

### 1.3 区间选择

ECG 只使用 PASS 窗口：相邻 PASS 窗口合并为连续候选区间，选择持续时间最长的一个；若长度相同，保留更早的区间。EDA 则保留首尾各去 5 s 后的完整时间轴，汇总所有 PASS 窗口的有效时长；局部非 PASS 样本通过 `valid_mask` 排除。

不会执行以下操作：

- ECG 不把两个被坏段隔开的区间拼接；EDA 虽可同时汇总多个有效段，但保留原时间坐标，不把坏段两端当作相邻样本；
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
7. 使用全实验共用的双相形态/主偏转规则，把修正峰对齐到 ±50 ms 内一致的定位点；稳定的前正后负 RS 形态优先取前面的明显正峰，其余保留主偏转判断，详见“双相 QRS 的统一对齐规则”；
8. 根据最终 R 峰计算 RR 间隔和瞬时心率。

`promac` 会组合多个成员检测器。个别成员（例如 Gamboa 或 Zong）不适用于某段数据时，NeuroKit 会在内部跳过它们；项目隐藏这类已被 `promac` 捕获的打印信息，但不会吞掉整个峰检测流程抛出的异常。

局部细化先确定记录级 QRS 极性，而不是固定使用正向最大值或逐搏选择绝对极值；这样既兼容反相 ECG，又避免双相 QRS 在正、负偏转之间跳变，保证 RR 使用一致的 fiducial point。

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

### 3.4 多指标 ECG 信号质量（Zhao–Zhang 2018）

每个 30 秒窗口对滤波 ECG 应用 Zhao–Zhang 单导联 SQI 的三指标启发式融合，综合：

- QRS 频带功率比例 `pSQI`；
- 峰度 `kSQI`；
- 基线相对功率 `basSQI`。

当前实现不包含原论文的双检测器一致性 `qSQI`，因此属于三指标适配，而非原论文完整方法。功率指标从同一个 Welch PSD 积分计算，`basSQI = 1 - P(0–1 Hz)/P(0–40 Hz)`，即先按 0–40 Hz 总功率归一化，避免 ECG 单位/整体增益改变基线分数；峰度采用 Pearson 定义（正态分布为 3），而不是 Fisher excess kurtosis。阈值及 simple-fusion 规则沿用 Zhao–Zhang/NeuroKit Zhao 规则：`pSQI` 按最高心率分段，`kSQI > 5`，`basSQI >= .95` 为 optimal、`.90–.95` 为 suspicious、低于 `.90` 为 unqualified；至少 2 项 unqualified（或 1 项 unqualified 且另 2 项 suspicious）判 Unacceptable，至少 2 项 optimal 且无 unqualified 判 Excellent，其余为 Barely acceptable。该实现和版本应随结果报告，且不能称为原论文完整验证方法。结果映射为：

```text
Excellent          → 继续检查 RR 指标
Barely acceptable  → REVIEW / BORDERLINE_SIGNAL_QUALITY
Unacceptable       → FAIL / UNACCEPTABLE_SIGNAL_QUALITY
```

“Barely acceptable”只触发复核而不直接排除，符合原论文要求重新评估边界质量的语义，也避免把特定设备或 QRS 形态差异误当伪迹。该 SQI 补足了仅凭修正后 RR 序列可能漏过规律性假峰或高频噪声的问题。依据：Zhao & Zhang (2018) 以及使用多 SQI、专家标注验证 ECG 可解释性的 Clifford et al. (2012)。

### 3.5 RR 合理比例

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

### 3.6 RR 局部异常比例

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

### 3.7 RR 坏段长度与短异常修正

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

### 双相 QRS 的统一对齐规则

自动峰精修使用所有实验共用的 `consistent_rs_v2` 规则，无被试或记录专属极性例外。先在候选峰 ±50 ms 内统计正、负偏转及其顺序；当至少 90% 心搏具有前正后负的稳定形态、正向偏转至少为负向偏转的 50%，且两者间隔为 8–50 ms 时，对齐前面的正峰。其余情况沿用整条记录的主偏转极性，避免强制把负向 QRS 改为正峰。相应阈值集中配置于 `project_config.py`。

This shared morphology heuristic is not medical Q/R/S identification. Complex shapes require manual review. Existing cached peaks must be explicitly reprocessed to apply a newer alignment version; manual edits remain authoritative.

## 4. ECG 窗口状态

判定顺序固定如下：

```text
若 ECG_SD 无效或 <= 1e-12
    → FAIL / FLAT_SIGNAL

否则，若 R峰数 < 根据40 bpm计算的最低峰数
    → FAIL / TOO_FEW_RPEAKS

否则，若存在超过3个连续异常RR或持续>5秒的坏段
    → FAIL / UNCORRECTABLE_RR_RUN

否则，若 RR_plausible_fraction < 0.90
      或 RR_outlier_fraction > 0.20
    → FAIL / IMPLAUSIBLE_RR

否则，若 Zhao–Zhang SQI = Unacceptable
    → FAIL / UNACCEPTABLE_SIGNAL_QUALITY

否则，若短异常RR比例 > 0.05
    → REVIEW / EXCESSIVE_RR_CORRECTION

否则，若 RR_plausible_fraction < 0.98
      或 RR_outlier_fraction > 0.05
    → REVIEW / BORDERLINE_RR

否则，若 Zhao–Zhang SQI = Barely acceptable
    → REVIEW / BORDERLINE_SIGNAL_QUALITY

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

无法达到 PASS/TRIMMED，但仍有 REVIEW 窗口，且不是“0 个 PASS 并已出现 FAIL”的情形时，片段标为 REVIEW。REVIEW 不自动生成最终特征，需要后续受限人工复核。

### 5.4 FAIL

若没有任何 PASS 窗且至少有一个 FAIL 窗，直接标为 `FAIL / NO_PASS_WINDOWS`。其余无法获得足够连续干净区间的情况标为：

```text
NO_SUFFICIENT_CONTINUOUS_CLEAN_INTERVAL
```

FAIL 不进入最终 ECG 处理。

## 6. ECG 最终复检

对初始状态为 PASS 或 TRIMMED 的片段：

1. 从源 FIF 重新读取 ECG；
2. 按 QC 起止时间裁剪；
3. 重新执行带通和陷波；
4. 重新运行配置的 R 峰检测流程；
5. 重新修正和细化主峰；
6. 对整个最终区间重新计算关键 QC；
7. 通过后才计算 HR/HRV。

最终复检同时保存两类 RR 坏段信息：检测器初始峰上的短异常候选比例（用于报告自动修正负担和频域资格），以及最终实际使用峰上的残留异常（用于判断是否仍有不可接受坏段）。人工峰替换自动峰时，残留异常以人工峰计算；修正比例以人工修改峰数/自动峰数计算。人工峰已受审阅者确认，因此不再用“初始检测器短异常 >5%”单独否决，但仍必须通过最终峰的 RR 合理性和坏段检查。

最终硬性排除条件为：

| 条件 | 结果 |
|---|---|
| ECG SD 非有限或 ≤1e-12 | 排除 |
| R 峰少于 30 | 排除 |
| RR 合理比例 <0.90 | 排除 |
| RR 局部异常比例 >0.20 | 排除 |
| Zhao–Zhang 三指标 SQI = Unacceptable | 排除 |
| 最终使用峰仍存在不可插值 RR 坏段 | 排除 |
| 无人工峰时，自动检测短异常 RR 比例 >0.05 | 排除 |

QRS 模板相关性仍写入最终报告，但不作为硬性排除条件。最终复检使用 RR 的 REVIEW/FAIL 分界作为硬排除门槛，而不是要求全部指标达到 PASS 分界。失败时记录为 `EXCLUDED_AFTER_FINAL_CHECK`。

### 6.1 HRV 指标资格

```text
有效时长 >= 120秒 → TimeDomain_Eligible（探索性短时结果）
有效时长 >= 300秒 → Standard_5min_HRV_Eligible
有效时长 >= 300秒 且短异常RR比例 <= 0.02 → FrequencyDomain_Eligible
```

**指南支持**：标准短时 HRV 优先使用至少 5 分钟且不同条件采用相同长度。当前项目为兼容短曲目，120–299 秒仅保留探索性时域指标并单独标记；不与标准 5 分钟结果混称。频域指标对伪迹更敏感，因此额外要求修正候选不超过 2%；2% 是依据伪迹敏感性实验设置的保守项目门槛，不是国际统一阈值。

## 7. EDA 初始预处理

### 7.1 单位校准

以下为设备文件头的格式示例，不对应仓库中提供的任何参与者记录：

```text
Dimension = uS
Physical range = 0–100
Samples per record = 100（即100 Hz）
```

程序依据 `EDA_INPUT_UNIT` 查表换算到 µS（支持 S、mS、uS/µS/μS、nS；未知单位会报错）。本数据集输入单位为 `uS`，因此自动得到 `EDA_SCALE_TO_US = 1.0`。校准结果、输入单位、换算系数和来源写入 NPZ 与 CSV。这属于**从原始元数据恢复单位并自动应用换算**，不是根据数据分布猜测单位。更换数据集时必须重新检查设备或原始文件头；若无元数据、设备传递函数或已知校准电阻，仅凭波形无法唯一确定换算比例。

### 7.2 处理步骤

对完整 FIF 片段执行：

1. 选择配置中的 GSR/EDA 通道并按元数据换算为 µS；
2. 在源采样率下保存原始有限样本掩码；
3. 在源采样率下对 NaN/Inf 做双向线性插值；
4. 重采样至 100 Hz，并把有限样本掩码映射到目标时间轴；
5. 使用 NeuroKit `eda_clean(method="neurokit")` 清理；
6. 使用 `eda_phasic(method="highpass")` 分解 tonic 和 phasic；
7. 使用 `eda_peaks(method="neurokit", amplitude_min=0.01)` 检测候选 SCR；
8. 只保留 `SCR_Amplitude >= 0.05 µS` 的候选。

**指南支持**：SPR 建议报告 SCR 最小幅度标准，常用范围为 0.01–0.05 µS；本项目选择较保守的 0.05 µS。100 Hz 对本数据不构成升采样，因为源 GSR 本身为 100 Hz。

## 8. EDA 窗口 QC 指标

EDA 指标从换算为 µS、修复缺失并重采样、但尚未经过 NeuroKit 清理的信号计算；有限样本比例仍来自插值前缺失掩码。

### 8.1 有限样本比例

```text
Finite_Fraction = 窗口内有限样本数 / 全部样本数
```

默认要求：

```text
Finite_Fraction >= 0.999
```

否则窗口 FAIL，理由为 `NONFINITE_OR_DROPOUT`。

### 8.2 设备量程与极低电导

检查每个样本是否位于：

```text
0 <= EDA <= 100 µS
```

只要存在低于 0 或高于 100 µS 的样本，窗口 FAIL，理由为 `OUTSIDE_DEVICE_RANGE`。两端依据本项目 BDF 物理量程。低于 0.05 µS 的电导比例另行报告；若至少 95% 样本低于此值，标记 `REVIEW / VERY_LOW_CONDUCTANCE`，不作为硬件超量程失败。0.05 µS 是特定 EDA QA 研究验证数据的范围下沿，不能据此断言其他设备或被试的低值无效。

### 8.3 最小可辨识窗口变化

每个 5 秒窗口计算原始 EDA 的峰峰值：

```text
Signal_Range = max(EDA) - min(EDA)
Signal_Range < 0.01 µS → REVIEW / BELOW_MIN_MEANINGFUL_VARIATION
```

SPR 指南给出的可计数 SCR 最低幅度范围为 0.01–0.05 µS，并特别要求结合设备噪声，避免把放大器噪声当成生理反应。这里采用其保守下限 `0.01 µS` 作为“窗口内没有达到文献所述最小可辨识变化”的**复核触发值**，而不是通用的平线或伪迹判据。低反应可能是真实生理状态，所以该规则只产生 REVIEW；记录是否可用仍由连续 PASS 时长和人工复核决定。更换设备或量程时必须重新验证此映射。

该补充规则用于捕获平滑插值或放大器微噪声造成的伪变化：此类数据的每个相邻样本可能都不相同，因而仅靠“完全相等样本比例”会漏检。

### 8.4 平线比例

先计算相邻样本一阶差分 `d`，再计算稳健尺度：

```text
scale = 1.4826 × median(abs(d - median(d)))
flat_tolerance = max(机器精度, 0.01 × scale)
Flatline_Fraction = mean(abs(d) <= flat_tolerance)
```

默认仅作为 REVIEW：

```text
Flatline_Fraction > 0.95 → REVIEW / POSSIBLE_FLATLINE_OR_QUANTIZATION
```

低反应或短时间稳定 EDA 可能是真实生理现象，因此不再因平线比例直接 FAIL。该 0.95 阈值是**项目待验证**的诊断警报。

### 8.5 最大绝对斜率

对相邻 µS 样本计算：

```text
Slope = diff(EDA) × sampling_rate
Max_Abs_Slope = max(abs(Slope))
```

若：

```text
Max_Abs_Slope > 10 µS/s → FAIL / NONPHYSIOLOGICAL_SLOPE
```

该规则来自经专家评价的 ambulatory EDA QA 流程，属于**文献启发并按本设备适配**。每个硬伪迹窗口前后再扩展 5 秒并标记 `ADJACENT_TO_ARTIFACT`，用于隔离电极接触变化的过渡段。

## 9. EDA 窗口与片段状态

窗口按以下优先级判定：

```text
Finite_Fraction < 0.999
    → FAIL
否则存在 <0 或 >100 µS 的样本
    → FAIL
否则 Max_Abs_Slope > 10 µS/s
    → FAIL
否则至少 95% 样本 <0.05 µS
    → REVIEW
否则 Signal_Range < 0.01 µS
    → REVIEW
否则 Flatline_Fraction > 0.95
    → REVIEW
否则
    → PASS
```

片段状态为：

- **PASS**：全部窗口均为 PASS，统一去掉首尾各 5 秒后仍不少于 120 秒；
- **TRIMMED**：完整片段不满足 PASS，但首尾各去 5 秒后，所有 PASS 窗口的总有效时长仍不少于 120 秒；分析区间保留完整原时间轴，局部坏段由 `valid_mask` 排除；
- **REVIEW**：有效总时长不足 120 秒，但存在 REVIEW 窗口；
- **FAIL**：有效总时长不足 120 秒，且没有可供复核的 REVIEW 窗口。

EDA 的 5 秒首尾裁剪用于减小清理、重采样和分解的边缘效应。

## 10. EDA 最终复检

PASS/TRIMMED 片段从源 FIF 重新执行：

```text
完整原时间轴上记录缺失并重采样 → 按 PASS 窗建立 valid_mask
→ 仅为连续滤波/分解而插值坏段 → EDA清理 → tonic/phasic分解
→ SCR检测 → 0.05 µS筛选
→ 按自动或人工记录级区间裁切并用 valid_mask 筛选 SCR
→ 最终特征和最终复核图
```

最终硬性排除条件为：

```text
经 valid_mask 筛选后的有效总时长 < 120 s
或最终有效样本仍存在设备量程外值
或最终相邻有效样本的 Max_Abs_Slope > 10 µS/s
```

最终报告保留有限比例、设备量程越界比例、低电导比例、平线比例、最大斜率和 SCR 数。低电导及平线警报均为 REVIEW，不自动排除真实低反应信号。

### 10.1 EDA 局部伪迹与非连续有效数据

EDA 不再因中间出现短伪迹而只保留最长连续 PASS 段。完成统一的首尾 5 秒裁剪后，
记录级判定使用全部 PASS 窗口的总有效时长：全部窗口通过时标为 PASS；存在局部坏窗、
但 PASS 窗口总时长不少于 120 秒时标为 TRIMMED，并保留完整原时间轴。

最终处理从 `eda_window_qc.csv` 建立逐采样点 `valid_mask`。非 PASS 窗口只为滤波和
tonic/phasic 分解进行线性修复，修复值本身不进入统计：

整条记录先完成修复、clean 和 tonic/phasic 分解，之后才按人工或自动分析区间裁切，
以避免裁切边界产生新的滤波差异，并保证 QC 图的 repaired-clean 与 final 图逐点一致。

- tonic、phasic 和 raw 汇总只使用 `valid_mask=True` 的样本；
- tonic slope 使用有效样本的原始时间坐标，不把分离的好段首尾拼接；
- SCR 峰位于坏窗或缓冲区时删除，SCR rate 的分母为有效时长而非总经过时长；
- 最终 NPZ 同时保存 `raw`、`repaired_raw` 和 `valid_mask`；报告保存有效时长、有效比例
  及为分解而修复的比例。

该策略依据 SPR 建议的“检查并排除包含伪迹的 portions/segments”，避免将一个局部机械
伪迹扩大为对更早或更晚全部有效数据的删除。插值仅用于使滤波和分解保持连续，不把插值
片段当作真实观测。对于必须连续的时序分析，仍应按 `valid_mask` 分段，不能跨坏段拼接。

## 11. 输出与可追溯性

人工点击“完成审阅”后，程序使用 `review/quality.py` 的同一套采样点区间解析规则生成图中底色、人工记录级决定和审核 FIF；不因有没有新增人工区间、有没有修改 R 峰，或原自动标签是 PASS/TRIMMED/REVIEW/FAIL 而切换判定路径。人工结果可以与保留作审计证据的原自动结果不同。

- 全程有效且满足本条件最低时长：PASS。
- ECG 仅部分有效：取最长连续 PASS 区间，满足最低时长则为 TRIMMED，不拼接不连续段。
- EDA 先剔除固定边缘，再按有效样本总时长判断；固定边缘内侧全有效时为 PASS，其余满足时长的情况为 TRIMMED。
- 有效时长不足且仍有待复核区间：REVIEW；有效时长不足且没有待复核区间：FAIL。

先成功写出审核 FIF，再保存完成状态。历史版本的已完成决定不会因为页面刷新而静默改写；历史一致性修复应保留旧决定与 FIF 备份。

“结束本次审核并更新汇总”按记录增量处理。检查点包含审核时间、输入文件状态、窗口与人工区间、处理参数和处理代码的签名；未变更且产物完整的记录跳过重算。特征表先写入、处理报告后写入，并用临时文件替换，避免中断时将未保存的特征误认为已完成。旧版本中审核时间匹配且产物完整的检查点可直接继承。页面显示当前记录、完成数量和阶段；缺少产物或输入发生变化时会重新处理。

主要 QC 输出为：

| 文件 | 内容 |
|---|---|
| `ecg/audit/ecg_window_qc.csv` | 每个 ECG 窗口的全部指标、状态和理由 |
| `ecg/audit/ecg_rr_bad_runs.csv` | 每个连续异常 RR 坏段的边界、长度和是否可短段修正 |
| `ecg/audit/ecg_qc_report.csv` | ECG 片段状态与自动选择区间 |
| `ecg/audit/final_ecg_processing_report.csv` | ECG 最终重处理和复检结果 |
| `eda/audit/eda_window_qc.csv` | 每个 EDA 窗口的全部指标、状态和理由 |
| `eda/audit/eda_qc_report.csv` | EDA 片段状态与自动选择区间 |
| `eda/audit/final_eda_processing_report.csv` | EDA 最终重处理和复检结果 |
| `manual_qc_decisions.csv` | 记录级人工纳入、裁切、排除或二审决定 |
| `manual_qc_intervals.csv` | 波形图上的多区间 PASS/REVIEW/FAIL 人工标注（审计与复核用） |
| `ecg/audit/ecg_qc_effective_summary.csv` / `eda/audit/eda_qc_effective_summary.csv` | 合并记录级人工裁决并统计人工区间数量/时长的有效汇总 |
| `config_snapshot.json` | 该次运行使用的所有配置 |

`plots/01_qc` 合并展示预处理波形、事件标记与窗口状态，`plots/02_final` 展示最终纳入区间；不再重复保存信息重叠的预处理图片。QC 结果应结合窗口表与波形图复核，不能只查看最终标签。

## 12. 当前已知限制

1. 阈值尚未通过双人盲审标注集校准；当前适合方法开发，不宜直接称为已验证标准。
2. ECG 主峰仍来自单一 ProMAC 集成流程；Zhao–Zhang SQI 增加了频谱、峰度和基线质量，但 NeuroKit2 版本未实现原论文的双检测器一致性 `qSQI`，因此必须保留波形、红点和 HR 图供人工确认。
3. QRS 模板相关性仅作为诊断信息，尚未形成经过验证的人工复核阈值。
4. ECG 当前记录的峰修正比例只反映峰数量变化，不反映所有峰位置移动，尚未进入决策。
5. EDA 输入单位和校准来源已按实验登记；新增实验仍需核对原始文件头，不能只改实验名称。10 µS/s 及低电导 REVIEW 门槛仍需在各实验盲审验证。
6. 换用其他数据集时不能沿用 `EDA_SCALE_TO_US=1.0`，必须重新核对设备或原始文件头。
7. 当前界面保存审阅者和记录级决定；没有独立的第二审阅者盲法与冲突仲裁工作流，不能把 SECOND_REVIEW 状态视为已完成双人复核。
8. 自动窗口步长为 5 秒；人工拖选按波形采样点附近的时间标注，有效区间和 FIF 坏段按实际边界导出，不强行限制为 5 秒。
9. 网页默认拖选时间区间后选择 PASS/FAIL/取消；区间仅由横坐标确定。人工区间按最新标注覆盖旧标注，并替换相应底色。ECG 开启 R 峰编辑后，悬浮显示当前横坐标对应波形样本的候选红点，单击添加、双击已有峰删除；此交互添加精确样本，不再自动吸附局部极值。区间与峰共享撤销历史。修改暂存在会话中，统一点击“完成审阅”保存区间、峰和完成状态；未修改记录也通过同一按钮确认。审核完成状态与 PASS/TRIMMED/FAIL/REVIEW 的数据质量结论分开，已完成但仍为 REVIEW 的记录不会被自动纳入。ECG 仍只使用最长连续 PASS 段，不拼接不连续段；EDA 保留完整时间轴并覆盖 `valid_mask`。最终特征仍需运行 finalize；汇总通过“结束本次审核并更新汇总”统一刷新。

## 13. 正式使用前的校准建议

1. 从数据中抽取约 20%–30%，覆盖明显优质、边界和失败信号；
2. 隐藏被试、实验条件和自动标签，由两名审阅者独立标注；
3. 计算人工一致性，并核对每个自动指标的分布；
4. 只在校准集上调整阈值；
5. 将调整后的设置冻结为新的协议版本；
6. 用剩余数据运行一次，不再根据结果修改阈值；
7. 保存旧人工 QC 与新自动 QC 的状态转换表并进行敏感性分析。

## 14. 文献依据

1. Task Force of the European Society of Cardiology and the North American Society of Pacing and Electrophysiology. *Heart rate variability: standards of measurement, physiological interpretation and clinical use* (1996). 标准短时 HRV、频域定义与记录长度基础。<https://www.escardio.org/static-file/Escardio/Guidelines/Scientific-Statements/guidelines-Heart-Rate-Variability-FT-1996.pdf>
2. American Physiological Society. *Guidelines for rigor and reproducibility of heart rate variability within human cardiovascular research* (2026). 支持至少 5 分钟、等长 epoch、自动检测结合人工复核，并报告 RR 删除/修正比例与持续时间。<https://journals.physiology.org/doi/full/10.1152/ajpheart.00041.2026>
3. Laborde S, Mosley E, Thayer JF. *Heart Rate Variability and Cardiac Vagal Tone in Psychophysiological Research* (2017). 支持 5 分钟优先、ECG/RR 人工或软件辅助伪迹校正及完整报告。<https://doi.org/10.3389/fpsyg.2017.00213>
4. *Heart Rate Variability Analysis: How Much Artifact Can We Remove?* (2020). 频域指标对少量 RR 删除或错位特别敏感，因此本项目将 2% 作为保守频域资格门槛，但该数值仍标为项目规则。<https://pubmed.ncbi.nlm.nih.gov/33017533/>
5. Society for Psychophysiological Research Ad Hoc Committee on Electrodermal Measures. *Publication recommendations for electrodermal measurements* (2012). 支持报告设备、校准、采样率、处理方法与 SCR 最小幅度；常见最小幅度为 0.01–0.05 µS。<https://doi.org/10.1111/j.1469-8986.2012.01384.x>
6. Taylor S et al. *Automatic Identification of Artifacts in Electrodermal Activity Data* (2015). 使用专家标注 5 秒 epoch，并强调 EDA 形态、突变和运动伪迹及跨设备阈值不可直接迁移。<https://pmc.ncbi.nlm.nih.gov/articles/PMC5413200/>
7. Kleckner IR et al. *Simple, Transparent, and Flexible Automated Quality Assessment Procedures for Ambulatory EDA Data* (2018). 经专家一致性验证的设备范围、±10 µS/s 斜率和伪迹过渡缓冲规则来源。<https://pmc.ncbi.nlm.nih.gov/articles/PMC5880745/>
8. Zhao Z, Zhang Y. *SQI Quality Evaluation Mechanism of Single-Lead ECG Signal Based on Simple Heuristic Fusion and Fuzzy Comprehensive Evaluation* (2018). 单导联 ECG 的 QRS 频带、峰度、基线及检测一致性多指标质量框架。<https://doi.org/10.3389/fphys.2018.00727>
9. Clifford GD et al. *Signal quality indices and data fusion for determining clinical acceptability of electrocardiograms* (2012). 在多位专家标注数据上验证多 SQI 融合，而非仅依赖 RR 规则。<https://doi.org/10.1088/0967-3334/33/9/1419>
10. Kelsey M et al. *Simulation of ambulatory electrodermal activity and the handling of low-quality segments* (2023). 比较低质量片段处理方式，并显示局部去除后插值可在一定伪迹密度下改善多类 EDA 特征恢复。<https://pubmed.ncbi.nlm.nih.gov/37863009/>
11. Makowski D et al. *NeuroKit2: A Python toolbox for neurophysiological signal processing* (2021). 当前 ECG/EDA 检测、清理、分解与 HRV 计算的主要软件实现依据。<https://doi.org/10.3758/s13428-020-01516-y>
