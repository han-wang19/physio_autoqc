# ECG CSV 列字典

适用于当前代码导出的 ECG 表。HRV 表位于 `ecg/results/`；`ecg/audit/` 保存过程记录。也可从 `results/reviewed_fif/` 自行提取特征。空值不等于零，通常表示不适用、未计算或处理失败，需结合状态和原因列判断。

## 时间、单位与状态

所有 `start_s`、`end_s`、`analysis_start_s`、`analysis_end_s` 都是相对于原始音乐记录起点的秒数，不是相对于裁剪后数组的时间。时间段按左闭右开 `[start_s, end_s)` 理解。RR 长度用毫秒；`*_fraction` 为 0–1 比例；布尔值用 True/False。

表的阶段不同，状态含义也不同。自动表的 REVIEW 不是“人工审核还没完成”。`qc_status` 描述纳入方式，是否最终成功要看 audit/ 最终处理报告的 `processing_status`。有效区间表和逐 RR 表只导出最终 INCLUDED 的记录。

## ecg_rr_intervals.csv：动态 HRV 的直接输入

每行是最终 R 峰中两个相邻峰之间的一个 RR，不是一个采样点或一个 QC 窗口。

| 列 | 含义与单位 |
| --- | --- |
| `record_id` | 记录编号，例如 A4_music5 表示被试 A4、曲目 5 |
| `rr_index` | 当前最终有效区间内的 RR 顺序，从 0 开始；不是原始采样点索引 |
| `start_s` | 该 RR 左侧 R 峰在原始记录中的时间，秒 |
| `end_s` | 该 RR 右侧 R 峰在原始记录中的时间，秒 |
| `rr_ms` | 两峰间距，毫秒；等于 `(end_s - start_s) * 1000`，允许浮点误差 |
| `rr_flagged` | 最终 RR 的固定算法残留异常提示；True 为触发，False 为未触发 |

### rr_flagged 到底表示什么

在最终 R 峰上，满足任一条件即标 True：

1. 由当前 RR 换算的心率低于 40 或高于 180 次/分，即 `60000 / rr_ms` 超出配置范围。
2. 当前 RR 与局部 RR 中位数的相对差超过 20%。局部中位数使用当前 RR、前两个和后两个 RR，边缘处使用实际存在的值，当前 RR 本身也参与中位数计算。

例如局部中位数是 800 ms，当前 RR 是 1100 ms，相对差是 37.5%，因此标 True，即使其换算心率仍在范围内。

这是 `ecg_qc._rr_artifact_summary` 的同一规则，最终处理后保存，不在每个动态窗口重新计算。阈值可能随实验配置变化，最终 NPZ 的 `rr_quality_hr_limits`、`rr_quality_local_deviation`、`rr_quality_rule` 保存实际使用值；以上是当前项目默认配置。

True **不表示人工确认错误、不表示已经修正，也不表示该 RR 已被排除**；False 也不是人工逐搏担保。它表达算法仍注意到的局部异常。即使记录整体 INCLUDED，也可能存在少量这种提示。RR 本身是观测到的峰间距，不能仅凭 HRV 列名里的 NN 就假定全部为无异常正常心搏。

动态提取只需汇总每窗已有标记。采用保留、修正或暂不计算含标记窗口的策略前，应固定分析方案；不能删行后让跨缺口的 RR 参与相邻差计算。

## ecg_valid_segments.csv：最终可用区间

每行是最终通过处理的一条记录的连续有效区间。目前每条 ECG 选择一个连续区间。

| 列 | 含义与单位 |
| --- | --- |
| `record_id` | 记录编号 |
| `start_s` | 最终信号起点在原始记录中的时间，秒 |
| `end_s` | 最终信号起点加实际数组时长，秒；末个 R 峰通常早于此时间 |
| `duration_s` | 最终数组长度除以采样率，秒；裁剪取整可能使它与原请求时长相差少量采样点 |
| `qc_status` | PASS/TRIMMED 为自动全段/裁剪纳入，MANUAL_PASS/MANUAL_TRIMMED 为人工全段/裁剪纳入 |
| `n_rr` | 最终 RR 数量，等于最终 R 峰数量减一 |
| `n_rr_flagged` | 其中 rr_flagged=True 的数量；没有删除这些行 |
| `frequency_domain_eligible` | 整段是否达到项目频域时长及有效修正比例要求；不是任意子窗口的资格 |

## final_ecg_hrv_features.csv：整段特征

每行使用整段最终 R 峰计算，不是动态滑窗结果。特征来自当前安装的 NeuroKit2；过短或不适用的指标可能为空。

| 列 | 含义与单位 |
| --- | --- |
| `record_id` | 记录编号 |
| `qc_status` | 该段纳入方式，同有效区间表 |
| `analysis_start_s` | 最终裁剪分析起点，秒 |
| `analysis_end_s` | 最终裁剪终点，秒；FIF 最后采样点及不包含右端点的裁剪可能产生毫秒级边界差；精确覆盖范围以有效区间表为准 |
| `usable_duration_s` | 实际处理数组时长，秒 |
| `standard_5min_hrv` | 实际时长是否至少 300 秒；不单独保证频域质量合格 |
| `mean_hr` | 代码中最终峰经 NeuroKit ecg_rate 计算后取平均，次/分；不保证与 `60000 / HRV_MeanNN` 完全相同 |
| `HRV_MeanNN` | RR 平均值，ms |
| `HRV_SDNN` | RR 样本标准差，ms |
| `HRV_RMSSD` | 相邻 RR 差的均方根，ms |
| `HRV_SDSD` | 相邻 RR 差的样本标准差，ms |
| `HRV_CVNN` | SDNN / MeanNN，无量纲 |
| `HRV_CVSD` | RMSSD / MeanNN，无量纲 |
| `HRV_MedianNN` | RR 中位数，ms |
| `HRV_MadNN` | RR 中位绝对偏差，按 NeuroKit mad 的实现缩放，ms |
| `HRV_MCVNN` | MadNN / MedianNN，无量纲 |
| `HRV_IQRNN` | RR 四分位距，ms |
| `HRV_SDRMSSD` | SDNN / RMSSD，无量纲 |
| `HRV_Prc20NN` | RR 第 20 百分位，ms |
| `HRV_Prc80NN` | RR 第 80 百分位，ms |
| `HRV_pNN20` | 相邻 RR 绝对差 >20 ms 的百分数；当前 NeuroKit 实现用 RR 数量作分母 |
| `HRV_pNN50` | 相邻 RR 绝对差 >50 ms 的百分数；同上 |
| `HRV_MinNN` | 最小 RR，ms |
| `HRV_MaxNN` | 最大 RR，ms |
| `HRV_SDANN1`、`HRV_SDANN2`、`HRV_SDANN5` | 分别以 1/2/5 分钟子段平均 RR 计算其标准差，ms；NeuroKit 要求至少三个子段，不足则为空 |
| `HRV_SDNNI1`、`HRV_SDNNI2`、`HRV_SDNNI5` | 分别以 1/2/5 分钟子段计算段内 RR 标准差后取平均，ms；同样有最少子段要求 |
| `HRV_HTI` | RR 个数 / RR 直方图最高柱计数，无量纲，受分箱设置影响 |
| `HRV_TINN` | RR 直方图三角插值的基线宽度，ms |
| `HRV_ULF` | 0–0.0033 Hz 频带谱积分；本项目短记录不适合解释 ULF，可能为空 |
| `HRV_VLF` | 0.0033–0.04 Hz 频带谱积分；短记录解释受限 |
| `HRV_LF` | 0.04–0.15 Hz 频带谱积分 |
| `HRV_HF` | 0.15–0.4 Hz 频带谱积分 |
| `HRV_VHF` | 0.4–0.5 Hz 频带谱积分 |
| `HRV_TP` | NeuroKit 计算的各频带总谱积分 |
| `HRV_LFHF` | LF/HF，无量纲；不要直接视为交感/副交感活动的独立测量 |
| `HRV_LFn` | 当前安装的 NeuroKit 实现为 LF/TP，无量纲；不是直接假定 LF/(LF+HF) |
| `HRV_HFn` | 当前安装的 NeuroKit 实现为 HF/TP，无量纲 |
| `HRV_LnHF` | HF 的自然对数 |

频域函数当前沿用 NeuroKit `normalize=True` 默认设置，因此这些谱积分不能直接标为未归一化的绝对功率 ms²。若研究需要绝对功率，应在特征提取方案中明确统一 normalize 参数并重新生成相应结果，不能只改列名或单位。

## audit/：过程表是什么

| 表 | 每行代表什么 |
| --- | --- |
| `ecg_preprocessing_report.csv` | 一条记录的预处理结果 |
| `ecg_window_qc.csv` | 一个原始自动滑窗，通常为 30 秒窗、5 秒步长；末窗对齐终点 |
| `ecg_rr_bad_runs.csv` | 检测器原始峰所识别的一段连续异常 RR，并非最终 RR 表中的固定标记 |
| `ecg_qc_report.csv` | 一条记录的原始自动决定 |
| `ecg_qc_effective_summary.csv` | 一条记录的自动结果叠加人工决定；保留大量自动指标作对照，不是全部指标重算 |
| `ecg_qc_effective_intervals.csv` | 人工覆盖后具有相同状态的连续时间段，不再是固定长度滑窗 |
| `final_ecg_processing_report.csv` | 一条记录最终处理成功、未纳入或失败的状态及诊断指标 |

以下列按名称查阅；同名时间、record_id 与前文一致。异常统计必须结合所在表判断来源，自动阶段统计不能替代最终残留统计。

| 列 | 含义与单位 |
| --- | --- |
| `status` | 预处理表为 OK/ERROR 等处理状态；窗口表及有效区间表为 PASS/REVIEW/FAIL |
| `reason` | 当前表所属阶段的原因代码或错误信息；有效汇总中保留的 reason 主要仍来自自动阶段 |
| `source_file` | 原始 FIF 路径 |
| `sfreq_hz` | ECG 采样率，Hz |
| `initial_peak_count` | 初始检测峰数 |
| `fixed_peak_count` | RR 峰修正后的峰数 |
| `final_peak_count` | 初始预处理最后对齐后的峰数；此处 final 不代表最终审核阶段 |
| `ecg_sd` | 当前记录/窗口滤波 ECG 的标准差，沿用数组幅度单位，不自动标为 µV |
| `ecg_range` | 初始预处理 ECG 的幅度极差，沿用数组单位 |
| `rpeaks` | 自动窗口内的峰数量 |
| `minimum_expected_rpeaks` | 根据窗口时长及配置最低心率推算的最低峰数 |
| `template_correlation` | 当前窗口/最终区间心搏与中位 QRS 模板相关系数的中位数，无量纲 |
| `rr_plausible_fraction` | 心率落在配置合理范围内的 RR 比例 |
| `rr_outlier_fraction` | 超过局部 RR 相对偏离阈值的比例 |
| `signal_quality` | Zhao–Zhang SQI 类别：Excellent/Barely acceptable/Unacceptable |
| `short_rr_correction_fraction` | 自动窗口中与短异常段相关的异常 RR 数 / 窗口 RR 数；不是人工修正比例 |
| `uncorrectable_rr_run` | 自动窗口是否与不满足短段修正规则的异常段重叠 |
| `overlapping_bad_run_count` | 与自动窗口重叠的异常连续段数量 |
| `run_number` | 自动异常连续段在记录内的编号 |
| `abnormal_rr_count` | 该连续异常段中异常 RR 的数量 |
| `correctable_short_run` | 是否满足最多连续 3 个且总长不超过 5 秒的短段规则；不是已修正的证明 |
| `severe_long_run` | 该异常段是否超过配置的严重长段时长，当前 >30 秒 |
| `auto_status` | 原始记录级自动 PASS/TRIMMED/REVIEW/FAIL |
| `automatic_status` | 有效汇总中保存的原始 auto_status 副本 |
| `effective_status` | 人工决定覆盖后的 PASS/TRIMMED/REVIEW/FAIL，不单独代表最终处理成功 |
| `manual_status` | INCLUDE_FULL/INCLUDE_TRIMMED/EXCLUDE/SECOND_REVIEW；空白表示尚无人工决定 |
| `pass_window_fraction` | 原始自动 PASS 窗口数 / 全部窗口数；不是人工后有效时间比例 |
| `rpeak_count` | 所属自动记录/最终区间的峰数量，须结合表的阶段读取 |
| `excellent_sqi_window_fraction` | 自动 Excellent SQI 窗口比例 |
| `borderline_sqi_window_count` | 自动 Barely acceptable SQI 窗口数量 |
| `unacceptable_sqi_window_count` | 自动 Unacceptable SQI 窗口数量 |
| `detected_rpeak_count` | 所属阶段未经峰修正的检测峰数 |
| `corrected_rpeak_count` | 所属阶段用于后续计算的峰数；最终处理时可能包含人工峰 |
| `short_bad_rr_count` | 所属阶段检测器原始峰中属于可修正短异常段的 RR 数 |
| `short_bad_rr_fraction` | 上述短异常 RR 数 / 原始检测 RR 总数 |
| `longest_bad_run_s` | 检测器原始峰中最长连续异常段时长，秒 |
| `total_bad_duration_s` | 检测器原始峰中异常连续段总时长，秒 |
| `time_domain_eligible` | 时长是否达到当前记录配置的项目最低值，默认 120 秒；条件门槛可在私人配置中调整，不是所有 HRV 指标有效性的保证 |
| `standard_5min_hrv_eligible` | 时长是否至少 300 秒 |
| `manual_pass_interval_count`、`manual_review_interval_count`、`manual_fail_interval_count` | 对应标签的已保存人工标注条数；可能有历史覆盖重叠，不是有效连续段数 |
| `manual_pass_duration_s`、`manual_review_duration_s`、`manual_fail_duration_s` | 对应已保存人工标注长度之和，秒；覆盖重叠可能重复计时，不等于最终可用时长 |
| `processing_status` | INCLUDED=最终处理成功；NOT_INCLUDED=QC 决定不纳入；EXCLUDED_AFTER_FINAL_CHECK=尝试处理但最终检查失败 |
| `manual_reviewed_at` | 本次最终处理所对应的人工审核版本时间，用于防止旧结果覆盖新决定 |
| `used_manual_rpeaks` | 是否实际使用人工编辑后保存的峰；False 不代表未人工审核记录 |
| `manual_rpeak_count` | 实际采用的人工峰在最终区间内的峰数；未采用则为 0 |
| `manual_rpeak_edit_count` | 完整记录的人工峰相对自动峰的编辑数量估计，取删除数与新增数的较大值；并非单个区间内精确编辑次数 |
| `manual_rpeak_edit_fraction` | 上述完整记录编辑数 / 最终区间自动峰数（分母至少为 1），不是逐 RR 异常比例，可能超过 1 |
| `effective_rpeak_correction_fraction` | 有人工峰时采用上述人工编辑比例；否则采用最终区间原始检测短异常 RR 比例，供频域资格判断 |
| `residual_bad_rr_count` | 最终实际采用 R 峰上仍触发异常规则的 RR 数，与 results/ 逐 RR 表的标记数一致 |
| `residual_bad_rr_fraction` | 上述数量 / 最终 RR 数 |
| `residual_longest_bad_run_s` | 最终 R 峰上残留异常连续段的最长时长，秒 |
| `peak_correction_fraction` | 检测器峰修正前后峰数差的绝对值 / 初始峰数；只反映数量变化，不表示全部位置变化 |

`frequency_domain_eligible`、`analysis_start_s`、`analysis_end_s`、`usable_duration_s` 与前文单位一致，但自动报告对应自动区间，最终报告对应最终区间。有效汇总虽更新人工区间和部分时长字段，SQI、异常统计及频域标志仍可能来自自动阶段；使用现成 HRV 时读取 results/ 中三张最终表。
