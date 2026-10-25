# 实验配置与通道提取

公共 `project_config.py` 保存共用处理参数，`config_profiles.py` 负责读取和校验实验配置。真实实验放入被忽略的 `project_config.local.json`。公共模板 `project_config.example.json` 仅含 DEMO 与相对路径。

## 私人配置

| 顶层字段 | 含义 |
| --- | --- |
| default_experiment | 默认实验，必须在 experiments 中登记 |
| output_root | 结果目录，建议 outputs 或项目外部的私人目录 |
| experiments | 以实验名为索引的配置对象，可登记多个实验 |

相对路径相对于 JSON 所在目录。`PHYSIO_CONFIG` 可指定其他 JSON；指定文件不存在会报错。`PHYSIO_EXPERIMENT` 或批处理的 `--experiment` 选择当前实验。

| 实验字段 | 含义 |
| --- | --- |
| DATA_DIR、FILE_PATTERN | 来源生理 FIF 的目录和匹配规则 |
| ECG_CHANNEL、EDA_CHANNEL | FIF 内通道名 |
| ECG_INPUT_UNIT | auto 使用校准标记、原始 BDF 或原生 FIF 单位；V/mV/uV 明确指定尚未校准的输入数值单位 |
| RAW_BDF_DIR、RAW_BDF_PATTERN | 可选原始 BDF 依据；模式支持 `{subject}` 或 `{subject:02d}` |
| ACQUISITION_SYSTEM | 默认 generic；确认是 Neuracle 时可填 neuracle，启用其特定文件头兼容规则 |
| EDA_INPUT_UNIT | S、mS、uS/µS/μS、nS，根据来源核实 |
| EDA_DEVICE_RANGE_US | 已换算成 µS 的设备有效范围，不能从样本分布猜测 |
| EDA_CALIBRATION_SOURCE | 输入单位和量程依据 |
| MIN_VALID_SECONDS_BY_CONDITION | 例如 `{"baseline":{"ecg":60,"eda":60}}`；空对象使用默认 120 秒 |
| REVIEW_PORT | 本机审核端口；并行实验使用不同端口 |

原生 V FIF 不要求补造 BDF。历史 FIF 数值与元数据不符时，必须提供实际单位或原始依据；已校准数据不重复缩放。各实验的条件门槛独立复制，未知字段、错误单位和非法范围会明确报错。

## 从已切分 FIF 提取

可选提取脚本要求来源通道名为 GSR、ECG，保留原数值、采样率和事件，不进行滤波或峰检测：

```powershell
python scripts/extract_physio_channels.py --format split-fif --source data/source_segments --destination data/MY_STUDY/physiology --subjects 1 2 3
```

输入如 Exp1_resting_raw.fif、Exp1_pass1_raw.fif、Exp1_pass2_raw.fif；输出为 A1_resting_physio_raw.fif、A1_active_physio_raw.fif、A1_passive_physio_raw.fif。只选择请求的被试，不会将 Exp10 误认为 Exp1。

## 从 BDF 和事件表提取

TIME/DURATION 为相对于来源 BDF 起点的秒数，LABEL 使用 rest1/resting、pass1、pass2。第二段静息使用 rest2，并开启 `--two-resting`。[examples](../examples/README.md) 仅为虚构格式示例。

```powershell
python scripts/extract_physio_channels.py --format bdf-events --source data/source_bdf --events-dir data/source_events --destination data/MY_STUDY/physiology --subjects 1 2 3 --bdf-pattern 'sub-A{subject:02d}_eeg.bdf' --events-pattern 'A{subject:02d}_events.tsv' --two-resting --dry-run
```

确认映射后去掉 `--dry-run`。CSV 用逗号分隔，TSV 用制表符；省略 events-dir 时使用来源目录。

- pass1→active、pass2→passive、rest1/resting→resting、rest2→resting2。
- 切分为最近采样点的左闭右开区间，清单保留请求时间与实际采样点边界。
- 先读取完整的两个通道再切分；MNE 会将混合采样率通道对齐到所选通道的最高采样率。EDA 预处理随后按配置重采样。
- 目标文件存在时拒绝覆盖；缺少条件、缺少通道、越界或重复条件会报错。
- 保留事件表中的条件重叠，不自动移动边界；清单的 overlap_with_other_conditions_s 记录重叠时长。
- FIF、extraction_manifest.csv 和来源信息属于私人产物，不提交 GitHub。

输出会回读校验波形、采样率、长度和有效 annotation。提取完成后仍需运行 `main.py --experiment MY_STUDY --stage preprocess` 和 `--stage qc`。
