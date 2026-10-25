# Physio AutoQC

用于 ECG / EDA 的本地科研流程：预处理、窗口质检、人工区间与 R 峰复核、审核后 FIF 导出，以及最终特征计算。使用 Python 3.10+。

仓库只提供代码、方法说明和虚构格式示例，**不包含参与者数据、真实事件表、人工审核记录或研究结果**。运算和网页默认在使用者自己的电脑上运行。默认阈值需要结合自己的设备与实验验证。

## 1. 安装

在项目根目录使用独立 Python 环境：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

macOS / Linux 激活环境使用 `source .venv/bin/activate`。

## 2. 配置自己的实验

复制公共模板为私人配置文件：

```powershell
Copy-Item project_config.example.json project_config.local.json
```

macOS / Linux 使用 `cp project_config.example.json project_config.local.json`。修改其中的实验名称、`DATA_DIR`、`FILE_PATTERN`、通道名、输入单位、设备量程和 `REVIEW_PORT`。不要把真实路径写回公共模板。

相对路径相对于 JSON 所在目录。`PHYSIO_CONFIG` 可指定另一份私人 JSON，`PHYSIO_EXPERIMENT` 可选择默认实验。没有私人配置时，程序提供一个不含数据的 DEMO，预期来源目录为 `data/DEMO/physiology/`。

私人配置和 `data/`、`outputs/` 均被 Git 忽略，不随代码提交。数据也可放在项目外部。完整字段见 [实验配置](docs/EXPERIMENTS.md)。

### 输入文件约定

每条源 FIF 应包含所配置的 ECG、GSR/EDA 通道，命名为：

```text
A1_baseline_physio_raw.fif
A1_task_physio_raw.fif
A2_task_physio_raw.fif
```

`A<被试编号>_<条件>_physio_raw.fif` 中条件可自行命名，`music1` 等也支持。ECG 原生 FIF 应使用 MNE 的 V 单位；其他数值单位需在私人配置中明确指定。EDA 单位和设备范围需单独核对，不能根据幅值猜单位。

完整 BDF 与 CSV/TSV 时间表可先使用 [通道提取脚本](docs/EXPERIMENTS.md#从-bdf-和事件表提取)；该可选脚本目前要求来源通道名为 ECG/GSR。[examples/](examples/README.md) 只包含手工编写的格式示例。

### 使用完全合成的数据试运行

不创建私人配置时：

```powershell
python scripts/make_demo_data.py
python main.py --experiment DEMO --stage preprocess
python main.py --experiment DEMO --stage qc
python scripts/start_review.py --experiment DEMO
```

生成的人工 FIF 位于被忽略的 `data/`；它只演示软件流程，不是算法有效性验证。已有私人配置时，可通过 `PHYSIO_CONFIG` 指向公共模板，并使用单独目录运行演示，避免混入自己的实验。

## 3. 预处理与自动质检

将 DEMO 换成 JSON 中登记的实验名称：

```powershell
python main.py --experiment DEMO --stage preprocess
python main.py --experiment DEMO --stage qc
```

ECG 使用滤波、峰检测和统一的形态对齐；EDA 进行单位转换、清理、tonic/phasic 分解和 SCR 检测。自动 QC 生成 PASS、TRIMMED、REVIEW、FAIL。预处理失败会记录到报告并终止后续流程，修复原因后再重跑。

默认最低有效时长为 120 秒。短基线可在私人配置中设置条件门槛，但不自动赋予标准 5 分钟或频域 HRV 资格，见 [QC 方法](docs/QC_METHODS.md)。

<a id="manual-review"></a>

## 4. 人工核查入口

```powershell
python scripts/start_review.py --experiment DEMO
```

终端会输出实验地址。示例配置为 **[http://127.0.0.1:8501](http://127.0.0.1:8501)**；保持终端运行，在同一页面左侧切换 ECG/EDA。页面顶部显示实验名称。可用 `--port 8510` 临时更换端口，按 `Ctrl+C` 停止服务。

地址仅在服务启动时可用。`127.0.0.1` 只供本机访问；把代码放到 GitHub 不会自动提供远程审核服务。多个实验同时运行时应设置不同端口。

审核步骤：选择记录 → 填写审阅者 → 标注 PASS/FAIL、编辑 R 峰 → 点击“完成审阅”。自动 PASS 也需要人工确认；人工完成状态与数据质量状态分开。

- **完成审阅**：保存区间、峰和记录决定，并导出审核终版 FIF。
- **结束本次审核并更新汇总**：同步最终特征、图片和有效汇总。人工暂定纳入仍可能被最终检查排除，原因见处理报告。

## 5. 后续分析入口

```text
outputs/<实验名称>/
├─ config_snapshot.json
├─ manual_qc_decisions.csv、manual_qc_intervals.csv
├─ ecg/
│  ├─ results/     reviewed_fif/、HRV 特征和 RR 表
│  ├─ internal/    预处理、人工峰、最终 NPZ
│  ├─ audit/       自动/有效 QC、最终处理报告
│  └─ plots/       各阶段图片
└─ eda/           同样分为 results、internal、audit、plots
```

日常从 `results/` 开始。ECG 审核 FIF 只有一个 ECG 通道，R 峰为零时长 `R_PEAK` annotation；ECG/EDA 的不可用区间均合并为不重叠的 `BAD_segment`。保留区间存于元数据，不拼接坏段两侧的信号或 RR。

审核 FIF 不表示已通过最终特征检查，现有流程是否纳入以 `audit/final_*_processing_report.csv` 为准。内部 NPZ 仍被程序依赖，不要随意删除。

- [ECG FIF 读取示例](docs/ECG_FIF_USAGE.md)
- [EDA FIF 读取示例](docs/EDA_FIF_USAGE.md)
- [ECG CSV 列字典](docs/ECG_CSV_COLUMNS.md)
- [单位与设备兼容](docs/ECG_UNIT_CHECK.md)

也可运行 `python main.py --experiment DEMO --stage finalize`。默认 `--stage all` 连续预处理、QC 和最终计算，不等待人工审核，适合自动初筛。

## 6. 开发与发布

```powershell
python -m unittest discover -s tests -q
node --check review/chart/chart.js
python scripts/check_public_files.py
```

Node 只用于可选的 JavaScript 语法检查，不是网页运行依赖。项目结构见 [维护说明](docs/PROJECT_STRUCTURE.md)，发布操作见 [GitHub 发布流程](docs/PUBLISHING.md)。科研方法与局限见 QC 方法文档。

许可为 [MIT](LICENSE)，第三方库声明见 [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES.md)。测试环境及验证边界见 [验证说明](docs/VALIDATION.md)。
