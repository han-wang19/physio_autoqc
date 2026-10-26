# Physio AutoQC

一个在自己电脑上运行的 ECG / EDA 预处理与人工复核工具。自动检测后，在网页中检查波形、修改 R 峰和坏段，保存带标记的 FIF，也可导出特征表。

## 第一次使用

需要 Python 3.10+。在项目目录中安装依赖：

```powershell
python -m pip install -r requirements.txt
```

复制配置模板，填写自己的实验名称和数据文件夹：

```powershell
Copy-Item project_config.example.json project_config.local.json
```

主要修改 `default_experiment`、`experiments` 中的实验名和 `DATA_DIR`；核对 ECG / EDA 通道名、单位及 EDA 量程。其余选项可先保持默认，字段解释见[配置说明](docs/EXPERIMENTS.md)。macOS / Linux 复制文件使用 `cp`。

每个输入 FIF 包含 ECG 和 GSR/EDA 通道，文件名例如 `A1_active_physio_raw.fif`、`A2_resting_physio_raw.fif`。ECG 原生 FIF 使用 V 单位；其他输入单位须正确配置。原始 BDF 可按[提取说明](docs/EXPERIMENTS.md)先提取和切分。

## 每次使用

```powershell
python app.py
```

按提示选择实验，再选择：

- **预处理 + 自动质检**：首次处理数据时使用，完成后启动人工复核网页。
- **直接打开人工复核**：已有处理结果时使用。

<a id="manual-review"></a>

打开终端显示的网址。在网页左侧选择 ECG / EDA 和记录，检查波形、编辑坏段和 R 峰后点击 **完成审阅**。记录按 A1、A2、…、A10 排序，无需补零重命名。

审核结束时，点击 **结束本次审核并更新汇总**，更新最终特征和报告。终端需要保持运行，按 `Ctrl+C` 停止服务。网页默认仅在本机可访问，GitHub 提供的是程序代码。

更新时会显示当前记录和处理阶段，只处理新增审核、修改或输入发生变化的记录。首次计算最终特征需要一些时间；再次点击会跳过已完成的部分，中断后也可继续。

## 结果在哪里

`outputs/<实验名>/ecg/results/` 和 `outputs/<实验名>/eda/results/`。

其中 `reviewed_fif/` 保存审核后的单通道波形：R 峰用零时长 `R_PEAK` 标记，坏段用区间 `BAD_segment` 标记。最终是否纳入分析，以对应 `audit/` 下的处理报告为准。其他目录由程序使用，请保留。

自己的数据、本地配置和输出默认不会提交到 GitHub。仓库不附带真实研究数据。

## 先试一下

新下载、尚未创建本地配置时，可以用合成数据体验：

```powershell
python scripts/make_demo_data.py
python app.py --experiment DEMO --action prepare
```

演示数据只用于了解操作，也可能被质量规则判定为不纳入。

更多说明：[详细操作与命令](docs/ADVANCED_USAGE.md) · [QC 方法](docs/QC_METHODS.md) · [FIF 读取](docs/ECG_FIF_USAGE.md) · [发布到 GitHub](docs/PUBLISHING.md) · [验证范围](docs/VALIDATION.md)

许可证：[MIT](LICENSE)。第三方声明见 [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES.md)。
