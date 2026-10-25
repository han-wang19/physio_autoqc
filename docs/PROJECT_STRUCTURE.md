# 项目结构与维护

## 入口与职责

| 路径 | 职责 |
| --- | --- |
| `app.py` | 日常统一入口：选择实验、预处理与质检、打开人工复核 |
| `main.py`、`project_config.py` | 批处理入口与实验配置 |
| `scripts/start_review.py` | 按实验名和固定端口启动人工核查网页 |
| `scripts/extract_physio_channels.py` | 从切分 FIF 或 BDF/events.csv 提取双通道生理数据 |
| `processing/ecg_*.py`、`processing/eda_*.py` | 各信号的预处理、质检、最终处理与分析导出 |
| `processing/common.py` | 文件路径、输出布局及公共辅助函数 |
| `review_app.py` | Streamlit 启动、组件注册、页面级错误提示 |
| `review/workflow.py` | 审核队列、筛选、草稿、撤销与保存流程 |
| `review/data.py` | 报告与 NPZ 读取、窗口状态整理；不依赖 Streamlit |
| `review/figures.py` | ECG/EDA 图表、状态底色与人工边界；EDA 依赖按需加载 |
| `review/manual_qc.py`、`review/outputs.py` | 人工决定持久化、有效结果及输出同步 |
| `review/chart/index.html` | 图表组件页面骨架 |
| `review/chart/chart.css` | 图表浮层、菜单与焦点样式 |
| `review/chart/chart.js` | Plotly 渲染、选区、R 峰编辑与窗口适配 |
| `review/chart/plotly.min.js` | 本地 Plotly 依赖，支持离线使用 |
| `reports/generate_qc_summary.py`、`reports/qc_coverage.py` | 汇总与覆盖图 |
| `tests/test_review_regressions.py` | 回归验证，使用临时目录隔离研究输出 |
| `docs/`、`docs/QC_METHODS.md` | 结构说明、字段字典与方法说明 |
| `outputs/<实验>/` | 中间数据、人工审核记录和分析产物 |

批处理与网页启动入口保持原路径。处理算法放在 `processing/`，审核代码和图表资源集中在 `review/`，汇总图放在 `reports/`，测试与说明各自归档。草稿与撤销留在 `review/workflow.py`。导入 `review_app` 不启动页面，测试直接导入实际函数，不再解析入口源码。前端按结构、样式、行为分文件维护；调整图表时不要直接编辑压缩的 Plotly 库。

## 想修改什么，就看哪里

- 改数据目录、实验名称或阈值：`project_config.py`。
- 改处理步骤：`processing/` 下对应信号的 `*_preprocess.py`、`*_qc.py`、`*_finalize.py`。
- 改网页标题、启动和汇总按钮：`review_app.py`。
- 改审核队列、记录搜索、撤销和保存：`review/workflow.py`。
- 改报告读取或窗口状态整理：`review/data.py`。
- 改曲线、坐标轴、状态底色：`review/figures.py`。
- 改鼠标交互：`review/chart/chart.js`；改菜单外观：`review/chart/chart.css`。
- 查最终数据：`outputs/<实验>/`；查列的含义：`docs/ECG_CSV_COLUMNS.md`。

## 运行与验证

```powershell
python main.py --stage qc
python scripts/start_review.py --experiment DEMO
python -m reports.generate_qc_summary
python -m unittest discover -s tests -q
node --check review/chart/chart.js
```

网页支持按记录编号搜索及查看当前信号的审核进度。搜索只筛选队列，不修改审核结果。修改仍需点击“完成审阅”保存；之后点击“结束本次审核并更新汇总”同步分析产物。

`outputs/` 与 `__pycache__/` 已由 Git 忽略。研究输出不要为了整理目录而手动迁移或删除；路径管理集中在 `processing/common.py`。

ECG 与 EDA 输出均分为 `results/`（分析结果）、`internal/`（NPZ 等内部依赖）、`audit/`（过程报告）和 `plots/`（图片），由 `processing/common.py` 统一管理路径。人工核查地址和各实验的独立启动命令见 [README](../README.md#manual-review)。
