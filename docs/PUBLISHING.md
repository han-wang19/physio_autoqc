# 只发布流程，不发布数据

## 发布内容

提交 Python/网页代码、requirements、通用文档、project_config.example.json、许可证和手写格式示例。examples 中的事件表是虚构示例，不对应真实记录。合成演示数据由脚本现场生成，FIF 本身不进入仓库。

私人 project_config.local.json、真实输入、事件表、人工审核结果、输出图表、NPZ/FIF 和 .local_docs 笔记都被忽略。个人配置可以放在仓库外，通过 PHYSIO_CONFIG 指定。

`.gitignore` 不会自动停止跟踪已经提交过的文件。发布前必须核对当前 Git 列表及历史。工具不自动删除数据、修改历史或推送。

## 检查与提交

```powershell
python -m unittest discover -s tests -q
python scripts/check_public_files.py
git status --short
git diff --stat
```

审查通过后按明确列表暂存代码及文档，避免将私人目录强制加入：

```powershell
git add -u
git add .gitignore README.md requirements.txt project_config.py config_profiles.py project_config.example.json main.py review_app.py processing review reports scripts tests docs examples
# 若已选定许可证，同时添加 LICENSE、THIRD_PARTY_NOTICES.md
git diff --cached --name-only
git diff --cached --stat
python scripts/check_public_files.py
git commit -m "Prepare reusable physiology QC workflow"
```

本地已经存在 origin 时无需重新添加。没有远程仓库时，先在 GitHub 创建一个空仓库，再设置自己的 URL：

```powershell
git remote -v
git remote add origin https://github.com/YOUR-USER/physio_autoqc.git
git push -u origin main
```

`remote add` 只在没有 origin 时执行。确认推送目标、分支和权限后才推送；代码审查及配置整理本身不等于完成发布。

## 检查历史

```powershell
git rev-list --objects --all
```

当前候选文件检查不代表历史也已清理。若旧提交含私人配置路径，更新工作区不能从历史中抹去它们；若希望完全不公开这些历史信息，可将审查过的纯代码快照放进一个全新的仓库。不要直接删除或重写仍需保留的原仓库历史。

## 使用者的流程

下载代码 → 安装依赖 → 复制示例 JSON 为私人配置 → 设置自己的数据、单位与端口 → 预处理 → QC → 本机人工审核 → 生成最终特征。运行不需要维护者的数据或电脑。

许可证确定软件允许如何使用，第三方代码保留其原许可。文件检查不能替代人工审查，也不会判断科研阈值是否适合新的设备和人群。
