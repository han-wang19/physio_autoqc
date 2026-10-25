# ECG 单位与设备兼容

MNE 标准电压数据使用 V。处理前确认实际样本与元数据是否一致，不能仅根据幅值推断单位。

`processing/units.py` 的判断顺序：

1. 已有 voltage_calibrated=true、unit=V 标记：不重复换算。
2. 私人配置明确指定 ECG_INPUT_UNIT 为 V/mV/uV：将尚未校准的输入数值换算为 V。
3. 找到匹配的原始 BDF：核对文件头。标准 uV/mV/V 已由 MNE 导入转换为 V，不重复转换。
4. 无原始 BDF 时，按原生 FIF 的 V 元数据处理。

历史 FIF 数值并非 V 却被写成 V 时，必须配置实际数值单位或提供原始依据，不能把错误标签视为校准结果。

## Neuracle 的 uV/mm 兼容解释

部分 Neuracle BDF 在电压物理单位字段使用 uV/mm。MNE 可能不识别并保留标定后数值。仅在配置设备为 neuracle 且原始头匹配时，项目将该数值按 µV 解释，乘 1e-6 转为 V。

这是基于厂商输入量程、物理/数字范围标定及官方读取实现的交叉解释，不是厂商对该字符串的逐字保证，也不是所有设备的通用规则。不要额外乘任意纸面毫米或显示增益。

- [厂商 FAQ](https://www.neuracle.cn/newsinfo/311703.html)
- [厂商 BDF 读取实现](https://github.com/neuracle/neuracle-eegfile-reader/blob/master/read_bdf.m)

NPZ、报告和审核 FIF 会记录来源及换算系数。旧缓存没有校准标记时，导出仍按来源解析进行兼容换算。

显示缩放与单位校准不同。`raw.plot(scalings={'ecg':'auto'}, clipping=None)` 仅调整显示，不修改数值单位。
