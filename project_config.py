"""公共处理参数；私人实验配置从 project_config.local.json 加载。"""

from pathlib import Path
import os
from copy import deepcopy
from config_profiles import ExperimentConfig, load_profiles

PROJECT_ROOT = Path(__file__).resolve().parent
EXPERIMENTS, DEFAULT_EXPERIMENT, OUTPUT_ROOT, CONFIG_SOURCE = load_profiles(PROJECT_ROOT)
EXPERIMENT_NAME = os.environ.get('PHYSIO_EXPERIMENT', DEFAULT_EXPERIMENT)
os.environ.setdefault('NUMBA_CACHE_DIR', str(PROJECT_ROOT / '.cache' / 'numba'))

DATA_DIR: Path
FILE_PATTERN: str
RAW_BDF_DIR: Path
RAW_BDF_PATTERN: str
ACQUISITION_SYSTEM: str
ECG_CHANNEL: str
ECG_INPUT_UNIT: str
EDA_CHANNEL: str
EDA_INPUT_UNIT: str
EDA_DEVICE_RANGE_US: tuple[float, float]
MIN_VALID_SECONDS_BY_CONDITION: dict[str, dict[str, float]] = {}
EDA_CALIBRATION_SOURCE: str
EDA_SCALE_TO_US: float
REVIEW_PORT: int

# 要处理的信号
RUN_ECG = True
RUN_EDA = True

# 绘图
SAVE_PLOTS = True
SHOW_PLOTS = False
PLOT_DPI = 180

# ECG 预处理
ECG_BANDPASS_HZ = (0.5, 40.0)
ECG_NOTCH_HZ = 50.0
ECG_PRIMARY_DETECTOR = "promac"
# Candidate detectors can land on either side of a biphasic QRS. This window
# is used to infer one record-level dominant polarity and align all beats to it.
ECG_LOCAL_REFINEMENT_MS = 50.0
# Shared morphology rule for a substantial R lobe followed by a deeper S lobe.
ECG_RS_MIN_POSITIVE_RATIO = 0.5
ECG_RS_MIN_CONSISTENCY = 0.9
ECG_RS_SEPARATION_MS = (8., 50.)

# 通用 QC 窗口
QC_WINDOW_SECONDS = 30.0
QC_STEP_SECONDS = 5.0

# ECG 自动 QC
ECG_MIN_SIGNAL_SD = 1e-12
ECG_PLAUSIBLE_HR_BPM = (40.0, 180.0)
ECG_PASS_RR_FRACTION = 0.98
ECG_FAIL_RR_FRACTION = 0.90
ECG_PASS_RR_OUTLIER_FRACTION = 0.05
ECG_FAIL_RR_OUTLIER_FRACTION = 0.20
ECG_RR_LOCAL_DEVIATION = 0.20
ECG_MAX_SHORT_RR_CORRECTION_FRACTION = 0.05
# 频域 HRV 对伪迹更敏感；超过 2% 仅保留时域结果。
ECG_MAX_FREQUENCY_RR_CORRECTION_FRACTION = 0.02
ECG_MAX_SHORT_BAD_RUN_SECONDS = 5.0
ECG_MAX_CONSECUTIVE_CORRECTABLE_RR = 3
ECG_LONG_BAD_RUN_SECONDS = 30.0
ECG_MIN_RPEAKS = 30
# 只有全部窗口PASS才把完整片段标为PASS；否则尝试选择连续干净区间。
ECG_MIN_PASS_WINDOW_FRACTION = 1.0
ECG_MIN_TIME_DOMAIN_SECONDS = 120.0
# 标准短时 HRV 至少 5 分钟；条件特定的纳入门槛另见私人实验配置。
ECG_STANDARD_HRV_SECONDS = 300.0
ECG_MIN_FREQUENCY_DOMAIN_SECONDS = 300.0

# EDA 预处理与自动 QC
# 输入单位与校准来源随实验选择，未知单位会在选择实验时失败。
EDA_UNIT_TO_US = {
    "S": 1_000_000.0,
    "mS": 1_000.0,
    "uS": 1.0,
    "µS": 1.0,
    "μS": 1.0,
    "nS": 0.001,
}
EDA_OUTPUT_UNIT = "uS"
EDA_TARGET_SFREQ = 100.0
EDA_EDGE_TRIM_SECONDS = 5.0
EDA_SCR_DETECTOR_MIN = 0.01
EDA_SCR_ABSOLUTE_MIN_US = 0.05
EDA_QC_WINDOW_SECONDS = 5.0
EDA_QC_STEP_SECONDS = 5.0
EDA_ARTIFACT_BUFFER_SECONDS = 5.0
EDA_MIN_FINITE_FRACTION = 0.999
EDA_LOW_CONDUCTANCE_REVIEW_US = 0.05
EDA_REVIEW_LOW_CONDUCTANCE_FRACTION = 0.95
EDA_MAX_ABS_SLOPE_US_PER_S = 10.0
# Boucsein et al. (2012, Psychophysiology, doi:10.1111/j.1469-8986.2012.01384.x)
# 建议可计数 SCR 的最低幅度为 0.01–0.05 µS，并强调结合设备噪声判断。
# 这里取保守下限 0.01 µS 作为 5 秒窗口“缺少可辨识变化”的人工复核触发值，
# 不是把无 SCR 等同于伪迹；平线可能是真实低反应，因此只标 REVIEW，不直接 FAIL。
EDA_MIN_WINDOW_RANGE_REVIEW_US = 0.01
EDA_REVIEW_FLATLINE_FRACTION = 0.95
EDA_MIN_VALID_SECONDS = 120.0
EDA_MIN_PASS_WINDOW_FRACTION = 1.0


def select_experiment(name):
    """Select a registered dataset for both batch processing and the review UI."""
    global EXPERIMENT_NAME, DATA_DIR, FILE_PATTERN, RAW_BDF_DIR, RAW_BDF_PATTERN
    global ACQUISITION_SYSTEM, ECG_CHANNEL, ECG_INPUT_UNIT, EDA_CHANNEL, EDA_INPUT_UNIT
    global EDA_DEVICE_RANGE_US, MIN_VALID_SECONDS_BY_CONDITION, EDA_CALIBRATION_SOURCE
    global EDA_SCALE_TO_US, REVIEW_PORT
    if name not in EXPERIMENTS:
        raise ValueError('Unknown experiment: {}. Available: {}'.format(name, ', '.join(EXPERIMENTS)))
    profile = EXPERIMENTS[name]
    scale = EDA_UNIT_TO_US[profile['EDA_INPUT_UNIT']]
    EXPERIMENT_NAME = name
    DATA_DIR, FILE_PATTERN = profile['DATA_DIR'], profile['FILE_PATTERN']
    RAW_BDF_DIR, RAW_BDF_PATTERN = profile['RAW_BDF_DIR'], profile['RAW_BDF_PATTERN']
    ACQUISITION_SYSTEM = profile['ACQUISITION_SYSTEM']
    ECG_CHANNEL, EDA_CHANNEL = profile['ECG_CHANNEL'], profile['EDA_CHANNEL']
    ECG_INPUT_UNIT = profile['ECG_INPUT_UNIT']
    EDA_INPUT_UNIT, EDA_SCALE_TO_US = profile['EDA_INPUT_UNIT'], scale
    EDA_DEVICE_RANGE_US = profile['EDA_DEVICE_RANGE_US']
    MIN_VALID_SECONDS_BY_CONDITION = deepcopy(profile['MIN_VALID_SECONDS_BY_CONDITION'])
    EDA_CALIBRATION_SOURCE = profile['EDA_CALIBRATION_SOURCE']
    REVIEW_PORT = profile['REVIEW_PORT']


def minimum_valid_seconds(signal, record_id):
    """Condition-specific inclusion floor; standard HRV eligibility stays separate."""
    if signal not in ('ecg', 'eda'):
        raise ValueError('Unknown signal: ' + str(signal))
    condition = str(record_id).split('_', 1)[-1]
    fallback = ECG_MIN_TIME_DOMAIN_SECONDS if signal == 'ecg' else EDA_MIN_VALID_SECONDS
    return MIN_VALID_SECONDS_BY_CONDITION.get(condition, {}).get(signal, fallback)


select_experiment(EXPERIMENT_NAME)
