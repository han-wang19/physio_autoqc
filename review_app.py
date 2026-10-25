"""Streamlit entry point. Run with: streamlit run review_app.py."""
import atexit
from pathlib import Path
import os

# Keep JIT caches writable when the Python installation is read-only.
os.environ.setdefault('NUMBA_CACHE_DIR', str(Path(__file__).resolve().parent / '.cache' / 'numba'))

import streamlit as st
import streamlit.components.v1 as components
import project_config as cfg

from review.manual_qc import refresh_dirty_summaries
from review.data import (
    load_report, load_record, load_window_qc, load_auto_pass_intervals,
)
from review.figures import ecg_figure, eda_figure
from review.workflow import run_review


def main():
    st.set_page_config(page_title=f"{cfg.EXPERIMENT_NAME} | Physio AutoQC Review", layout="wide")
    st.caption(f"实验名称：{cfg.EXPERIMENT_NAME}")
    st.sidebar.caption(f"实验名称：{cfg.EXPERIMENT_NAME}")
    st.title("ECG / EDA 人工复核")
    chart = components.declare_component(
        "physio_review_chart", path=str(Path(__file__).parent / "review" / "chart"),
    )
    if not getattr(atexit, "_physio_qc_registered", False):
        atexit.register(refresh_dirty_summaries)
        atexit._physio_qc_registered = True
    if st.sidebar.button("结束本次审核并更新汇总", type="primary"):
        with st.spinner("正在更新人工审核记录的 QC 图、最终特征及汇总图……"):
            refresh_dirty_summaries(update_records=True)
        st.sidebar.success("人工审核记录的 QC 图、最终处理结果、有效汇总表和覆盖图已更新；最终检查未通过的记录见处理报告。")
    try:
        run_review(load_report, load_record, load_window_qc, load_auto_pass_intervals,
                   ecg_figure, eda_figure, chart)
    except FileNotFoundError as error:
        st.error(str(error))
        st.stop()


if __name__ == "__main__":
    main()
