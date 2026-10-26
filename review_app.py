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
from review.data import load_report, load_record
from review.figures import ecg_figure, eda_figure
from review.fif_export import warm_review_export
from review.workflow import run_review


def main():
    st.set_page_config(page_title=f"{cfg.EXPERIMENT_NAME} | Physio AutoQC Review", layout="wide")
    warm_review_export()
    st.markdown('''<style>
        .stMainBlockContainer {padding-top: 4rem; padding-bottom: 1rem;}
        </style>''', unsafe_allow_html=True)
    record_header = st.empty()
    record_header.markdown(f"**实验名称：{cfg.EXPERIMENT_NAME}** · ECG / EDA 复核")
    chart = components.declare_component(
        "physio_review_chart", path=str(Path(__file__).parent / "review" / "chart"),
    )
    if not getattr(atexit, "_physio_qc_registered", False):
        atexit.register(refresh_dirty_summaries)
        atexit._physio_qc_registered = True
    if st.sidebar.button("结束本次审核并更新汇总", type="primary"):
        progress_area = st.sidebar.empty()

        def show_progress(signal, record_id, done, total, stage):
            text = f'{signal.upper()} · {done}/{total} 已完成 · {record_id}：{stage}' if total else f'{signal.upper()} · {stage}'
            progress_area.progress(done / total if total else 1., text=text)

        progress_area.progress(0., text='检查需要更新的记录……')
        try:
            result = refresh_dirty_summaries(update_records=True, progress=show_progress)
        except Exception as error:
            progress_area.empty()
            st.sidebar.error(f'更新未完成：{error}。已完成记录已保存，再次点击可接着处理。')
        else:
            progress_area.empty()
            if result['summaries']:
                st.sidebar.success('已更新 {} 条，跳过 {} 条未变更记录，汇总已刷新。'.format(result['updated'], result['skipped']))
            else:
                st.sidebar.success('当前结果已是最新，无需更新。')
    try:
        run_review(load_report, load_record, ecg_figure, eda_figure, chart, header=record_header)
    except FileNotFoundError as error:
        st.error(str(error))
        st.stop()


if __name__ == "__main__":
    main()
