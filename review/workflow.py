"""Draft-based QC review: one history and one completion action per record."""
from datetime import datetime

import numpy as np
import pandas as pd
import streamlit as st
import project_config as cfg

from review.manual_qc import (
    INTERVAL_COLUMNS, read_manual_decisions, read_manual_intervals,
    load_manual_peak_times, save_manual_intervals, save_manual_peaks,
    save_manual_decision,
)


def candidate_statuses(report, completed):
    """Show saved review outcomes, falling back to automatic QC for pending records."""
    statuses = dict(zip(report.record_id.astype(str), report.auto_status))
    status_map = {'INCLUDE_FULL': 'PASS', 'INCLUDE_TRIMMED': 'TRIMMED',
                  'EXCLUDE': 'FAIL', 'SECOND_REVIEW': 'REVIEW'}
    for item in completed.itertuples(index=False):
        statuses[str(item.record_id)] = status_map.get(item.manual_status, 'REVIEW')
    return statuses


def snapshot(draft):
    return ([dict(row) for row in draft['intervals']], draft['peaks'].copy())


def apply_edit(draft, intervals, peaks):
    if intervals == draft['intervals'] and np.array_equal(peaks, draft['peaks']):
        return False
    draft['undo'] = (draft['undo'] + [snapshot(draft)])[-50:]
    draft['intervals'] = intervals
    draft['peaks'] = np.asarray(peaks, dtype=int)
    draft['dirty'] = True
    return True


def run_review(load_report, load_record, load_windows, load_pass, ecg_figure,
               eda_figure, chart):
    just_finished = st.session_state.pop('review_just_finished', None)
    if just_finished:
        st.success('{}：区间、R 峰和审核状态已保存，已人工审核。'.format(just_finished))
    signal = st.sidebar.radio('信号', ['ecg', 'eda'], format_func=str.upper)
    report = load_report(signal)
    decisions = read_manual_decisions()
    completed = decisions[(decisions.signal == signal) & (
        decisions.reason.astype(str).str.startswith('REVIEW_COMPLETED')
        | decisions.reason.eq('ACCEPTED_AUTOMATIC_RESULT'))]
    reviewed = set(completed.record_id.astype(str))
    total = len(report)
    done = int(report.record_id.astype(str).isin(reviewed).sum())
    st.sidebar.progress(done / total if total else 0., text='已审核 {} / {} 条'.format(done, total))
    queue = st.sidebar.radio('审核队列', ['待人工审核', '已人工审核', '全部记录'])
    is_reviewed = report.record_id.astype(str).isin(reviewed)
    choices = report[~is_reviewed] if queue == '待人工审核' else report[is_reviewed] if queue == '已人工审核' else report
    st.sidebar.caption('待人工审核包含所有尚未点击完成审阅的记录，包括自动 PASS。修改暂存于当前会话；点击完成审阅后保存并移入已人工审核。')
    search = st.sidebar.text_input('查找记录', placeholder='输入被试或曲目编号').strip()
    if search:
        choices = choices[choices.record_id.astype(str).str.contains(search, case=False, regex=False)]
    if choices.empty:
        st.info('当前队列为空，请切换左侧筛选。')
        return
    statuses = candidate_statuses(report, completed)
    rid = st.sidebar.selectbox('记录', choices.record_id.astype(str).tolist(),
        format_func=lambda value: '{} | {} | {}'.format(value, statuses[value], '已人工审核' if value in reviewed else '待审核'))
    row = report[report.record_id.astype(str).eq(rid)].iloc[0]
    data = load_record(signal, rid)
    sfreq = float(data['sfreq'])
    size = data['filtered'].size if signal == 'ecg' else data['clean'].size
    duration = size / sfreq
    source_rule = str(data.get('peak_alignment_rule', 'dominant_v1'))
    key = 'review_draft_{}_{}'.format(signal, rid)
    if key not in st.session_state:
        times = load_manual_peak_times(rid) if signal == 'ecg' else None
        peaks = np.rint(times * sfreq).astype(int) if times is not None else data.get('primary_peaks', np.array([], dtype=int))
        st.session_state[key] = dict(intervals=read_manual_intervals(signal, rid).to_dict('records'),
            peaks=np.unique(peaks[(peaks >= 0) & (peaks < size)]), undo=[], dirty=False, event=None,
            source_rule=source_rule)
    draft = st.session_state[key]
    if signal == 'ecg' and draft.get('source_rule', 'dominant_v1') != source_rule:
        if draft['dirty']:
            st.warning('自动 R 峰对齐规则已更新；当前仍显示你的未保存草稿。请先处理草稿，再重新打开记录核对。')
        else:
            times = load_manual_peak_times(rid)
            peaks = np.rint(times * sfreq).astype(int) if times is not None else data['primary_peaks']
            draft.update(peaks=np.unique(peaks[(peaks >= 0) & (peaks < size)]),
                         source_rule=source_rule, undo=[])
    reviewer = st.text_input('审阅者', key='reviewer', placeholder='完成审阅前填写')
    c1, c2, c3 = st.columns([1, 1, 2])
    peak_edit = c1.toggle('编辑 R 峰', key='peak_mode_'+key) if signal == 'ecg' else False
    if c2.button('↶ 撤销', disabled=not draft['undo'], key='undo_'+key):
        draft['intervals'], draft['peaks'] = draft['undo'].pop()
        draft['dirty'] = True
        st.rerun()
    st.caption('拖选时间段后选择 PASS / FAIL / 取消。滚轮缩放时间轴，两指横滑平移；按住 Shift 拖动也可平移。' +
               (' 开启 R 峰编辑：悬浮显示候选红点，单击添加，双击已有红点删除。' if signal == 'ecg' else ''))
    st.caption('自动结果：{} | {}。记录总时长：{:.1f} 秒。'.format(row.auto_status, row.reason, duration))
    if pd.notna(row.analysis_start_s) and pd.notna(row.analysis_end_s):
        auto_start, auto_end = float(row.analysis_start_s), float(row.analysis_end_s)
        st.caption('自动选定分析区间：{:.1f}–{:.1f} 秒（共 {:.1f} 秒）。'.format(
            auto_start, auto_end, auto_end-auto_start))
    else:
        st.caption('自动结果未选定分析区间。')
    st.caption('底色显示窗口级自动判断与当前人工修改合成的状态；绿色区段可能有多段，自动选定分析区间见上方。')
    intervals = pd.DataFrame(draft['intervals'], columns=INTERVAL_COLUMNS)
    bad, passes = load_windows(signal, rid), load_pass(signal, rid)
    start = float(row.analysis_start_s) if pd.notna(row.analysis_start_s) else 0.
    start = min(max(0., start), max(0., duration-20))
    if signal == 'ecg':
        fig, _, _, _ = ecg_figure(data, draft['peaks'], start, min(duration,start+20), rid,
                                  bad, passes, intervals, 1., 1., 'pan')
    else:
        fig = eda_figure(data, start, min(duration,start+30), rid, bad, passes, intervals, 1., 1., 'pan')
    event = chart(figure=fig.to_json(), signal=signal, mode='peak' if peak_edit else 'interval',
                  viewport=draft.get('viewport'),
                  hr_limits=list(cfg.ECG_PLAUSIBLE_HR_BPM), rr_deviation=cfg.ECG_RR_LOCAL_DEVIATION,
                  sfreq=sfreq, waveform=data['filtered'].tolist() if signal == 'ecg' else [],
                  peaks=draft['peaks'].tolist(), key='draft_chart_'+key, default=None)
    if event and draft['event'] != event.get('id'):
        draft['event'] = event.get('id')
        if isinstance(event.get('viewport'), dict):
            draft['viewport'] = event['viewport']
        if event['kind'] == 'interval' and event.get('label') in ('PASS','FAIL'):
            left, right = sorted(map(float,event['bounds']))
            left, right = max(0.,left), min(duration,right)
            if right-left >= .01:
                new = dict(signal=signal, record_id=rid, start_s=left, end_s=right,
                    label=event['label'], reviewer=reviewer.strip(), reason='', notes='',
                    reviewed_at=datetime.now().isoformat(timespec='microseconds'))
                apply_edit(draft, draft['intervals']+[new], draft['peaks'])
                st.rerun()
        elif event['kind'] == 'peak' and peak_edit:
            peaks = draft['peaks']
            index = int(event['index'])
            if event.get('action') == 'DELETE':
                new_peaks = peaks[peaks != index]
            else:
                new_peaks = np.unique(np.append(peaks, min(size-1,max(0,index))))
            if apply_edit(draft, draft['intervals'], new_peaks):
                st.rerun()
    if draft['dirty']:
        st.warning('有修改尚未保存，请点击「完成审阅」保存当前记录。')
    else:
        st.caption('当前无未保存修改。')
    finish_clicked = c3.button('完成审阅', type='primary', key='finish_'+key, disabled=not reviewer.strip())
    if finish_clicked and not just_finished:
        saved = []
        for item in draft['intervals']:
            saved.append(dict(item, reviewer=item.get('reviewer') or reviewer.strip()))
        save_manual_intervals(signal, rid, saved)
        if signal == 'ecg':
            save_manual_peaks(rid, data['primary_peaks'], draft['peaks'], sfreq)
        current = read_manual_decisions()
        found = current[(current.signal == signal) & current.record_id.eq(rid)]
        if saved and not found.empty:
            result = found.iloc[-1].to_dict()
            result['reason'] = 'REVIEW_COMPLETED_INTERVALS'
        else:
            result = dict(signal=signal, record_id=rid, auto_status=row.auto_status,
                manual_status={'PASS':'INCLUDE_FULL','TRIMMED':'INCLUDE_TRIMMED','FAIL':'EXCLUDE'}.get(row.auto_status,'SECOND_REVIEW'),
                analysis_start_s=row.analysis_start_s, analysis_end_s=row.analysis_end_s,
                reason='REVIEW_COMPLETED_AUTOMATIC', notes='')
        result.update(reviewer=reviewer.strip(), reviewed_at=datetime.now().isoformat(timespec='microseconds'))
        save_manual_decision(result)
        from review.fif_export import export_reviewed_signal
        export_reviewed_signal(signal, rid)
        draft.update(intervals=saved, dirty=False, undo=[])
        st.session_state['review_just_finished'] = rid
        st.toast('区间、R 峰和审核状态已保存。')
        st.rerun()
