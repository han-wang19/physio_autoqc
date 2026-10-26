"""Draft-based QC review: one history and one completion action per record."""
from datetime import datetime
from base64 import b64encode
from types import SimpleNamespace

import numpy as np
import pandas as pd
import streamlit as st
import project_config as cfg
from reports.record_labels import natural_record_key
from review.data import load_quality_rows
from review.quality import reviewed_disposition

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


def apply_chart_event(draft, event, signal, rid, size, sfreq, peak_edit):
    """Consume each chart edit once, including an edit arriving with Finish."""
    if (not event or event.get('record_id', rid) != rid
            or event.get('signal', signal) != signal or draft['event'] == event.get('id')):
        return False
    draft['event'] = event.get('id')
    if isinstance(event.get('viewport'), dict):
        draft['viewport'] = event['viewport']
    if event['kind'] == 'interval' and event.get('label') in ('PASS', 'FAIL'):
        left, right = sorted(map(float, event['bounds']))
        left, right = max(0., left), min(size / sfreq, right)
        if right - left >= .01:
            new = dict(signal=signal, record_id=rid, start_s=left, end_s=right,
                label=event['label'], reviewer='', reason='', notes='',
                reviewed_at=datetime.now().isoformat(timespec='microseconds'))
            return apply_edit(draft, draft['intervals'] + [new], draft['peaks'])
    elif event['kind'] == 'peak' and peak_edit:
        peaks, index = draft['peaks'], int(event['index'])
        new_peaks = (peaks[peaks != index] if event.get('action') == 'DELETE'
                     else np.unique(np.append(peaks, min(size - 1, max(0, index)))))
        return apply_edit(draft, draft['intervals'], new_peaks)
    return False


def complete_review(signal, rid, row, data, draft):
    """Persist edits and the reviewed FIF before marking this record complete."""
    saved = [dict(item, reviewer='') for item in draft['intervals']]
    sfreq = float(data['sfreq'])
    size = data['filtered' if signal == 'ecg' else 'clean'].size
    quality = load_quality_rows(signal, rid, size, sfreq, pd.DataFrame(saved, columns=INTERVAL_COLUMNS))
    disposition = reviewed_disposition(signal, rid, quality, size, sfreq)
    if draft['dirty']:
        save_manual_intervals(signal, rid, saved, derive=False)
        if signal == 'ecg':
            save_manual_peaks(rid, data['primary_peaks'], draft['peaks'], float(data['sfreq']))
    result = dict(signal=signal, record_id=rid, auto_status=row.auto_status, **disposition,
                  reason='REVIEW_COMPLETED_RESOLVED_QC', notes='')
    result.update(reviewer='', reviewed_at=datetime.now().isoformat(timespec='microseconds'))
    from review.fif_export import export_reviewed_signal
    export_reviewed_signal(signal, rid, data=data, peaks=draft['peaks'],
                           decision=SimpleNamespace(**result), quality_rows=quality)
    save_manual_decision(result)
    draft.update(intervals=saved, dirty=False, undo=[])


def run_review(load_report, load_record, ecg_figure, eda_figure, chart, header=None):
    just_finished = st.session_state.pop('review_just_finished', None)
    if just_finished:
        st.toast('保存成功。')
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
    rid = st.sidebar.selectbox('记录', sorted(choices.record_id.astype(str).tolist(), key=natural_record_key),
        format_func=lambda value: '{} | {} | {}'.format(value, statuses[value], '已人工审核' if value in reviewed else '待审核'))
    if header is not None:
        header.markdown(f'**实验名称：{cfg.EXPERIMENT_NAME}　｜　当前记录：{rid}　｜　{signal.upper()}**')
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
    active_record = (cfg.EXPERIMENT_NAME, signal, rid)
    if st.session_state.get('review_active_record') != active_record:
        draft.pop('viewport', None)
        st.session_state['review_active_record'] = active_record
    if signal == 'ecg' and draft.get('source_rule', 'dominant_v1') != source_rule:
        if draft['dirty']:
            st.warning('自动 R 峰对齐规则已更新；当前仍显示你的未保存草稿。请先处理草稿，再重新打开记录核对。')
        else:
            times = load_manual_peak_times(rid)
            peaks = np.rint(times * sfreq).astype(int) if times is not None else data['primary_peaks']
            draft.update(peaks=np.unique(peaks[(peaks >= 0) & (peaks < size)]),
                         source_rule=source_rule, undo=[])
    c1, c2, c3 = st.columns([1, 1, 2])
    peak_edit = c1.toggle('编辑 R 峰', key='peak_mode_'+key) if signal == 'ecg' else False
    if c2.button('↶ 撤销', disabled=not draft['undo'], key='undo_'+key):
        draft['intervals'], draft['peaks'] = draft['undo'].pop()
        draft['dirty'] = True
        st.rerun()
    finish_clicked = c3.button('完成审阅', type='primary', key='finish_'+key)
    if finish_clicked and not just_finished:
        try:
            apply_chart_event(draft, st.session_state.get('review_chart'),
                              signal, rid, size, sfreq, peak_edit)
            complete_review(signal, rid, row, data, draft)
        except Exception as error:
            st.error('保存未完成，请重试：{}'.format(error))
            return
        st.session_state['review_just_finished'] = rid
        st.rerun()
    status_text = '自动结果：{} | {}；总时长 {:.1f} 秒。'.format(row.auto_status, row.reason, duration)
    if pd.notna(row.analysis_start_s) and pd.notna(row.analysis_end_s):
        auto_start, auto_end = float(row.analysis_start_s), float(row.analysis_end_s)
        status_text += ' 自动分析区间 {:.1f}–{:.1f} 秒。'.format(auto_start, auto_end)
    else:
        status_text += ' 自动结果未选定分析区间。'
    intervals = pd.DataFrame(draft['intervals'], columns=INTERVAL_COLUMNS)
    quality = pd.DataFrame(load_quality_rows(signal, rid, size, sfreq, intervals))
    bad = quality[quality.status.ne('PASS')]
    passes = quality.loc[quality.status.eq('PASS'), ['start_s', 'end_s']]
    if signal == 'ecg':
        fig, _, _, _ = ecg_figure(data, draft['peaks'], 0., min(duration, 60.), rid,
                                  bad, passes, intervals, 1., 1., 'pan')
    else:
        fig = eda_figure(data, 0., min(duration, 60.), rid, bad, passes, intervals, 1., 1., 'pan')
    event = chart(figure=fig.to_json(), signal=signal, record_id=rid, mode='peak' if peak_edit else 'interval',
                  viewport=draft.get('viewport'),
                  hr_limits=list(cfg.ECG_PLAUSIBLE_HR_BPM), rr_deviation=cfg.ECG_RR_LOCAL_DEVIATION,
                  sfreq=sfreq,
                  waveform_b64=b64encode(np.asarray(data['filtered'], dtype='<f8').tobytes()).decode('ascii') if signal == 'ecg' else '',
                  peaks=draft['peaks'].tolist(), key='review_chart', default=None)
    st.caption('拖选区间后标记 PASS / FAIL；滚轮缩放，Shift + 拖动平移。' +
               (' 开启「编辑 R 峰」后单击添加、双击删除。' if signal == 'ecg' else ''))
    st.caption(status_text)
    if apply_chart_event(draft, event, signal, rid, size, sfreq, peak_edit):
        st.rerun()
    if draft['dirty']:
        st.warning('有修改尚未保存，请点击「完成审阅」保存当前记录。')
    else:
        st.caption('当前无未保存修改。')
