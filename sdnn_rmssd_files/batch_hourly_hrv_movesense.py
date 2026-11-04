#!/usr/bin/env python3
"""Batch Movesense ECG CSV -> SDNN and RMSSD for the first 5 minutes of each clock hour.

Original recordings are read only. Run with --help. This script uses the local
wall-clock timestamps written by the device; it does not infer a timezone.
"""
from pathlib import Path
import argparse
import csv
import importlib.metadata
import json
import sys
from datetime import datetime, timezone

DEVICE = "Movesense"
VERSION = "2.0.0"
WINDOW_SECONDS = 300
HOUR_NS = 3_600_000_000_000
SECOND_NS = 1_000_000_000


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', help='Folder of ECG CSVs, searched recursively')
    p.add_argument('--output', help='New output folder (must not already exist)')
    p.add_argument('--pattern', default='*.csv', help='Recursive filename pattern')
    p.add_argument('--fs', type=int, default=256, help='Actual acquisition rate, Hz (default 256)')
    p.add_argument('--ecg-column', help='Override the device-specific ECG column name')
    p.add_argument('--time-column', default='timestamp', help='Sample timestamp column')
    p.add_argument('--delimiter', default=',')
    p.add_argument('--skiprows', type=int, default=0, help='Metadata lines before header')
    p.add_argument('--rr-min-ms', type=float, default=300)
    p.add_argument('--rr-max-ms', type=float, default=2000)
    p.add_argument('--relative-threshold', type=float, default=.20)
    p.add_argument('--min-kept-percent', type=float, default=80)
    p.add_argument('--invert', action='store_true', help='Invert ECG before R-peak detection')
    p.add_argument('--list-only', action='store_true', help='List candidate files without analysis')
    a = p.parse_args()
    if not a.input:
        a.input = input(f'Paste the {DEVICE} ECG folder path, then press Return: ').strip().strip("\"'")
    a.input = Path(a.input).expanduser().resolve()
    if not a.input.is_dir():
        p.error(f'Input folder does not exist: {a.input}')
    if not a.list_only and not 100 <= a.fs <= 4000:
        p.error('Supported sampling rates: 100–4000 Hz; use the actual recording rate.')
    if DEVICE == 'Fibion' and a.fs != 256:
        p.error('The Fibion script requires the fixed 256 Hz acquisition rate.')
    if not 0 < a.rr_min_ms < a.rr_max_ms:
        p.error('RR limits must be positive and increasing.')
    if not 0 < a.relative_threshold < 1 or not 0 <= a.min_kept_percent <= 100:
        p.error('Invalid screening thresholds.')
    a.output = Path(a.output).expanduser().resolve() if a.output else (
        a.input / ('hourly_hrv_results_' + DEVICE.lower() + '_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f')))
    return a


def choose_columns(a, columns):
    expected = a.ecg_column or ('ECG' if DEVICE == 'Movesense' else 'ecg')
    if expected not in columns:
        if any(c.lower() == 'ecg' for c in columns):
            raise ValueError(f'Wrong device structure: {DEVICE} requires {expected!r}; found {columns}')
        return None, None
    if a.time_column not in columns:
        raise ValueError(f'Missing timestamp column {a.time_column!r}; found {columns}')
    return expected, a.time_column


def parse_time(values):
    """Return (absolute_ns, local_wall_clock_ns, display_strings)."""
    s = values.astype('string').str.strip()
    if DEVICE == 'Movesense':
        valid = s.str.fullmatch(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{1,9}', na=False)
        if not valid.all():
            raise ValueError('Expected Movesense timestamps like 2026-10-04 20:48:05.003906250')
        wall = pd.to_datetime(s, format='%Y-%m-%d %H:%M:%S.%f', errors='coerce')
        absolute = wall
    else:
        valid = s.str.fullmatch(r'\d{4}-\d{1,2}-\d{1,2}T\d{2}:\d{2}:\d{2}\.\d{1,3}(?:[+-]\d{2}:\d{2}|Z)', na=False)
        if not valid.all():
            raise ValueError('Expected Fibion timestamps like 2026-10-5T16:37:54.004-04:00')
        normalized = s.str.replace(
            r'^(\d{4})-(\d{1,2})-(\d{1,2})T',
            lambda m: f'{m[1]}-{int(m[2]):02d}-{int(m[3]):02d}T', regex=True)
        absolute = pd.to_datetime(normalized, format='ISO8601', utc=True, errors='coerce')
        # Strip the explicit offset to retain the device's local clock for hour anchoring.
        wall_text = normalized.str.replace(r'(?:[+-]\d{2}:\d{2}|Z)$', '', regex=True)
        wall = pd.to_datetime(wall_text, format='%Y-%m-%dT%H:%M:%S.%f', errors='coerce')
    if absolute.isna().any() or wall.isna().any():
        raise ValueError('Unparseable or missing timestamps; no rows were silently removed.')
    abs_ns = absolute.dt.as_unit('ns').astype('int64').to_numpy(dtype='int64')
    wall_ns = wall.dt.as_unit('ns').astype('int64').to_numpy(dtype='int64')
    return abs_ns, wall_ns


def hour_label(hour_ns):
    stamp = pd.Timestamp(hour_ns, unit='ns')
    return stamp.strftime('%Y-%m-%d %H:00:00')


def screen_rr(rr, a):
    valid = (rr >= a.rr_min_ms) & (rr <= a.rr_max_ms)
    med = pd.Series(np.where(valid, rr, np.nan)).rolling(
        11, center=True, min_periods=1).median().to_numpy()
    return valid & (np.abs(rr-med) <= a.relative_threshold * med)


def detect_window(a, x, absolute_ns, wall_ns, row_numbers, hour_ns):
    """Detect beats with padding; assign intervals by the second peak's wall time."""
    valid = np.isfinite(x)
    cut = np.zeros(len(x), dtype=bool)
    cut[0] = True
    dt = np.diff(absolute_ns) / SECOND_NS
    cut[1:] |= (~np.isfinite(dt)) | (dt <= 0) | (dt > 1.5/a.fs)
    cut[1:] |= ~valid[:-1] | ~valid[1:]
    edges = np.unique(np.r_[0, np.flatnonzero(cut), len(x)])
    all_rr, all_keep, all_rows, all_wall, rmssd_diffs = [], [], [], [], []
    failures = 0
    core_stop = hour_ns + WINDOW_SECONDS*SECOND_NS
    for left, right in zip(edges[:-1], edges[1:]):
        if not valid[left] or right-left < 10*a.fs:
            continue
        signal = x[left:right]
        if np.ptp(signal) == 0:
            continue
        try:
            clean = nk.ecg_clean(-signal if a.invert else signal,
                                  sampling_rate=a.fs, method='neurokit')
            _, info = nk.ecg_peaks(clean, sampling_rate=a.fs,
                                   method='neurokit', correct_artifacts=False)
            peaks = np.asarray(info['ECG_R_Peaks'], dtype=int) + left
            peaks = peaks[(peaks >= left+int(.5*a.fs)) &
                          (peaks < right-int(.5*a.fs))]
            rr = np.diff(peaks).astype(float) * 1000 / a.fs
            keep = screen_rr(rr, a)
            end_wall = wall_ns[peaks[1:]]
            belongs = (end_wall >= hour_ns) & (end_wall < core_stop)
            all_rr.extend(rr[belongs])
            all_keep.extend(keep[belongs])
            all_rows.extend(row_numbers[peaks[1:]][belongs])
            all_wall.extend(end_wall[belongs])
            # RMSSD uses only adjacent, retained intervals that both belong to the core.
            pair = keep[:-1] & keep[1:] & belongs[:-1] & belongs[1:]
            rmssd_diffs.extend(np.diff(rr)[pair])
        except Exception:
            failures += 1
    return (np.asarray(all_rr), np.asarray(all_keep, dtype=bool),
            np.asarray(all_rows, dtype=np.int64), np.asarray(all_wall, dtype=np.int64),
            np.asarray(rmssd_diffs, dtype=float), failures)


def analyze_hour(a, hour_ns, pieces, relative, output_writer, rr_writer):
    x = np.concatenate([p[0] for p in pieces])
    absolute_ns = np.concatenate([p[1] for p in pieces])
    wall_ns = np.concatenate([p[2] for p in pieces])
    rows = np.concatenate([p[3] for p in pieces])
    order = np.argsort(rows, kind='stable')
    x, absolute_ns, wall_ns, rows = x[order], absolute_ns[order], wall_ns[order], rows[order]
    unique = np.r_[True, rows[1:] != rows[:-1]]
    x, absolute_ns, wall_ns, rows = x[unique], absolute_ns[unique], wall_ns[unique], rows[unique]
    core = (wall_ns >= hour_ns) & (wall_ns < hour_ns + WINDOW_SECONDS*SECOND_NS)
    core_t = absolute_ns[core]
    core_x = x[core]
    dt = np.diff(core_t) / SECOND_NS
    gaps = dt[np.isfinite(dt) & (dt > 1.5/a.fs)]
    missing = int(np.sum(np.maximum(0, np.rint(gaps*a.fs)-1)))
    rr, keep, rr_rows, rr_wall, rmssd_diffs, failures = detect_window(
        a, x, absolute_ns, wall_ns, rows, hour_ns)
    retained = rr[keep]
    percent = 100*len(retained)/len(rr) if len(rr) else 0.
    coverage = float(np.sum(retained)/1000/WINDOW_SECONDS*100)
    tolerance_ns = int(1.5/a.fs*SECOND_NS)
    complete_start = len(core_t) and wall_ns[core][0] <= hour_ns+tolerance_ns
    complete_end = len(core_t) and wall_ns[core][-1] >= hour_ns+WINDOW_SECONDS*SECOND_NS-tolerance_ns
    flags = []
    if not (complete_start and complete_end): flags.append('INCOMPLETE_CLOCK_WINDOW')
    if len(core_x) < int(WINDOW_SECONDS*a.fs*.99): flags.append('LOW_SAMPLE_COUNT')
    if len(retained) < 30: flags.append('FEW_INTERVALS')
    if percent < a.min_kept_percent: flags.append('LOW_RETENTION')
    if coverage < 90: flags.append('LOW_RR_COVERAGE')
    if np.any(~np.isfinite(core_x)): flags.append('INVALID_ECG')
    if len(gaps): flags.append('TIMESTAMP_GAPS')
    if failures: flags.append('DETECTION_FAILURE')
    sdnn = float(np.std(retained, ddof=1)) if len(retained) > 1 else ''
    rmssd = float(np.sqrt(np.mean(rmssd_diffs**2))) if len(rmssd_diffs) else ''
    row = dict(file=relative, hour_start_local=hour_label(hour_ns),
        window_end_local=hour_label(hour_ns)[:11] + pd.Timestamp(
            hour_ns+WINDOW_SECONDS*SECOND_NS, unit='ns').strftime('%H:%M:%S'),
        fs_hz=a.fs, core_samples=len(core_x), window_seconds=WINDOW_SECONDS,
        sdnn_screened_ms=sdnn, rmssd_screened_ms=rmssd,
        rr_detected=len(rr), rr_kept=len(retained), rmssd_pairs=len(rmssd_diffs),
        percent_nn_kept=percent, retained_rr_coverage_percent=coverage,
        missing_samples_estimate=missing, timestamp_gaps=len(gaps),
        invalid_ecg_rows=int(np.sum(~np.isfinite(core_x))),
        detection_failures=failures, status='REVIEW' if flags else 'COMPUTED',
        flags=';'.join(flags))
    output_writer.writerow(row)
    for interval, accepted, row_number, end_wall in zip(rr, keep, rr_rows, rr_wall):
        rr_writer.writerow(dict(file=relative, hour_start_local=hour_label(hour_ns),
            rr_end_row=int(row_number), rr_end_local=pd.Timestamp(end_wall, unit='ns').isoformat(),
            rr_ms=float(interval), kept=int(accepted)))
    return row


def analyze_file(path, a, relative, output_writer, rr_writer):
    header = pd.read_csv(path, sep=a.delimiter, skiprows=a.skiprows, nrows=0,
                         encoding='utf-8-sig')
    ecg, time = choose_columns(a, list(header.columns))
    if ecg is None:
        return dict(file=relative, status='SKIPPED', message='No ECG column',
                    hours=0, computed_hours=0, review_hours=0)
    chunk_rows = max(WINDOW_SECONDS*a.fs, 100_000)
    pad_seconds = max(3., a.rr_max_ms/1000+1)
    pad_ns = int(pad_seconds*SECOND_NS)
    active, finished = {}, set()
    tail = None
    last_abs = last_wall = None
    row_offset = 0
    results = []
    reader = pd.read_csv(path, sep=a.delimiter, skiprows=a.skiprows,
                         encoding='utf-8-sig', usecols=[ecg, time], chunksize=chunk_rows)
    for df in reader:
        x = pd.to_numeric(df[ecg], errors='coerce').to_numpy(dtype=float)
        absolute_ns, wall_ns = parse_time(df[time])
        rows = np.arange(row_offset, row_offset+len(df), dtype=np.int64)
        row_offset += len(df)
        if last_abs is not None and (absolute_ns[0] <= last_abs or wall_ns[0] <= last_wall):
            raise ValueError('Duplicate/decreasing timestamps or a local-clock reversal; split DST-transition files before analysis.')
        if np.any(np.diff(absolute_ns) <= 0) or np.any(np.diff(wall_ns) <= 0):
            raise ValueError('Duplicate or decreasing timestamps; repair timestamp interpretation before analysis.')
        dt = np.diff(np.r_[[last_abs], absolute_ns] if last_abs is not None else absolute_ns)/SECOND_NS
        positive = dt[dt > 0]
        if len(positive) and abs(float(np.median(positive))*a.fs-1) > .15:
            raise ValueError('Timestamp spacing disagrees with --fs; verify recording settings.')
        current = (x, absolute_ns, wall_ns, rows)
        candidate_hours = np.unique((wall_ns[(wall_ns % HOUR_NS) < WINDOW_SECONDS*SECOND_NS] // HOUR_NS)*HOUR_NS)
        for hour_ns in candidate_hours:
            hour_ns = int(hour_ns)
            if hour_ns in finished or hour_ns in active:
                continue
            active[hour_ns] = []
            if tail is not None:
                mask = (tail[2] >= hour_ns-pad_ns) & (tail[2] < hour_ns+(WINDOW_SECONDS*SECOND_NS)+pad_ns)
                if np.any(mask): active[hour_ns].append(tuple(v[mask] for v in tail))
        for hour_ns, pieces in active.items():
            mask = (wall_ns >= hour_ns-pad_ns) & (wall_ns < hour_ns+WINDOW_SECONDS*SECOND_NS+pad_ns)
            if np.any(mask): pieces.append(tuple(v[mask] for v in current))
        for hour_ns in [h for h in active if wall_ns[-1] >= h+WINDOW_SECONDS*SECOND_NS+pad_ns]:
            results.append(analyze_hour(a, hour_ns, active.pop(hour_ns), relative,
                                        output_writer, rr_writer))
            finished.add(hour_ns)
        keep_tail = wall_ns >= wall_ns[-1]-pad_ns
        tail = tuple(v[keep_tail] for v in current)
        last_abs, last_wall = int(absolute_ns[-1]), int(wall_ns[-1])
    for hour_ns in sorted(active):
        results.append(analyze_hour(a, hour_ns, active[hour_ns], relative,
                                    output_writer, rr_writer))
    if row_offset == 0:
        raise ValueError('File contains no sample rows.')
    computed = sum(r['status'] == 'COMPUTED' for r in results)
    review = sum(r['status'] == 'REVIEW' for r in results)
    status = 'REVIEW' if review or not results else 'COMPUTED'
    message = ('No samples occurred in the first five minutes of any clock hour.'
               if not results else 'Automated screened RR; normal sinus beats not confirmed')
    return dict(file=relative, status=status, message=message, ecg_column=ecg,
                time_column=time, fs_hz=a.fs, sample_rows=row_offset,
                hours=len(results), computed_hours=computed, review_hours=review)


def writer(handle, fields):
    w = csv.DictWriter(handle, fieldnames=fields.split(), extrasaction='raise')
    w.writeheader()
    return w


def main():
    a = arguments()
    paths = sorted(p for p in a.input.rglob(a.pattern)
                   if p.is_file()
                   and not any(part.startswith(('sdnn_results_', 'hourly_hrv_results_'))
                               for part in p.relative_to(a.input).parts)
                   and a.output not in p.parents)
    if a.list_only:
        for p in paths: print(p)
        print(f'{len(paths)} candidates; files without the expected ECG column will be skipped.')
        return 0
    if not paths:
        print('No matching CSV files. Check the folder and --pattern.', file=sys.stderr)
        return 1
    a.output.mkdir(parents=True, exist_ok=False)
    manifest = {k:str(v) if isinstance(v, Path) else v for k,v in vars(a).items()}
    manifest.update(script_version=VERSION, device=DEVICE,
        window_definition='first 300 seconds of each local clock hour',
        created_utc=datetime.now(timezone.utc).isoformat(),
        packages={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','neurokit2']},
        python=sys.version, files=[str(p.relative_to(a.input)) for p in paths])
    (a.output/'run_settings.json').write_text(json.dumps(manifest, indent=2))
    counts = dict(COMPUTED=0, REVIEW=0, ERROR=0, SKIPPED=0)
    with (a.output/'file_summary.csv').open('w', newline='') as f, \
         (a.output/'hourly_hrv.csv').open('w', newline='') as h, \
         (a.output/'rr_audit.csv').open('w', newline='') as r:
        summaries = writer(f, 'file status message ecg_column time_column fs_hz sample_rows hours computed_hours review_hours')
        hours = writer(h, 'file hour_start_local window_end_local fs_hz core_samples window_seconds sdnn_screened_ms rmssd_screened_ms rr_detected rr_kept rmssd_pairs percent_nn_kept retained_rr_coverage_percent missing_samples_estimate timestamp_gaps invalid_ecg_rows detection_failures status flags')
        intervals = writer(r, 'file hour_start_local rr_end_row rr_end_local rr_ms kept')
        for i, path in enumerate(paths, 1):
            relative = str(path.relative_to(a.input))
            print(f'[{i}/{len(paths)}] {relative}', flush=True)
            try:
                result = analyze_file(path, a, relative, hours, intervals)
            except Exception as error:
                result = dict(file=relative, status='ERROR', message=f'{type(error).__name__}: {error}')
            summaries.writerow(result)
            counts[result['status']] += 1
            for handle in (f, h, r): handle.flush()
            print('  '+result['status']+'  hourly windows: '+str(result.get('hours', '')), flush=True)
    print(f'Results: {a.output}\n{counts}')
    print('Compare only matching, complete clock windows after reviewing flags and device clock alignment.')
    return 1 if counts['ERROR'] or counts['COMPUTED']+counts['REVIEW'] == 0 else 0


if __name__ == '__main__':
    if '--help' in sys.argv or '-h' in sys.argv or '--list-only' in sys.argv:
        raise SystemExit(main())
    try:
        import numpy as np
        import pandas as pd
        import neurokit2 as nk
    except ImportError as error:
        raise SystemExit(f'Missing dependency: {error}. Run: python -m pip install numpy pandas scipy neurokit2')
    raise SystemExit(main())
