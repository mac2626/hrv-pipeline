HOURLY SDNN + RMSSD: MOVESENSE AND FIBION
Version 2.0.0 | Prepared 2026-11-04

PURPOSE
These scripts calculate SDNN and RMSSD from the first five minutes of every
local clock hour present in each ECG CSV. For example, the 10:00 hour uses
10:00:00.000 through (but not including) 10:05:00.000. Windows are NOT
anchored to the first sample in a file.

FILES
  batch_hourly_hrv_movesense.py  Main analysis engine and Movesense entry point
  batch_hourly_hrv_fibion.py     Fibion entry point; keep beside the first file
  requirements_hourly_hrv.txt    Python package requirements

Keep the two Python files together in:
  $HOME/Desktop/CLEARS/HRV Data QC/hrv-pipeline

SETUP ON YOUR MAC
In Terminal:
  cd "$HOME/Desktop/CLEARS/HRV Data QC/hrv-pipeline"
  source "$HOME/python-datalogger-tool/.venv/bin/activate"
  python -m pip install -r requirements_hourly_hrv.txt

RUN, ONE DEVICE AT A TIME
  python batch_hourly_hrv_movesense.py
  python batch_hourly_hrv_fibion.py

Each script asks for that device's ECG folder and searches all subfolders.
Use a separate input folder for each device. You may instead specify a path:
  python batch_hourly_hrv_movesense.py --input "/full/path/to/Movesense folder"
  python batch_hourly_hrv_fibion.py --input "/full/path/to/Fibion folder"

DEVICE RULES RETAINED FROM THE EARLIER BATCH SCRIPTS
Movesense requires columns timestamp and uppercase ECG. Its timestamps are
naive local wall-clock values with up to nine fractional digits; no timezone
is assigned or inferred. Default acquisition rate is 256 Hz.

Fibion requires timestamp and lowercase ecg. Unpadded dates and explicit
timezone offsets are supported. Acquisition rate is fixed at 256 Hz despite
rounded millisecond timestamps. The explicit offset is used for elapsed-time
and gap checks, while the written local clock is used to select each hour.

OUTPUTS
Each run creates a new hourly_hrv_results_DEVICE_DATE_TIME folder inside the
input folder. Original CSVs are never changed.

hourly_hrv.csv
  One row per first-five-minute clock window. Key fields are:
  - hour_start_local: clock hour used for window selection
  - sdnn_screened_ms: sample SD of retained RR intervals (n-1 denominator)
  - rmssd_screened_ms: root mean square of successive RR differences
  - rmssd_pairs: number of adjacent retained interval pairs used
  - status and flags: automated data-quality review information

file_summary.csv
  One row per source CSV with counts of computed and review windows.

rr_audit.csv
  Every detected interval assigned to a clock window, its source row, end
  time, duration and retained/excluded status.

run_settings.json
  Input, settings, script/package versions and run timestamp.

METHOD
The ECG is cleaned with NeuroKit2 ecg_clean(method='neurokit'); R peaks are
found with ecg_peaks(method='neurokit', correct_artifacts=False). No samples,
beats or excluded intervals are interpolated.

The signal is split at nonfinite ECG values and timestamp discontinuities.
RR intervals crossing a split are never calculated. Intervals are retained
when they are 300-2000 ms inclusive and within 20% of a centered rolling
median of up to 11 intervals. These defaults match the earlier batch scripts.

SDNN is calculated from all retained intervals whose second R peak is in the
five-minute window. RMSSD is calculated only when both successive intervals
are retained, adjacent in the original detected sequence, and assigned to the
same five-minute window. It never bridges across a rejected interval or gap.

QUALITY FLAGS
REVIEW does not automatically mean unusable; inspect the flags and ECG.
Flags identify incomplete clock windows, low sample count, fewer than 30
retained intervals, <80% retention, <90% retained-RR coverage, invalid ECG,
timestamp gaps or detector failures. A window that starts after HH:00 or ends
before HH:05 is labeled INCOMPLETE_CLOCK_WINDOW.

COMPARISON RULES
Compare Movesense and Fibion only when hour_start_local identifies the same
real recording interval and both device clocks were synchronized. The script
does not infer a Movesense timezone and does not automatically merge devices.
Use only complete windows after reviewing flags and R-peak locations.

Automated screening does not prove that retained beats are normal sinus beats.
The screened labels therefore must not be interpreted as clinical validation.
