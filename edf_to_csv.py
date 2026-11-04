#!/usr/bin/env python3
"""
Convert any EDF file (e.g. a Movesense ECG recording) to CSV.

Usage:
    python3 edf_to_csv.py "path/to/recording.edf"
    python3 edf_to_csv.py "path/to/recording.edf" -o "path/to/output.csv"

By default the CSV is saved next to the .edf with the same name.

Output columns:
    timestamp  - clock time of each sample (EDF start time + sample offset)
    time_s     - seconds since the start of the recording
    <channel>  - one column per EDF channel (e.g. the ECG signal in mV)

Requires:  pip install pyedflib numpy pandas

The file is read and written in chunks, so very large recordings
(e.g. 27 hours of ECG at 256 Hz, ~25 million rows) do not fill up memory.
If channels have different sampling rates, one CSV is written per channel.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyedflib

CHUNK = 500_000  # samples per read/write
TS_FORMAT = "%Y-%m-%d %H:%M:%S.%f"


def build_chunk(start_dt, fs, first, n):
    """Return (timestamps, seconds) for samples first .. first+n-1."""
    secs = np.arange(first, first + n, dtype="float64") / fs
    offsets = pd.to_timedelta(np.round(secs * 1e6).astype("int64"), unit="us")
    return pd.Timestamp(start_dt) + offsets, secs


def write_csv(path, f, channels, fs, total, start_dt):
    """Stream the given channel indices into one CSV."""
    if os.path.exists(path):
        os.remove(path)
    names = [f.getLabel(i) for i in channels]
    for first in range(0, total, CHUNK):
        n = min(CHUNK, total - first)
        ts, secs = build_chunk(start_dt, fs, first, n)
        df = pd.DataFrame({"timestamp": ts, "time_s": np.round(secs, 6)})
        for i, name in zip(channels, names):
            df[name] = f.readSignal(i, start=first, n=n)
        df.to_csv(
            path,
            mode="w" if first == 0 else "a",
            header=(first == 0),
            index=False,
            date_format=TS_FORMAT,
        )
        done = min(first + n, total)
        print(f"\r  {done / total * 100:5.1f}%  ({done:,} / {total:,} rows)", end="", flush=True)
    print()
    print(f"Wrote {path}")


def convert(edf_path, out_path=None):
    base = os.path.splitext(edf_path)[0]
    out_path = out_path or base + ".csv"

    f = pyedflib.EdfReader(edf_path)
    try:
        n_ch = f.signals_in_file
        labels = f.getSignalLabels()
        rates = [float(f.getSampleFrequency(i)) for i in range(n_ch)]
        totals = [int(n) for n in f.getNSamples()]
        start_dt = f.getStartdatetime()

        print(f"Reading:    {edf_path}")
        print(f"Start time: {start_dt}")
        for i in range(n_ch):
            hours = totals[i] / rates[i] / 3600
            print(f"  ch{i}: {labels[i]!r}  {rates[i]:g} Hz  {totals[i]:,} samples  ({hours:.2f} h)")

        if len(set(rates)) == 1 and len(set(totals)) == 1:
            write_csv(out_path, f, list(range(n_ch)), rates[0], totals[0], start_dt)
        else:
            print("Channels differ in sampling rate/length: writing one CSV per channel.")
            for i in range(n_ch):
                safe = "".join(c if c.isalnum() else "_" for c in labels[i])
                write_csv(f"{base}_{safe}.csv", f, [i], rates[i], totals[i], start_dt)
    finally:
        f.close()

    print("SUCCESS!")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Convert an EDF file to CSV")
    p.add_argument("edf", help="path to the input .edf file")
    p.add_argument("-o", "--output", help="path for the output .csv (optional)")
    args = p.parse_args()
    if not os.path.isfile(args.edf):
        sys.exit(f"File not found: {args.edf}")
    convert(args.edf, args.output)
