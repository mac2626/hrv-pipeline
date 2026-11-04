#!/usr/bin/env python3
"""Batch Fibion ECG CSV -> SDNN and RMSSD for the first 5 minutes of each clock hour.

Keep this file beside batch_hourly_hrv_movesense.py. Original recordings are
read only. Fibion timestamps retain their local clock for hourly windowing,
while their explicit offsets are respected for elapsed-time and gap checks.
"""
import sys
import batch_hourly_hrv_movesense as pipeline

pipeline.DEVICE = "Fibion"
pipeline.__doc__ = __doc__


if __name__ == '__main__':
    if '--help' in sys.argv or '-h' in sys.argv or '--list-only' in sys.argv:
        raise SystemExit(pipeline.main())
    try:
        import numpy as np
        import pandas as pd
        import neurokit2 as nk
    except ImportError as error:
        raise SystemExit(
            f'Missing dependency: {error}. Run: python -m pip install numpy pandas scipy neurokit2')
    pipeline.np = np
    pipeline.pd = pd
    pipeline.nk = nk
    raise SystemExit(pipeline.main())
