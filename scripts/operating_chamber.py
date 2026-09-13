#!/usr/bin/env python3
"""Run from a source checkout; numerical and rendering inputs are frozen per output."""
import os
from pathlib import Path
import sys
import tempfile
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'chamber-cfd-mpl'))
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fluid_dynamic.operating.workflow import main
if __name__=='__main__':
    main(ROOT)
