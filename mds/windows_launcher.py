"""Absolute-path entry point for Windows shortcuts and file associations."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))

if __name__=='__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    from mds.app import main
    raise SystemExit(main())
