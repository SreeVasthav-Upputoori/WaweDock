"""Absolute-path entry point for source and pip Windows file associations."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mds.app import main

if __name__ == "__main__":
    raise SystemExit(main())
