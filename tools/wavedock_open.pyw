"""Explorer launcher; uses an absolute source path regardless of working directory."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mds.app import main

if __name__ == "__main__":
    raise SystemExit(main())
