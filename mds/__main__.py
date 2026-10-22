"""Allows python -m mds, including frozen repair worker processes."""
import multiprocessing
import sys

if __name__ == "__main__":
    multiprocessing.freeze_support()
    from mds.app import main
    sys.exit(main())
