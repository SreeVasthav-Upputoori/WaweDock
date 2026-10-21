"""Optional repair command; normal Windows startup registers automatically."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mds.windows_integration import register

if __name__ == "__main__":
    register()
    print("WaweDock Open with registration is up to date; defaults are unchanged.")
