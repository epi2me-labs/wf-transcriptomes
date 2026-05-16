#!/usr/bin/env python
"""Import path setup for wf-container workflow_glue tests."""

from pathlib import Path
import sys


# Add /host/bin so `import workflow_glue` resolves in CI containers.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
