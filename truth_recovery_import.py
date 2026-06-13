"""Importable shim for the hyphenated truth-recovery harness directory.

`truth-recovery/known_truth_validation.py` cannot be imported as a normal module
(the hyphen is not a valid identifier), so tests import `run_dist` from here.
"""
import importlib.util
from pathlib import Path

_path = Path(__file__).resolve().parent / "truth-recovery" / "known_truth_validation.py"
_spec = importlib.util.spec_from_file_location("known_truth_validation", _path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

run_dist = _mod.run_dist
DISTS = _mod.DISTS
