"""Per-process cache for correction files (test, used by utils/corrections_cached.py).

Why a separate module: with run_uproot_job + futures_executor (coffea 0.7) the
processor, including the functions stored in data/corrections*.coffea, is
unpickled again for every chunk. Functions built by running
`python utils/corrections*.py` are pickled by value, so a cache stored next to
them would start empty on every chunk. A module like this one is pickled by
reference (re-imported on the worker), so its cache lives as long as the
worker process and every chunk after the first reuses the parsed files.

Requires `utils/` to be importable on the workers (true for run.py, which runs
from decaf/analysis, and for condor jobs, which ship the decaf tarball).
"""
import correctionlib

_CACHE = {}


def get(path):
    """correctionlib.CorrectionSet.from_file(path), parsed once per process."""
    key = ('cset', path)
    if key not in _CACHE:
        _CACHE[key] = correctionlib.CorrectionSet.from_file(path)
    return _CACHE[key]


def cached(key, factory):
    """Generic per-process cache: factory() is called only the first time `key` is seen."""
    if key not in _CACHE:
        _CACHE[key] = factory()
    return _CACHE[key]


def clear():
    _CACHE.clear()


def info():
    """Keys currently cached (for checks / debugging)."""
    return sorted(_CACHE, key=str)
