"""Which stage files of a MATLAB output folder belong to the last run (run_manifest.json)."""
import json
from pathlib import Path


def has_stage(out_dir, stage):
    """True only if `stage` (e.g. '05_bcg_aas') was written by the run that wrote the manifest.
    A folder without a manifest (older runs) falls back to 'file exists', with a warning."""
    d = Path(out_dir)
    f = d / (stage + ".set")
    mf = d / "run_manifest.json"
    if mf.exists():
        return stage in json.loads(mf.read_text()).get("stages", []) and f.exists()
    if f.exists():
        print(f"WARNING: no run_manifest.json in {d}; {f.name} may be left over from an earlier run")
    return f.exists()
