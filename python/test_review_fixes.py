import numpy as np
import pytest
from pybv import write_brainvision

from run_manifest import has_stage
from run_pipeline import run

FS, TR, NVOL = 1000, 1800, 15


def _write(tmp_path, drop=None):
    rng = np.random.default_rng(0)
    n = TR * (NVOL + 2)
    on = TR + TR * np.arange(NVOL)
    x = rng.normal(scale=3, size=(7, n))
    ev = np.array([[o, 1, 128] for o in on])
    if drop is not None:
        ev = np.delete(ev, drop, axis=0)
    ev = np.vstack([ev, [[on[3] + 100, 7, 1]]])                 # one extra, non-volume event
    names = ["Fz", "Cz", "Pz", "Oz", "CW1", "CW2", "ECG"]
    write_brainvision(data=x * 1e-6, sfreq=FS, ch_names=names, fname_base="t", folder_out=tmp_path,
                      events=ev, unit="V")
    return tmp_path / "t.vhdr", on


def test_missing_marker_stops(tmp_path):
    vhdr, _ = _write(tmp_path, drop=6)
    with pytest.raises(SystemExit, match="deviate from the TR"):
        run(vhdr, tmp_path / "out", marker="S  1")


def test_output_keeps_events_and_channel_types(tmp_path):
    import mne
    vhdr, on = _write(tmp_path)
    run(vhdr, tmp_path / "out", marker="S  1")
    r = mne.io.read_raw_fif(tmp_path / "out" / "cwl_corrected_raw.fif", preload=False, verbose="error")
    assert r.get_channel_types() == ["eeg"] * 4 + ["misc", "misc", "ecg"]
    t = {round(a["onset"], 3): a["description"] for a in r.annotations}
    # the extra event sits 100 samples after the 4th volume marker; output starts at the first marker
    assert any(abs(k - (on[3] + 100 - on[0]) / FS) < 0.01 and "S  7" in v for k, v in t.items()), t
    assert sum("S  1" in v for v in t.values()) == NVOL


def test_manifest_gates_stage(tmp_path):
    (tmp_path / "05_bcg_aas.set").write_text("old")
    assert has_stage(tmp_path, "05_bcg_aas")                     # no manifest -> legacy fallback
    (tmp_path / "run_manifest.json").write_text('{"stages": ["04_cwl_corrected"]}')
    assert not has_stage(tmp_path, "05_bcg_aas")
