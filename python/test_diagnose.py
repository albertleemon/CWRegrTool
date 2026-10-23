import json
import shutil
import numpy as np
from test_pipeline import test_end_to_end
from diagnose import main


def test_rename_and_diagnose(tmp_path):
    test_end_to_end(tmp_path)
    # simulate a renamed recording whose header still names the old files
    d = tmp_path / "ren"; d.mkdir()
    for ext in ("vhdr", "vmrk", "eeg"):
        shutil.copy(tmp_path / f"t.{ext}", d / f"sub_x.{ext}")
    from run_pipeline import run
    run(d / "sub_x.vhdr", tmp_path / "out2", marker="S  1")
    main(tmp_path / "out2")
    assert (tmp_path / "out2" / "diagnostics.png").exists()
    s = json.loads((tmp_path / "out2" / "diagnostics.json").read_text())
    assert s["mean_var_removed_pct"] > 0
