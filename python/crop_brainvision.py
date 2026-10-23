"""Cut a short excerpt out of a large BrainVision recording (INT_16, multiplexed).

The excerpt starts `--pre` seconds before the first volume marker and ends one TR
after the N-th volume marker, so it is a valid input for run_pipeline.py but small
enough to share (30 volumes at 5000 Hz, 72 channels ~ 40 MB).

    python crop_brainvision.py in.vhdr out_dir --volumes 30
"""
import argparse
import re
from pathlib import Path

import numpy as np


def crop(vhdr, out_dir, volumes=30, pre=1.0, marker="R128"):
    vhdr = Path(vhdr)
    h = vhdr.read_text(encoding="utf-8", errors="replace")
    nch = int(re.search(r"NumberOfChannels=(\d+)", h).group(1))
    fs = 1e6 / float(re.search(r"SamplingInterval=([\d.]+)", h).group(1))
    assert "INT_16" in h and "MULTIPLEXED" in h, "only INT_16 multiplexed supported"
    data_name = re.search(r"DataFile=(.+)", h).group(1).strip()
    mrk_name = re.search(r"MarkerFile=(.+)", h).group(1).strip()
    # the header may point to an old file name (e.g. RS1.eeg): fall back to siblings
    eeg = vhdr.parent / data_name
    if not eeg.exists():
        eeg = vhdr.with_suffix(".eeg")
    vmrk = vhdr.parent / mrk_name
    if not vmrk.exists():
        vmrk = vhdr.with_suffix(".vmrk")
    m = vmrk.read_text(encoding="utf-8", errors="replace")

    vols = [int(x) for x in re.findall(rf"^Mk\d+=Response,{marker},(\d+),", m, re.M)]
    tr = int(np.median(np.diff(vols)))
    a = max(vols[0] - 1 - int(pre * fs), 0)                         # 0-based start
    last = vols[min(volumes, len(vols)) - 1]
    b = last - 1 + tr                                               # 0-based exclusive end

    mm = np.memmap(eeg, dtype="<i2", mode="r").reshape(-1, nch)
    seg = np.array(mm[a:b])
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = vhdr.stem + "_crop"
    seg.tofile(out_dir / f"{stem}.eeg")

    lines, k = [], 0
    for ln in m.splitlines():
        mt = re.match(r"^Mk\d+=([^,]*),([^,]*),(\d+),(.*)$", ln)
        if mt:
            pos = int(mt.group(3))
            if mt.group(1) == "New Segment":
                pos_new = 1
            elif a < pos <= b and not (mt.group(2) == marker and pos > last):  # 1-based
                pos_new = pos - a
            else:
                continue
            k += 1
            ln = f"Mk{k}={mt.group(1)},{mt.group(2)},{pos_new},{mt.group(4)}"
        elif ln.startswith("DataFile="):
            ln = f"DataFile={stem}.eeg"
        lines.append(ln)
    (out_dir / f"{stem}.vmrk").write_text("\n".join(lines) + "\n", encoding="utf-8")
    h2 = re.sub(r"DataFile=.+", f"DataFile={stem}.eeg", h)
    h2 = re.sub(r"MarkerFile=.+", f"MarkerFile={stem}.vmrk", h2)
    (out_dir / f"{stem}.vhdr").write_text(h2, encoding="utf-8")
    print(f"wrote {out_dir / stem}.*  {seg.shape[0]} samples, {seg.nbytes / 1e6:.1f} MB")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vhdr"); p.add_argument("out_dir")
    p.add_argument("--volumes", type=int, default=30)
    p.add_argument("--pre", type=float, default=1.0, help="seconds kept before the first volume")
    p.add_argument("--marker", default="R128")
    a = p.parse_args()
    crop(a.vhdr, a.out_dir, a.volumes, a.pre, a.marker)
