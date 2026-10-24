"""Heart rate from a Siemens PMU pulse log (*.puls, finger pulse oximeter).

File layout: a header (first values, then "5002 LOGVERSION... 6002"), the pulse wave
samples (400 Hz on VE11 / Prisma, 50 Hz on older VB systems; the rate is derived from the
log duration), "5000" inserted where the scanner detected a pulse beat, "5003", and a
footer with the scanner's own average ("PULS Freq Per: <bpm> <period ms>") and the log
start/stop times in ms since midnight (MDH clock = the clock of AcquisitionTime in the
BIDS bold.json).

    python read_siemens_puls.py sub-X.puls                       # whole log, 5-min bins
    python read_siemens_puls.py sub-X.puls --start 15:55:17.38 --seconds 640.8
"""
import argparse
import re

import numpy as np


def read_puls(path):
    txt = open(path).read()
    body = txt.split("5003")[0].split()
    i = body.index("6002") + 1 if "6002" in body else 4
    samples, trig = [], []
    for tkn in body[i:]:
        if not tkn.lstrip("-").isdigit():
            continue
        v = int(tkn)
        if v == 5000:
            trig.append(len(samples))
        elif v not in (5002, 6002):
            samples.append(v)
    times = {k: int(v) for k, v in re.findall(r"(Log\w+Time):\s+(\d+)", txt)}
    footer = re.search(r"PULS Freq Per:\s+(\d+)\s+(\d+)", txt)
    start, stop = times["LogStartMDHTime"], times["LogStopMDHTime"]
    fs = len(samples) / ((stop - start) / 1000)
    return np.array(samples), np.array(trig) / fs, fs, start, stop, footer


def clock(ms):
    return "%02d:%02d:%06.3f" % (ms // 3600000, (ms // 60000) % 60, (ms % 60000) / 1000)


def to_ms(hms):
    h, m, s = hms.split(":")
    return int((int(h) * 3600 + int(m) * 60 + float(s)) * 1000)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("puls")
    ap.add_argument("--start", help="window start, HH:MM:SS(.sss) scanner clock (bold.json AcquisitionTime)")
    ap.add_argument("--seconds", type=float, help="window length, e.g. n_volumes * TR")
    ap.add_argument("--bin", type=float, default=300, help="bin length in s for the whole-log summary")
    a = ap.parse_args()
    x, tb, fs, start, stop, footer = read_puls(a.puls)
    ibi = np.diff(tb)
    ok = (ibi > 0.3) & (ibi < 2.0)
    hr = 60 / ibi
    tm = tb[1:]
    print(f"log {clock(start)} - {clock(stop)} ({(stop - start) / 1000:.0f} s), pulse sampled at {fs:.1f} Hz, "
          f"{len(tb)} detected beats")
    if footer:
        print(f"scanner footer: PULS {footer.group(1)} bpm, period {footer.group(2)} ms")
    print(f"whole log: HR median {np.median(hr[ok]):.1f} bpm (IQR {np.percentile(hr[ok], 25):.1f}-"
          f"{np.percentile(hr[ok], 75):.1f})")
    if a.start:
        t0 = (to_ms(a.start) - start) / 1000
        t1 = t0 + (a.seconds or 600)
        sel = ok & (tm >= t0) & (tm < t1)
        print(f"window {a.start} + {t1 - t0:.1f} s: {sel.sum()} beats, HR median {np.median(hr[sel]):.1f} bpm "
              f"(IQR {np.percentile(hr[sel], 25):.1f}-{np.percentile(hr[sel], 75):.1f}, "
              f"min {hr[sel].min():.0f}, max {hr[sel].max():.0f})")
    else:
        for b0 in np.arange(0, (stop - start) / 1000, a.bin):
            sel = ok & (tm >= b0) & (tm < b0 + a.bin)
            if sel.sum() > 5:
                print(f"  {clock(start + int(b0 * 1000))[:8]}  HR {np.median(hr[sel]):5.1f} bpm  ({sel.sum()} beats)")


if __name__ == "__main__":
    main()
