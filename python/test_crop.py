import numpy as np
from crop_brainvision import crop
from run_pipeline import volume_onsets


def test_crop_roundtrip(tmp_path):
    import mne
    nch, fs, tr = 3, 5000, 9000
    n = 60000 + tr * 8
    rng = np.random.default_rng(0)
    x = rng.integers(-300, 300, size=(n, nch)).astype("<i2")
    x.tofile(tmp_path / "RS1.eeg")
    ch = "\n".join(f"Ch{i+1}=C{i},,0.5,µV" for i in range(nch))
    (tmp_path / "a.vhdr").write_text(
        f"Brain Vision Data Exchange Header File Version 1.0\n[Common Infos]\nCodepage=UTF-8\n"
        f"DataFile=RS1.eeg\nMarkerFile=a.vmrk\nDataFormat=BINARY\nDataOrientation=MULTIPLEXED\n"
        f"NumberOfChannels={nch}\nSamplingInterval=200\n[Binary Infos]\nBinaryFormat=INT_16\n"
        f"[Channel Infos]\n{ch}\n", encoding="utf-8")
    vols = [30001 + tr * i for i in range(8)]
    mk = "\n".join(f"Mk{i+2}=Response,R128,{p},1,0" for i, p in enumerate(vols))
    (tmp_path / "a.vmrk").write_text(
        "Brain Vision Data Exchange Marker File, Version 1.0\n[Common Infos]\nDataFile=RS1.eeg\n"
        f"[Marker Infos]\nMk1=New Segment,,1,1,0,20180429\n{mk}\n", encoding="utf-8")
    crop(tmp_path / "a.vhdr", tmp_path / "o", volumes=4, pre=1.0)
    raw = mne.io.read_raw_brainvision(tmp_path / "o" / "a_crop.vhdr", preload=True, verbose="error")
    on = volume_onsets(raw, "R128")
    assert len(on) == 4 and on[0] == 5000 and set(np.diff(on)) == {tr}
    # data identical to the source at the same position (resolution 0.5 uV, in volts)
    src = x[30000 - 5000 + 7, 0] * 0.5e-6
    assert abs(raw.get_data()[0, 7] - src) < 1e-12
