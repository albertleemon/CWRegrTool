# 在本机运行（Python 版）

## 准备
1. 下载本仓库分支 `claude/awesome-clarke-obslcf`（`git clone` 后 `git checkout claude/awesome-clarke-obslcf`，或在 GitHub 页面下载 ZIP）。
2. Python 3.10+，安装依赖：
   ```
   pip install -r python/requirements.txt
   ```
3. 把 `.vhdr`、`.vmrk`、`.eeg` 三个文件放在同一个文件夹。
   - 你的 `.vhdr` 里写的是 `DataFile=RS1.eeg`。如果 `.eeg` 已经改名（例如 `sub-BOR5A_task-rest_run-01_eeg.eeg`），脚本会自动改用与 `.vhdr` 同名的 `.eeg` 和 `.vmrk`，不用手动改。

## 先用小片段试跑（约 1–2 分钟）
```
python python/crop_brainvision.py sub-BOR5A_task-rest_run-01_eeg.vhdr crop/ --volumes 40
python python/run_pipeline.py crop/sub-BOR5A_task-rest_run-01_eeg_crop.vhdr out_crop/ --exclude 65 66
python python/diagnose.py out_crop/
```

## 整条记录（建议内存 ≥ 16 GB，预计 5–15 分钟）
```
python python/run_pipeline.py sub-BOR5A_task-rest_run-01_eeg.vhdr out/ --exclude 65 66
python python/diagnose.py out/
```
没有头动参数时自动使用 moving average 模板；有 SPM 的 `rp_*.txt` 时加 `--rp rp_run01.txt`。

## 输出
| 文件 | 内容 |
|---|---|
| `out/cwl_corrected_raw.fif` | 校正后的数据（500 Hz，所有通道；可用 MNE 读取） |
| `out/cwl_corrected.npz` | 校正后数据、CWL 回归减去的伪影、回归前的 EEG |
| `out/pipeline_log.json` | 使用的参数、通道、volume 数，以及**削波通道占比** |
| `out/diagnostics.png`、`out/diagnostics.json` | 回归前后谱、10 s 示例、每通道去除的方差 |

**把 `pipeline_log.json`、`diagnostics.png`、`diagnostics.json` 发给我就可以**（都很小）。

## 需要注意
- 记录长度约 3.47 M 样本，最后一个 volume 之后仍有约 12 s，不会触发“最后一个 artifact 超出记录”。
- 输出只保留第一个到最后一个 volume 之间的部分，约 5 分 20 秒；校正范围外没有模板，仍含梯度伪影，所以被剪掉。
- 如果 `pipeline_log.json` 里 `clipped_channels_fraction` 不为空，说明有通道的梯度伪影超出放大器量程。这一版 Python 没有削波处理，请把这个字段一并发给我。

## 头动参数：用 fMRIPrep 的 confounds.tsv 代替 SPM 的 rp 文件
`--rp` 直接接受 fMRIPrep 的 `*_desc-confounds_timeseries.tsv`（读取 `trans_x/y/z`、`rot_x/y/z`，单位 mm 与弧度，与 SPM rp 相同）：
```
python python/run_pipeline.py "$IN" $OUT/out --exclude 65 66 \
   --rp sub-BOR5A_task-rest_run-01_desc-confounds_timeseries.tsv --threshold 0.15
```
MATLAB 版只读 rp 文本文件，先转换：`python python/fmriprep_to_rp.py confounds.tsv rp_run01.txt`，再设 `opts.rp_file`、`opts.threshold`。

对 sub-BOR5A 这个 run：tsv 有 356 行，与 `R128` 数量一致；相邻 volume 的最大平移位移只有 0.293 mm，所以论文的 0.3 mm 阈值下没有任何头动事件，会退化为 moving average。阈值 0.15 mm 约有 5 个事件，0.1 mm 有 18 个。
