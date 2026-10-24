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

对 sub-BOR5A 这个 run：tsv 有 356 行，与 `R128` 数量一致；相邻 volume 的最大平移位移只有 0.293 mm，所以论文的 0.3 mm 阈值下没有任何头动事件，会退化为 moving average。阈值 0.15 mm 有 6 个事件，0.1 mm 有 18 个。


## 可选：梯度残余的 OBS 去除（Niazy 2005，FASTR）
梯度模板减法之后，若高频仍有层频率谐波残峰（例如 sub-BOR5A 的 84.5、101、117.75 Hz），可加残差 PCA：
对每个通道，把模板减法后的残差高通到 70 Hz，按 volume 做 PCA，去掉前 N 个主成分。只拟合 70 Hz 以上的残差，40 Hz 以下的脑电不受影响。
**这一步不在 van der Meer et al. (2016) 的方法里，默认关闭；使用时请在方法部分注明 Niazy et al. (2005)。**
```
CWL_OBS=4 ./matlab/run_headless.sh "$IN" $OUT/ma_obs4 65,66          # MATLAB
python python/run_pipeline.py "$IN" $OUT/py_obs4 --exclude 65 66 --obs-npc 4   # Python
```


## 完整流程一条命令（推荐）
```
export MATLAB_BIN=/media/bianca/slurm/Matlab/R2024b/bin/matlab
conda activate cwl            # 或你装了依赖的 Python 环境
./run_subject.sh RAW.vhdr OUT_DIR 65,66            # 默认含 OBS（CWL_OBS=4）
```
| 阶段 | 文件 | 方法 |
|---|---|---|
| 02 梯度模板减法（5000 Hz，全部通道） | `02_gradient_corrected.set` | Allen 2000；有头动参数时 Moosmann 2009 |
| 02b OBS 残差去除（5000 Hz，全部通道） | `02b_obs_corrected.set` | Niazy 2005（FASTR 的 OBS），70 Hz 以上 |
| 03 带通 0.33–125 Hz + 500 Hz | `03_filtered.set` | van der Meer 2016 |
| 04 CWL 回归（EEG 通道） | `04_cwl_corrected.set` ← 最终结果 | van der Meer 2016 |

诊断：`stage_spectra.png`（原始、02、02b、04、随机相位 CWL 对照的平均谱）、`diagnostics_py.png`、`diagnostics.json`、`residual_peaks.png/json`。
磁盘：02 和 02b 各约 1 GB（5000 Hz），确认结果后可删除。


## 可选参数（运行 `run_subject.sh` 前用环境变量设置）
| 参数 | 默认 | 作用 |
|---|---|---|
| `CWL_OBS` | 4 | OBS 主成分个数；0 = 跳过 OBS |
| `CWL_OBS_HP` | 70 | OBS 只拟合该频率以上的残差；60 可同时处理 67 Hz 谐波 |
| `CWL_LOWPASS` | 关 | 回归前只对碳线环信号做低通（8 阶，零相位），例如 90：回归不再把环路自身的高频噪声带进 EEG，但该频率以上的伪影也不再被去除 |
| `CWL_RP` | 无 | SPM rp 文件或 fMRIPrep confounds.tsv，用于头动信息模板 |
| `CWL_THRESHOLD` | 0.3 | 头动阈值（mm） |

例：`CWL_OBS_HP=60 CWL_LOWPASS=90 ./run_subject.sh RAW.vhdr OUT 65,66`
每次运行的设置写在 `OUT/settings.txt`，也记录在 `04_cwl_corrected.set` 的 `EEG.etc.cwl_pipeline` 里。
Python 版对应：`--obs-npc`、`--obs-hp`、`--cwl-lowpass`。
