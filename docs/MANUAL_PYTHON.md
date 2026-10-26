# Python 版使用手册（v0.3.0）

适用：`python/` 文件夹。MATLAB 版流程见 `docs/MANUAL_MATLAB.md`。

先说清楚两个版本的分工，避免混淆：

| | MATLAB 版 | Python 版 |
|---|---|---|
| 校正流程 | `cwl_run_pipeline`，每个阶段存 `.set` | `python/run_pipeline.py`，只存最终结果 |
| CWL 回归 | **原作者代码** | 按论文和原 MATLAB 代码移植（`cwregression.py`） |
| 真实数据 | 已在 sub-BOR5A、sub-kvhDR 上运行 | **尚未在真实数据上运行**，合成数据端到端测试通过 |
| 评估与报告脚本 | – | 评估/报告脚本（`evaluate_*`、`denoise_report.py` 等）**读的是 MATLAB 输出的 `.set`**，不读 Python 流程的输出 |

所以：想要去噪报告、头皮图和与磁体外的比较，请用 MATLAB 版流程做校正，再用本手册第 6、7 节的脚本做评估。Python 流程适合：不方便用 MATLAB、做交叉验证，或二次开发。

---

## 1. 环境

```bash
python -m venv cwl-env && source cwl-env/bin/activate     # 或 conda create -n cwl python=3.11
pip install -r python/requirements.txt
```
依赖：`numpy scipy mne>=1.6 matplotlib pymatreader`；测试额外需要 `pytest pybv`。

内存：72 通道 × 5 分钟 × 5000 Hz，双精度约 0.9 GB，运算时有多份拷贝，建议 ≥ 16 GB；整条记录预计 5–15 分钟。

## 2. 数据准备

- `.vhdr`、`.vmrk`、`.eeg` 放在同一文件夹。若 `.vhdr` 里 `DataFile=` 写的是旧文件名，脚本会自动改用与 `.vhdr` 同名的 `.eeg` 和 `.vmrk`。
- 通道：CWL = `CW1`–`CW6`，ECG = `ECG`，未使用的通道（如 65/66）用 `--exclude` 排除。
- 头动参数：SPM `rp_*.txt` 或 fMRIPrep `*_desc-confounds_timeseries.tsv`（读取 `trans_x/y/z`、`rot_x/y/z`，mm 与弧度）。`--rp` 两种都接受。
- 体积标记默认 `R128`（BrainVision 的 `Response,R128` 被 MNE 读为 `Response/R128`，脚本取最后一段匹配）。
- 记录在最后一个 volume 之后至少还有 1 个 TR。

## 3. 先用小片段试跑（1–2 分钟）

```bash
python python/crop_brainvision.py sub-XXX_task-rest_run-01_eeg.vhdr crop/ --volumes 40
python python/run_pipeline.py crop/sub-XXX_task-rest_run-01_eeg_crop.vhdr out_crop/ --exclude 65 66
python python/diagnose.py out_crop/
```
`crop_brainvision.py` 只支持 INT_16 + MULTIPLEXED，输出一个可以分享的小文件（30 个 volume 约 40 MB）。

## 4. 完整运行：`run_pipeline.py`

```bash
python python/run_pipeline.py RAW.vhdr OUT/ --exclude 65 66 \
       --rp sub-XXX_task-rest_run-01_desc-confounds_timeseries.tsv --threshold 0.15 \
       --obs-npc 4 --obs-hp 60 --artifact-lowpass 90
python python/diagnose.py OUT/
```

### 4.0 速查：必需参数与可选参数

**必需**：`run_pipeline.py` 只有两个位置参数 ——

| 参数 | 说明 |
|---|---|
| `vhdr` | 原始 BrainVision `.vhdr` 路径 |
| `out_dir` | 输出目录（不存在会自动创建；建议放在 BIDS 之外） |

**强烈建议设置**（技术上可选）
- `--exclude 65 66`：排除未使用的通道，否则它们会被当作 EEG 目标通道。
- `--rp ...`：有头动文件时给，否则用 moving average 模板（`--threshold` 只在给了 `--rp` 时起作用）。

**其余全部可选**：

| 用途 | 参数 | 默认 |
|---|---|---|
| 数据识别 | `--marker`、`--cwl-regex`、`--ecg` | `R128`、`^CW`、`ECG` |
| 梯度校正 | `--rp`、`--threshold`、`--templates`、`--span PRE POST` | 无、0.3、25、`0 TR-1` |
| 回归 | `--target-fs` | 500 |
| 可选改进（不是论文原方法，默认关闭） | `--obs-npc`、`--obs-hp`、`--artifact-lowpass`、`--bcg-aas`、`--bcg-hr` | 0、70、关、关、无 |
| 安全检查 | `--no-marker-check`（关闭检查） | 检查开启 |

**用法示例**
```bash
# 1) 最小用法
python python/run_pipeline.py RAW.vhdr OUT/

# 2) 论文原方法 + 排除未用通道
python python/run_pipeline.py RAW.vhdr OUT/ --exclude 65 66

# 3) 加头动信息
python python/run_pipeline.py RAW.vhdr OUT/ --exclude 65 66 --rp confounds.tsv --threshold 0.15

# 4) 推荐改进组合
python python/run_pipeline.py RAW.vhdr OUT/ --exclude 65 66 --obs-npc 4 --obs-hp 60 --artifact-lowpass 90

# 5) 再加心跳锁定 AAS
python python/run_pipeline.py RAW.vhdr OUT/ --exclude 65 66 --obs-npc 4 --obs-hp 60 --artifact-lowpass 90 --bcg-aas --bcg-hr 64
```

### 4.1 命令行参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `vhdr` `out_dir` | – | 输入头文件、输出目录 |
| `--rp` | 无 | 头动文件；省略则用 moving average 模板 |
| `--threshold` | 0.3 | 头动阈值（mm）。示例数据最大帧间位移 0.293 mm，0.3 下无事件；0.15 约 6 个，0.1 约 18 个 |
| `--templates` | 25 | 每个 volume 的模板数 |
| `--marker` | `R128` | 体积标记 |
| `--cwl-regex` | `^CW` | CWL 通道正则 |
| `--ecg` | `ECG` | ECG 通道名 |
| `--exclude` | 无 | 要排除的通道名，如 `65 66` |
| `--span PRE POST` | `0 TR-1` | 伪影范围（样本） |
| `--target-fs` | 500 | 回归采样率 |
| `--obs-npc` | 0 | OBS 主成分数（0 = 关；推荐 4；不在原论文中） |
| `--obs-hp` | 70 | OBS 高通（Hz），60 可同时处理 67 Hz 谐波 |
| `--artifact-lowpass` | 无 | 对 CWL 估计的伪影做低通（Hz），推荐 90 |
| `--bcg-aas` | 关 | 回归后再做心跳锁定 AAS（Allen 1998） |
| `--bcg-hr` | 无 | 预期心率（bpm），帮助 R 波检测 |
| `--no-marker-check` | 关 | 默认遇到漏标、重复、间隔偏离 TR 超过 1% 会报错；此项关闭检查 |

### 4.2 流程

1. 读取（MNE），电压换算成 µV。
2. 体积标记检查。
3. 梯度模板减法（全部通道，含 CWL）；可选 OBS。
4. 剪掉校正范围之外的首尾。
5. 0.33–125 Hz 带通（全部通道）→ 降采样到 500 Hz。
6. CWL 回归（EEG 通道）；可选估计伪影低通。
7. 可选 BCG AAS。

### 4.3 输出

| 文件 | 内容 |
|---|---|
| `cwl_corrected_raw.fif` | 校正后数据（500 Hz，全部通道，MNE 可读）。保留标记事件（时间重新对齐到输出起点）、通道类型（EEG / ECG / 其余 misc）和电极位置 |
| `cwl_corrected.npz` | `data`、回归前的 EEG（`pre_cwl`）、减去的伪影（`artifact`）、`ch_names`、`fs`、`eeg_idx`、`cwl_idx` |
| `pipeline_log.json` | 参数、通道、volume 数、TR、头动事件数、**削波通道占比**（`clipped_channels_fraction`，不为空说明梯度伪影超出量程；Python 版没有削波处理）、输出样本数 |
| `diagnostics.png/.json` | `diagnose.py` 生成：回归前后谱、10 s 示例、每通道去除方差 |

### 4.4 与 MATLAB 版的行为差异

- CWL 回归窗口在记录首尾保持满权重，边缘也被校正；MATLAB 原工具首尾不校正，所以输出会被裁掉。两者在满权重内部区域一致（相对 RMS 差约 1e-6）。
- 降采样用 `scipy.signal.resample_poly`，MATLAB 版用 EEGLAB 风格 Kaiser 窗滤波；在无混叠频段内差异很小（合成测试最大差约 5e-4）。
- 没有移植 Bergen 的削波/NaN 处理和基线校正。
- 不生成逐阶段的 02/02b/03 文件。

## 5. 作为库调用

```python
import sys; sys.path.insert(0, "python")
from run_pipeline import run
out, names = run("RAW.vhdr", "OUT", rp="confounds.tsv", threshold=0.15,
                 exclude=("65", "66"), obs_npc=4, obs_hp=60.0, artifact_lowpass=90.0)

from cwregression import cw_regression, bandpass          # 单独使用 CWL 回归
corrected, artifact = cw_regression(eeg, cwl, fs, window_s=4.0, delay_s=0.021,
                                    filter_band=None, return_artifact=True, artifact_lowpass=90)
```
模块：`gradient_correction.py`（`rp_weighting_matrix`、`moving_average_matrix`、`gradient_correct`、`obs_residual`）、`cwregression.py`、`bcg_aas.py`（`r_peaks`、`bcg_aas`）、`metrics.py`（论文式 4、5、6，EO/EC 频带）。

## 6. 评估与报告脚本（对 MATLAB 输出运行）

输入是 MATLAB 版的输出目录 `OUT_DIR`（里面有 `03_filtered.set`、`04_cwl_corrected.set`，可选 `02`、`02b`、`05`，以及 `run_manifest.json`）。只有清单里列出的阶段才会被采用；没有清单的旧目录会警告。

最省事：`./evaluate_subject.sh OUT_DIR OUTSIDE.vhdr`（在仓库根目录，一次运行下面 6.1–6.4 全部，输出同时写入 `evaluation.txt`）。

### 6.0 速查：评估与报告脚本的必需与可选参数

| 脚本 | 必需 | 可选 |
|---|---|---|
| `denoise_report.py` | `OUT_DIR` | `--outside`（要有相对磁体外的指标就必须给）、`--vhdr`、`--hr`、`--ecg`、`--marker`、`--slice-freq`、`--crop` |
| `evaluate_outside.py` | `OUT_DIR`、`OUTSIDE.vhdr` | `--crop START END` |
| `compare_inside_outside.py` | `INSIDE.set`、`OUTSIDE.vhdr`、`OUT.png` | `--channels`、`--start`、`--seconds`、`--also` |
| `check_bcg_residual.py` | `OUT_DIR` | `--outside`、`--ecg`、`--hr` |
| `check_removed_signal.py` | `OUT_DIR` | `--outside`、`--channels`、`--k` |
| `plot_stage_spectra.py` | `OUT_DIR`（且里面要有 `diagnostics_data.mat`） | `--vhdr`、`--title` |
| `analyze_residual_peaks.py` | `OUT_DIR` | `--confounds`、`--bold`、`--fmin` |
| `check_volume_alignment.py` | `RAW.vhdr`、`OUT_DIR` | `--channels`、`--marker`、`--k` |
| `plot_fig3.py` | `OUT.png`、至少一个 `--data`（或 `--mat`） | `--outside`、`--title` |
| `fmriprep_to_rp.py` | `confounds.tsv`、`rp.txt` | – |
| `read_siemens_puls.py` | `sub.puls` | `--start`、`--seconds`、`--bin` |

最常用的组合：
```bash
# 全部评估（仓库根目录）：两个参数都必需
./evaluate_subject.sh OUT_DIR OUTSIDE_run-00.vhdr
#   可选环境变量：CWL_RAW_VHDR=RAW.vhdr（加 raw 阶段）、CWL_BCG_HR=64、PYTHON=...
# 只要去噪报告
python python/denoise_report.py OUT_DIR --outside OUTSIDE_run-00.vhdr
# 去噪报告，完整可选项
python python/denoise_report.py OUT_DIR --outside OUTSIDE.vhdr --vhdr RAW.vhdr --hr 64 --slice-freq 16.84 --crop 10 300
```

### 6.1 `denoise_report.py`：去噪报告

```bash
python python/denoise_report.py OUT_DIR --outside OUTSIDE.vhdr [--vhdr RAW.vhdr] [--hr 64] \
       [--ecg ECG] [--marker R128] [--slice-freq 16.84] [--crop START END]
```
| 参数 | 说明 |
|---|---|
| `--outside` | 磁体外参考（BrainVision）；没有它就没有"相对磁体外"的指标 |
| `--vhdr` | 本次运行的原始记录，加入未校正的 raw 阶段 |
| `--hr` | 预期心率（bpm） |
| `--marker` | raw 阶段用于对齐的体积标记 |
| `--slice-freq` | 层频率（Hz），用于"层频率谐波残峰"指标，默认 16.84 |
| `--crop` | 只使用磁体外记录的这一段（秒） |

输出：`denoise_spectra.png`（每一步的平均谱 + 相对磁体外的 dB 差）、`denoise_indicators.png/.json/.csv`（指标，含义见 `docs/MANUAL_MATLAB.md` 第 10 节）、`denoise_topomaps.png`（每一步的头皮分布图：1–30 Hz RMS、心跳锁定 RMS、相对磁体外的 dB）。头皮图需要至少 16 个通道有标准 10-05 位置，否则跳过。所有阶段在相同样本和相同 EEG 通道上比较，先经过与流程相同的 0.33–125 Hz 带通和 500 Hz 降采样。

### 6.2 其他评估脚本

| 脚本 | 用法 | 作用 |
|---|---|---|
| `evaluate_outside.py` | `evaluate_outside.py OUT_DIR OUTSIDE.vhdr [--crop 10 300]` | 论文式 4 谱距离（全频带和 1–30 Hz）、频段 dB 差、alpha 地形；并检查参考是否确实在磁体外（CWL 1–20 Hz 功率比） |
| `compare_inside_outside.py` | `compare_inside_outside.py INSIDE.set OUTSIDE.vhdr OUT.png [--channels ...] [--start 60] [--seconds 10] [--also 03_filtered.set]` | 示例轨迹、平均谱、每通道 1–30 Hz RMS |
| `check_bcg_residual.py` | `check_bcg_residual.py OUT_DIR [--outside OUTSIDE.vhdr] [--ecg ECG] [--hr BPM]` | 心跳锁定平均（−0.2…+0.8 s）和残差 RMS；写 `bcg_residual.png/.json`、`bcg_topomap.png` |
| `check_removed_signal.py` | `check_removed_signal.py OUT_DIR [--outside OUTSIDE.vhdr] [--channels ...] [--k 4]` | CWL 去除成分与邻近电极的相干性（是否误去脑电）；灵敏度有限，隐藏在大伪影里的小损失检测不到 |
| `plot_stage_spectra.py` | `plot_stage_spectra.py OUT_DIR [--vhdr RAW.vhdr]` | 各阶段平均谱 + 随机相位对照（需要 `cwl_diagnose` 生成的 `diagnostics_data.mat`） |
| `plot_diagnostics_mat.py` | `plot_diagnostics_mat.py OUT_DIR` | 由 `diagnostics_data.mat` 画 `diagnostics_py.png` |
| `plot_fig3.py` | `plot_fig3.py fig3.png --data "列/行=路径" ... --outside OUTSIDE` | 论文图 3 风格的方法对比图 |

### 6.3 排查类脚本

| 脚本 | 用法 | 作用 |
|---|---|---|
| `analyze_residual_peaks.py` | `analyze_residual_peaks.py OUT_DIR [--confounds TSV] [--bold BOLD.json] [--fmin 40]` | 残峰是否与体积锁定、是否在 1/TR 网格上，每 volume 和每通道分布，与头动相关 |
| `check_volume_alignment.py` | `check_volume_alignment.py RAW.vhdr OUT_DIR [--channels ...] [--marker R128] [--k 25]` | 原始数据上逐 volume 的时间错位和削波检查 |
| `fmriprep_to_rp.py` | `fmriprep_to_rp.py confounds.tsv rp_run01.txt` | 把 fMRIPrep tsv 转成 SPM 样式 rp 文件，并打印帧间位移统计 |
| `read_siemens_puls.py` | `read_siemens_puls.py sub.puls [--start HH:MM:SS.sss] [--seconds N] [--bin 300]` | 读 Siemens PMU 脉搏日志，得到心率，可作为 `--hr` |

## 7. 测试

```bash
python -m pytest python -q
```
覆盖：端到端合成数据、梯度校正与指标、CWL 回归、裁剪、诊断；以及新增的回归测试——漏标记会报错、输出保留事件和通道类型、运行清单会挡住旧的阶段文件。

与 MATLAB 版的数值对照：`python python/make_crosscheck_data.py crosscheck_data.mat`，再在 MATLAB 里运行 `cwl_crosscheck`（见 MATLAB 手册第 11 节）。

## 8. 常见问题

| 现象 | 处理 |
|---|---|
| `volume interval(s) deviate from the TR` | 标记有漏标、重复或多余；核对 `.vmrk`，确认无误用 `--no-marker-check` |
| `found N 'R128' markers` | 标记名不对，用 `--marker` 指定（MNE 里形如 `Stimulus/S  1` 时给 `S  1`） |
| `last artifact exceeds the recording` | 最后一个 volume 之后不足 1 个 TR，去掉最后一个标记或裁剪 |
| `no CWL channels match` | 调整 `--cwl-regex` |
| 评估时警告没有 `run_manifest.json` | 旧目录，里面的 02b/05 可能是旧结果；用独立目录重跑 MATLAB 流程 |
| 画图时 `FT_Render_Glyph ... raster overflow` | 旧版 `denoise_report.py` 的已知问题，已修复；更新到最新分支 |
| 内存不足 | 先用 `crop_brainvision.py` 试小片段；整条记录需要 ≥ 16 GB |
| `pipeline_log.json` 中 `clipped_channels_fraction` 不为空 | 有通道的梯度伪影超出放大器量程，Python 版没有削波处理；请改用 MATLAB 版或先检查数据 |

## 9. 局限

- Python 流程尚未在真实数据上运行，结果请与 MATLAB 版交叉验证。
- OBS、估计伪影低通、BCG AAS 都不是论文原方法，使用时请在方法部分注明（Niazy 2005；Allen 1998）。
- 90 Hz 以上的伪影（氦泵开启时）不会被低通选项去除。
- 评估脚本的磁体外参考不是纯净真值（50 Hz 工频峰、<1 Hz 漂移、状态不同）。
