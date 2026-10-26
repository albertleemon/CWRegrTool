# MATLAB 版使用手册（v0.3.0）

适用：`matlab/` 文件夹（独立版，不需要 EEGLAB）。评估与报告部分是 Python 脚本，见第 9 节；Python 版流程见 `docs/MANUAL_PYTHON.md`。

---

## 1. 它做什么

对同步 EEG–fMRI 记录做伪影校正（van der Meer et al., NeuroImage 2016 的碳线环 CWL 回归），并在此基础上提供可选的改进步骤和评估：

```
BrainVision 原始数据（5000 Hz，EEG + ECG + CWL）
  ├─ 02  梯度伪影模板减法（全部通道，含 CWL；有头动参数时用 rp-informed 模板）
  ├─ 02b OBS 残差 PCA（可选，Niazy 2005，不在原论文中）
  ├─ 03  0.33–125 Hz 带通 + 降采样到 500 Hz（全部通道同一滤波器）
  ├─ 04  CWL 回归（只作用于 EEG 通道；可选对"估计出的伪影"做低通）
  └─ 05  心跳锁定 BCG AAS（可选，Allen 1998，默认关闭）
```

每个阶段保存成 EEGLAB 格式（`.set` + `.fdt`），最终结果是 `04_cwl_corrected.set`（开了 BCG AAS 时是 `05_bcg_aas.set`）。

## 2. 环境要求

| 项目 | 要求 |
|---|---|
| MATLAB | 已在 R2024b 测试 |
| 工具箱 | **Signal Processing Toolbox**（`butter`、`filtfilt`、`resample`、`pwelch`） |
| 不需要 | EEGLAB、bva-io、firfilt 插件、Statistics Toolbox、网络 |
| 内存 | 72 通道 × 5 分钟 × 5000 Hz，建议 ≥ 16 GB |
| 磁盘 | 02、02b 各约 1 GB（5000 Hz），确认结果后可删 |
| 评估/报告（可选） | Python 3.10+：`numpy scipy mne pymatreader matplotlib`（见 `python/requirements.txt`） |

## 3. 安装

1. 取得本仓库分支 `claude/awesome-clarke-obslcf`（`git clone` 后 `git checkout claude/awesome-clarke-obslcf`，或下载 ZIP）。
2. `matlab/` 文件夹自带全部依赖，可以单独复制到任何位置。
3. 在 MATLAB 中：
   ```matlab
   cd /path/to/CWRegrTool/matlab
   setup_cwl_pipeline      % 只设置路径，并检查 Signal Processing Toolbox
   ```
   如果路径里有旧的 EEGLAB 或 CWRegrTool，不会冲突：本文件夹会放在路径最前面。
4. 若 `setup_cwl_pipeline` 报 Signal Processing Toolbox 函数失败，多半是 MATLAB 路径被破坏，运行 `fix_matlab_path` 或 `debug_cwl_env` 查看。

## 4. 运行前检查（强烈建议）

1. `R128` 标记数 = 实际 EPI volume 数（含 dummy），相邻间隔 = TR（本流程会自动检查，见第 8 节）。
2. 通道：CWL = `CW1`–`CW6`，ECG = `ECG`，未使用的通道（如 65/66）放进 `exclude`。
3. 若用头动参数：fMRIPrep 的 `*_desc-confounds_timeseries.tsv`（读取 `trans_x/y/z`、`rot_x/y/z`）或 SPM `rp_*.txt`，行数应等于 volume 数。
4. 记录在最后一个 volume 之后至少还有 1 个 TR 的数据。
5. BrainVision 头文件里的通道单位（µV / mV / V）会被自动换算成 µV。若 `DataFile=` 写的是旧文件名，会自动改用与 `.vhdr` 同名的 `.eeg` / `.vmrk`。
6. **每组参数用独立的输出目录。**

## 5. 三种运行方式

### 5.1 在 MATLAB 里直接调用

```matlab
cd /path/to/CWRegrTool/matlab
setup_cwl_pipeline
opts = struct('exclude', {{'65','66'}}, 'obs_npc', 4, 'obs_hp', 60, ...
              'cwl_lowpass', 90, 'out_dir', '/path/to/out');
EEG = cwl_run_pipeline('/path/to/sub-XXX_task-rest_run-01_eeg.vhdr', opts);
cwl_diagnose(opts.out_dir);      % 回归前后对比 + 随机相位对照
```
示例脚本：`example_BOR5A.m`、`example_NB.m`。

### 5.2 命令行（不开界面）

```bash
export MATLAB_BIN=/usr/local/MATLAB/R2024b/bin/matlab     # 按你的安装位置
CWL_OBS=4 CWL_OBS_HP=60 CWL_LOWPASS=90 ./matlab/run_headless.sh RAW.vhdr OUT_DIR 65,66
```
日志写入 `OUT_DIR/matlab_run.log`，出错时退出码非 0。

### 5.3 一条命令跑完整流程（推荐）

```bash
export MATLAB_BIN=/path/to/matlab
conda activate cwl                      # 装了 Python 依赖的环境
CWL_OBS=4 CWL_OBS_HP=60 CWL_LOWPASS=90 ./run_subject.sh RAW.vhdr OUT_DIR 65,66
./evaluate_subject.sh OUT_DIR OUTSIDE_run-00.vhdr    # 与磁体外记录比较，见第 9 节
```
`run_subject.sh` 依次做：校正（MATLAB）→ `cwl_diagnose`（含随机相位对照）→ 各阶段平均谱 → 图和残峰分析。

## 6. 参数

### 6.0 速查：必需参数与可选参数

**必需（不给就无法运行）**

| 入口 | 必需项 | 说明 |
|---|---|---|
| `cwl_run_pipeline(vhdr, opts)` | `vhdr` | 原始 BrainVision `.vhdr` 的完整路径；`opts` 可以省略（全用默认值） |
| `run_headless.sh RAW.vhdr OUT_DIR [excl]` | `RAW.vhdr`、`OUT_DIR` | 两个位置参数；第 3 个（排除通道，默认 `65,66`）可省略 |
| `run_subject.sh RAW.vhdr OUT_DIR [excl]` | 同上 | 同上 |
| `evaluate_subject.sh OUT_DIR OUTSIDE.vhdr` | `OUT_DIR`、`OUTSIDE.vhdr` | `OUT_DIR` 里必须已有 `03_filtered.set` 和 `04_cwl_corrected.set` |
| `denoise_report.py OUT_DIR` | `OUT_DIR` | 其余全是可选；但没有 `--outside` 就没有"相对磁体外"的指标 |

**强烈建议设置**（技术上可选，但不设会出问题或写到不该写的地方）
- `exclude`：数据里未使用的通道（你的数据是 `{'65','66'}`）。不排除的话，它们会被当作 EEG 目标通道参与回归和统计。
- `out_dir`：默认写到 `<vhdr 目录>/preprocessed`，也就是你的 BIDS 数据目录里；请改到 BIDS 之外，并且每组参数用独立目录。
- `rp_file`：有 fMRIPrep 头动文件时建议给，否则用 moving average 模板。

**其余全部可选**（不给就用默认值），按用途分组：

| 用途 | 字段 | 默认 |
|---|---|---|
| 数据识别 | `marker`、`cwl_regex`、`ecg_label` | `'R128'`、`'^CW'`、`'ECG'` |
| 梯度校正 | `rp_file`、`threshold`、`n_templates`、`span` | 无、0.3 mm、25、`[0 TR-1]` |
| 滤波与回归 | `band`、`target_fs`、`window_s`、`delay_s` | `[0.33 125]`、500、4、0.021 |
| 可选改进（均不是论文原方法，默认关闭） | `obs_npc`、`obs_hp`、`cwl_lowpass`、`bcg_aas`、`bcg_hr` | 0、70、0、0、空 |
| 安全检查 | `check_markers` | 1（开） |

**用法示例**
```matlab
% 1) 最小用法：只给文件，其余默认（输出在数据目录旁的 preprocessed，不推荐）
EEG = cwl_run_pipeline('/path/run-01_eeg.vhdr');

% 2) 论文原方法（无任何改进）
opts = struct('exclude', {{'65','66'}}, 'out_dir', '/path/out_paper');
EEG = cwl_run_pipeline(vhdr, opts);

% 3) 加头动信息
opts.rp_file = '/path/..._desc-confounds_timeseries.tsv';  opts.threshold = 0.15;

% 4) 推荐改进组合（OBS + 估计伪影低通）
opts.obs_npc = 4;  opts.obs_hp = 60;  opts.cwl_lowpass = 90;

% 5) 再加心跳锁定 AAS（输出 05_bcg_aas.set）
opts.bcg_aas = 1;  opts.bcg_hr = 64;
```
```bash
# 命令行：环境变量全部可选
CWL_OBS=4 CWL_OBS_HP=60 CWL_LOWPASS=90 CWL_RP=confounds.tsv CWL_THRESHOLD=0.15 \
    ./run_subject.sh RAW.vhdr OUT_DIR 65,66
```

### 6.1 `opts` 字段（`cwl_run_pipeline`，全部可选）

| 字段 | 默认 | 说明 |
|---|---|---|
| `out_dir` | `<vhdr 目录>/preprocessed` | 输出目录 |
| `marker` | `'R128'` | 体积标记类型 |
| `rp_file` | `[]` | SPM `rp_*.txt` 或 fMRIPrep `*.tsv`；为空则用 moving average 模板 |
| `threshold` | 0.3 | 头动阈值（mm，帧间平移） |
| `n_templates` | 25 | 每个 volume 的模板数 |
| `cwl_regex` | `'^CW'` | CWL 通道名的正则 |
| `ecg_label` | `'ECG'` | ECG 通道名 |
| `exclude` | `{}` | 既不是 EEG 目标也不是 CWL 的通道，如 `{'65','66'}` |
| `span` | `[0 TR-1]` | 相对标记的伪影范围（样本） |
| `band` | `[0.33 125]` | 带通（Hz） |
| `target_fs` | 500 | 回归采样率 |
| `window_s` | 4 | 回归窗长（s） |
| `delay_s` | 0.021 | 延迟嵌入半宽（s） |
| `obs_npc` | 0 | OBS 主成分个数，0 = 关；推荐 4（不在原论文中） |
| `obs_hp` | 70 | OBS 拟合的高通频率（Hz）；60 可同时处理 67 Hz 谐波 |
| `cwl_lowpass` | 0 | 对"估计出的伪影"做低通（Hz，8 阶零相位），0 = 关；推荐 90 |
| `bcg_aas` | 0 | 1 = 回归后再做心跳锁定 AAS，输出 `05_bcg_aas.set` |
| `bcg_hr` | 空 | 预期心率（bpm），帮助 R 波检测；可来自 Siemens `.puls` |
| `check_markers` | 1 | 1 = 体积标记间隔与 TR 偏差超过 1% 就报错；0 = 关闭检查 |

### 6.2 环境变量（`run_headless.sh` / `run_subject.sh`）

| 变量 | 对应 opts | 默认 |
|---|---|---|
| `MATLAB_BIN` | – | `matlab` |
| `CWL_RP` | `rp_file` | 无 |
| `CWL_THRESHOLD` | `threshold` | 0.3 |
| `CWL_OBS` | `obs_npc` | `run_subject.sh` 默认 4；`run_headless.sh` 默认 0 |
| `CWL_OBS_HP` | `obs_hp` | 70 |
| `CWL_LOWPASS` | `cwl_lowpass` | 关 |
| `CWL_BCG_AAS` | `bcg_aas` | 0 |
| `CWL_BCG_HR` | `bcg_hr` | 无 |
| `PYTHON` | – | `python`（`run_subject.sh` 里的 Python 步骤用） |

### 6.3 推荐设置

`CWL_OBS=4 CWL_OBS_HP=60 CWL_LOWPASS=90`。含义：OBS 去掉梯度残余的层频率谐波；CWL 估计的伪影在 90 Hz 以上不扣，避免把环路噪声带进 EEG。氦泵开着的数据应关闭 `CWL_LOWPASS`（90 Hz 以上的氦泵伪影不会被去除）。

## 7. 输出

每次运行先删除输出目录里旧的阶段文件，结束时写 `run_manifest.json`。

| 文件 | 内容 |
|---|---|
| `02_gradient_corrected.set/.fdt` | 梯度校正后，5000 Hz，全部通道（已剪掉首尾未校正部分） |
| `02b_obs_corrected.set/.fdt` | 开了 OBS 才有 |
| `03_filtered.set/.fdt` | 带通 + 500 Hz（CWL 回归前） |
| `04_cwl_corrected.set/.fdt` | **最终结果**（已裁成 CWL 实际校正的区间） |
| `05_bcg_aas.set/.fdt` | 开了 BCG AAS 才有 |
| `run_manifest.json` | 本次实际生成的阶段、输入、主要参数。评估脚本据此判断哪些阶段属于本次运行 |
| `settings.txt` | `run_subject.sh` 记录的设置 |
| `matlab_run.log` | 完整日志 |
| `rp_from_confounds.txt` | 用 fMRIPrep tsv 时转换出的 rp 文件 |

`.set` 可用 EEGLAB（`pop_loadset`）、MNE（`read_raw_eeglab`）或本包的 `cwl_loadset` 打开。

`EEG.etc.cwl_pipeline` 里记录了复现所需信息：输入 `.vhdr`、marker、全部原始标记位置、原始采样率和点数、裁剪起点（`first_kept_raw_sample`、`output_start_raw_s`）、`trimmed_to`、头动文件/阈值/事件数、模板数、回归参数（窗长、延迟、taper、方法、目标采样率）、通道正则、ECG 标签、排除列表、`bcg_hr`、OBS 参数与 `obs_frac`、`cwl_lowpass`、创建时间。

时间对应：输出第 1 个样本对应原始记录的 `output_start_raw_s` 秒（原始记录以 0 s 为起点）。

## 8. 运行时的安全检查

| 情形 | 行为 |
|---|---|
| 体积标记漏标、重复、间隔偏离 TR 超过 1% | 报错 `cwl:markers`，列出位置和间隔 |
| 头动行数与标记数不符 | 警告 |
| 通道单位为 V / mV / nV | 自动换算为 µV；未知单位警告 |
| 重跑到同一目录 | 先删旧阶段文件，写新的 `run_manifest.json` |
| CWL 回归没有校正任何样本（记录短于窗长） | 报错并说明 |
| 回归前 NaN | 报错 |

CWL 回归的有效区间取自原工具返回的窗口权重（权重 > 0），不依赖数值是否发生变化：原工具不校正开头约半个窗和结尾约 1–3 个窗，输出自动裁掉这些部分。

## 9. 诊断与评估

### 9.1 MATLAB 内置

```matlab
cwl_diagnose('/path/to/out')     % 写 diagnostics.png/.json/diagnostics_data.mat
```
内容：回归前后平均谱、CWL 谱、10 s 示例、每通道去除方差、频段功率比，以及随机相位对照（过拟合下限）。

`cwl_eoec_metrics.m` 实现论文的 EO/EC 指标（谱距离、Jaccard、地形距离），需要同一被试磁体外的睁眼和闭眼段。

### 9.2 与磁体外记录比较（Python 脚本，对 MATLAB 输出运行）

```bash
./evaluate_subject.sh OUT_DIR OUTSIDE_run-00.vhdr
# 可选：CWL_RAW_VHDR=RAW.vhdr 在去噪报告中加入未校正的 raw 阶段；CWL_BCG_HR=64
```
写入 `OUT_DIR`：

| 文件 | 内容 |
|---|---|
| `inside_vs_outside.png` / `.channels.csv` | 示例轨迹、平均谱、每通道 1–30 Hz RMS |
| `outside_comparison.png` / `.json` | 谱距离（论文式 4）、频段差、alpha 地形 |
| `bcg_residual.png` / `.json`、`bcg_topomap.png` | 心跳锁定残差（03、04、[05]、磁体外）和头皮图 |
| `removed_signal.json` | CWL 去除成分与邻近电极的相干性 |
| `denoise_spectra.png` | **每一步**（raw/02/02b/03/04/05/磁体外）的平均谱 + 相对磁体外的 dB 差 |
| `denoise_indicators.png/.json/.csv` | 去噪指标（第 10 节） |
| `denoise_topomaps.png` | 每一步的头皮分布图 |
| `evaluation.txt` | 全部终端输出 |

只生成去噪报告：
```bash
python python/denoise_report.py OUT_DIR --outside OUTSIDE.vhdr [--vhdr RAW.vhdr] [--hr 64] [--slice-freq 16.84] [--marker R128] [--crop START END]
```

其他诊断：`analyze_residual_peaks.py`（残峰是否与体积锁定、与层频率对应）、`check_volume_alignment.py`（原始数据逐 volume 时间错位和削波）、`plot_fig3.py`（论文图 3 风格）。用法见 `docs/MANUAL_PYTHON.md` 第 7 节。

## 10. 如何看去噪指标

| 指标 | 含义 | 好的方向 |
|---|---|---|
| `rms_1_30_uV` | 1–30 Hz RMS（通道中位数） | 接近磁体外 |
| `spectral_distance_pct`、`..._1_30Hz_pct` | 论文式 4，与磁体外平均谱的相对距离 | 越小越好；全频带易被高频底噪和 50 Hz 拉高，1–30 Hz 更敏感 |
| `band_power`、`dB_vs_outside_*` | 7 个频段功率及相对磁体外的 dB | 越接近 0 dB 越好 |
| `slice_peak_db` | 层频率谐波高出局部背景多少 dB | 越小越好（0 = 无） |
| `bcg_locked_rms_uV` / `bcg_residual_pct` | R 波锁定平均的 RMS 及占 03 的百分比 | 越小越好；磁体外值是下限 |
| `step_change_db_1_30` | 本步相对上一步的 1–30 Hz 功率变化 | 显示每一步的作用 |
| `cwl_variance_removed` | 1 − var(04)/var(03) | 仅作参考 |
| `log_spectrum_r` | 与磁体外对数谱的相关（1–45 Hz） | 越接近 1 越好 |
| `alpha_topo_r` | 8–13 Hz 功率地形与磁体外的相关 | 越接近 1 越好 |

解读时注意：
- 磁体外参考不是纯净真值：常有 50 Hz 工频峰和 <1 Hz 漂移，状态也与磁体内不同。
- 磁体外心跳锁定 RMS（约 1 µV）含心跳诱发电位，不是零。
- 头皮图在电极位置之外是插值外推，最外圈不是实测值。
- 去掉的方差比例高不等于去噪好，需结合随机相位对照和 alpha 是否被过度去除。

## 11. 与 Python 版的一致性检查

```bash
python python/make_crosscheck_data.py crosscheck_data.mat
```
```matlab
setup_cwl_pipeline; cwl_crosscheck('crosscheck_data.mat')
```
8 项：模板矩阵、梯度校正、OBS（去除比例与输出波形）、R 波检测、BCG AAS、CWL 回归（完整权重内部区域）、5000→500 Hz 降采样。任何 FAIL 都会报错。**不覆盖**：BrainVision / .set 读写、事件、带通、CWL 回归的首尾样本、整条流程端到端。

## 12. 文件清单

| 文件 | 作用 |
|---|---|
| `setup_cwl_pipeline.m` | 设置路径并检查 Signal Processing Toolbox |
| `cwl_run_pipeline.m` | 主流程 |
| `cwl_load_brainvision.m` | 读取 BrainVision（含通道单位） |
| `cwl_saveset.m` / `cwl_loadset.m` | 写 / 读 .set + .fdt |
| `cwl_resample.m` | 降采样（与 EEGLAB `pop_resample` 相同算法） |
| `cwl_gradient_correct.m` | 梯度模板减法 |
| `cwl_obs_residual.m` | OBS 残差去除 |
| `cwl_sosfiltfilt.m` | 零相位二阶节滤波 |
| `cwl_cwregression.m` | 调用原作者的 CWL 回归，并返回窗口权重 |
| `cwl_bcg_aas.m` | 心跳锁定 AAS |
| `cwl_diagnose.m`、`cwl_eoec_metrics.m` | 诊断与论文指标 |
| `cwl_crosscheck.m` | 与 Python 版对照 |
| `run_headless.sh`、`run_headless_job.m` | 无界面运行 |
| `debug_cwl_env.m`、`fix_matlab_path.m` | 环境检查 / 修复路径 |
| `vendor/` | 第三方代码（未修改）：CWRegrTool（MIT）、BergenToolboxModified（GPL v2）、firfilt（GPL） |

## 13. 常见问题

| 现象 | 处理 |
|---|---|
| `cwl:markers` 报错 | 看提示中的位置和间隔，核对 `.vmrk`；确认无误可设 `opts.check_markers = 0` |
| Signal Processing Toolbox 检查失败 | MATLAB 路径损坏：`fix_matlab_path`，或 `debug_cwl_env` |
| 评估时警告"没有 run_manifest.json" | 这是旧目录，里面的 02b/05 可能是旧结果；用独立目录重跑 |
| 磁盘不够 | 确认结果后删除 02、02b（各约 1 GB） |
| 评估脚本读不到 `.set` | 需要装 `mne`、`pymatreader`；`PYTHON` 环境变量指向对应环境 |
| `denoise_report` 的头皮图被跳过 | 至少 16 个通道要有标准 10-05 位置 |
| 缺 `diagnostics_data.mat` | 只是谱图少一条随机相位对照线 |

## 14. 局限

- 原论文方法之外的步骤（OBS、估计伪影低通、BCG AAS）都是可选的、默认关闭或需要显式设置的改进，使用时请在方法部分注明（Niazy 2005；Allen 1998）。
- 与 TR 锁相的诱发活动可能被 OBS 部分去除，任务态实验需注意。
- 90 Hz 以上的伪影（氦泵开启时）不会被去除。
- 没有实现 Bergen 的基线校正和削波处理（检查数据无削波即可）。
- 未实现论文里的 OBS / OBS-ICA 对比方法和 EEG–BOLD 分析。
- 覆盖范围之外的电极（如左后颞边缘）的残余 BCG 较大，可用 BCG AAS 进一步处理，但会同时去掉心跳锁相的脑活动。

## 参考文献
- Allen, P. J., et al. (1998/2000). *NeuroImage*.
- Moosmann, M., et al. (2009). *NeuroImage*, 45(4), 1144–1150.
- Niazy, R. K., et al. (2005). *NeuroImage*, 28(3), 720–737.
- van der Meer, J. N., et al. (2016). *NeuroImage*, 125, 880–894.

（DOI 请在投稿前核对。）
