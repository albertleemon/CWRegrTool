# CWL 流程：独立 MATLAB 版

这个文件夹自带全部依赖，可以单独复制到任何位置使用。

**只需要：** MATLAB（已在 R2024b 上测试）+ Signal Processing Toolbox。
**不需要：** EEGLAB、bva-io、firfilt 插件、Statistics Toolbox、git 或网络。

## 使用

```matlab
cd /path/to/matlab        % 本文件夹
setup_cwl_pipeline        % 只设置路径，不下载任何东西
opts = struct('exclude', {{'65','66'}}, 'obs_npc', 4, 'obs_hp', 60, 'cwl_lowpass', 90, ...
              'out_dir', '/path/to/out');
EEG = cwl_run_pipeline('/path/to/sub-XXX_task-rest_run-01_eeg.vhdr', opts);
cwl_diagnose(opts.out_dir);                     % 回归前后对比 + 随机相位对照
```

不开界面运行：`./run_headless.sh RAW.vhdr OUT 65,66`，参数用环境变量传入（`CWL_OBS`、`CWL_OBS_HP`、`CWL_LOWPASS`、`CWL_BCG_AAS`、`CWL_RP` 等，见 `../docs/RUN_LOCALLY.md`）。

## 输出
`02_gradient_corrected`、`02b_obs_corrected`、`03_filtered`、`04_cwl_corrected`、`05_bcg_aas`（后者可选），每个都是 `.set` + `.fdt`。格式就是 EEGLAB 数据集，EEGLAB（`pop_loadset`）、MNE（`read_raw_eeglab`）和 `cwl_loadset` 都能打开。

## 文件

| 文件 | 作用 |
|---|---|
| `setup_cwl_pipeline.m` | 设置路径并检查 Signal Processing Toolbox |
| `cwl_run_pipeline.m` | 主流程：读入 → 梯度校正 → OBS（可选）→ 滤波与降采样 → CWL 回归 → BCG AAS（可选） |
| `cwl_load_brainvision.m` | 读取 BrainVision（.vhdr/.vmrk/.eeg），结果与 bva-io 的 `pop_loadbv` 相同 |
| `cwl_saveset.m` / `cwl_loadset.m` | 写 / 读 EEGLAB 数据集（.set + .fdt） |
| `cwl_resample.m` | 降采样，算法与 EEGLAB `pop_resample` 相同（Kaiser 窗抗混叠滤波 + `resample`） |
| `cwl_gradient_correct.m` | 梯度伪影模板减法（Bergen `CorrectionMatrix` 的无界面实现） |
| `cwl_obs_residual.m` | OBS 去除梯度残差（Niazy 2005） |
| `cwl_sosfiltfilt.m` | 零相位二阶节滤波 |
| `cwl_cwregression.m` | 调用原作者的 CWL 回归（与 `pop_cwregression` 的计算相同） |
| `cwl_bcg_aas.m` | 可选：CWL 后的心跳锁定 AAS（Allen 1998） |
| `cwl_diagnose.m`、`cwl_eoec_metrics.m` | 诊断与论文指标 |
| `cwl_crosscheck.m` | 与 Python 版逐模块对照（需 `python/make_crosscheck_data.py`） |
| `debug_cwl_env.m`、`fix_matlab_path.m` | 环境检查 / 修复损坏的 MATLAB 路径 |
| `run_headless.sh`、`run_headless_job.m` | 无界面运行 |

### `vendor/`：随附的第三方代码（均未修改）
| 文件夹 | 来源 | 许可 |
|---|---|---|
| `vendor/CWRegrTool` | van der Meer et al. (2016) 的 CWL 回归核心：`run_cwregression`、`m_do_taperedhann`、`+misc`、`+tools` 等 | MIT |
| `vendor/BergenToolboxModified` | 头动模板选择 `m_rp_info`、`m_moving_average`（Moosmann 2009） | GPL v2 |
| `vendor/firfilt` | 抗混叠滤波器设计 `firws`、`firwsord`、`windows`（EEGLAB firfilt） | GPL |

## 运行安全检查与记录
- 体积标记间隔与 TR 偏差超过 1% 时（漏标、重复、多余标记）直接报错；确认无误可设 `opts.check_markers = 0`（Python：`--no-marker-check`）。
- 每次运行先删除输出目录里旧的阶段文件（02/02b/03/04/05），结束时写 `run_manifest.json`（本次实际生成的阶段）；评估脚本只采用清单中列出的阶段。仍建议每组参数使用独立输出目录。
- BrainVision 通道单位（V/mV/µV/nV）按头文件换算为 µV。
- CWL 回归的有效区间取自原工具返回的窗口权重（权重 > 0），不再依赖数值是否变化。
- `EEG.etc.cwl_pipeline` 记录输入 .vhdr、原始标记位置、裁剪起点（`first_kept_raw_sample`、`output_start_raw_s`）、头动文件/阈值/事件数、模板数、回归参数等。

## 与旧版（依赖 EEGLAB）的一致性
在同一份合成 BrainVision 记录（16 通道，5000 Hz，40 个 volume，含梯度伪影、心冲击伪影和 CWL 信号）上，分别运行旧版（EEGLAB + bva-io + firfilt）和本版，参数为 OBS 4 / 60 Hz、低通 90 Hz、BCG AAS 开：
02、02b、03、04、05 五个输出的数据差值均为 **0**，通道、采样率、事件（含裁剪处的 boundary）也相同。
如果 MATLAB 路径里还有旧的 EEGLAB 或 CWRegrTool，也不会冲突：`setup_cwl_pipeline` 会把本文件夹放在路径最前面，那些同名文件内容也相同。


完整使用手册：`../docs/MANUAL_MATLAB.md`（参数、输出、评估与指标解读、常见问题）。去噪报告由 `python/denoise_report.py` 生成（MATLAB 流程之后运行，见手册第 9 节）。
