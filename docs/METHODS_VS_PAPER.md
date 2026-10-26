# 本流程与原论文方法的区别和改进

对照文献：van der Meer, J. N., et al. (2016). *Carbon-wire loop based artifact correction outperforms post-processing EEG/fMRI corrections—A validation of a real-time simultaneous EEG/fMRI correction method.* NeuroImage, 125, 880–894。（作者和 DOI 请投稿前自行核对。）

> 说明：下面对论文的描述只包含本仓库文档（`docs/EEG_fMRI_CWL_v2.2.md`、`docs/IMPROVEMENTS.md`、`docs/PIPELINE.md`）和代码里已经记录并使用的内容，没有重新逐段核对论文原文。具体数值（窗长、阈值、模板数、报告的降幅等）写进论文或基金文本前，请对照原文再确认一遍。

---

## 1. 原论文做了什么

### 1.1 问题

同步 EEG–fMRI 的 EEG 里叠加了三类大伪影：

| 伪影 | 来源 | 特点 |
|---|---|---|
| 梯度伪影 | MRI 梯度切换 | 幅度极大，与扫描严格同步，每个 volume 形状几乎相同 |
| BCG（心冲击）伪影 | 心跳带动头和电极在磁场中运动 | 与心跳相关，但每次形状不同，不与扫描同步 |
| 氦泵 / 振动伪影 | 氦压缩机振动传到电极和导线 | 宽频、窄峰，与扫描不同步 |

传统做法（AAS、OBS、OBS-ICA 等）是在记录之后用模板或成分分解去除 BCG，需要 ECG 的 R 波，并假设伪影形状在一段时间内稳定。

### 1.2 核心想法

BCG 和氦泵伪影的共同本质是**电极和导线在磁场中的微小运动**。碳线环（carbon-wire loop，CWL）不接触头皮，只感应运动，不含脑电。因此把 EEG 对 CWL 做线性回归：能被 CWL 线性解释的部分就是运动伪影，减掉即可。这样既不需要 ECG 的 R 波，也不需要模板，原则上可以实时进行。

### 1.3 方法（论文的具体配置）

1. **梯度伪影**：先用 Bergen EEG-fMRI 工具箱（EEGLAB 插件，图形界面）做模板减法，采用头动信息引导的模板选择（Moosmann et al., 2009）：每个模板用邻近 25 个伪影片段，帧间位移超过 0.3 mm 时重置模板缓冲区。EEG 和 CWL 通道一起校正。
2. **滤波**：EEG 和 CWL 同一个 0.33–125 Hz 带通。
3. **CWL 回归**：4 s Hann 窗，相邻窗重叠 50%（窗权重之和为 1）；每个 CWL 取 ±21 ms 的时间延迟，回归量由延迟后的 6 个 CWL 组成；每窗对每个 EEG 通道做最小二乘，窗间叠加后得到整条伪影轨迹，从 EEG 中减去。窗口很短，所以系数随头姿和氦泵状态自适应。
4. **采样率**：500 Hz（每个 CWL 约 21 个延迟，共约 126 个回归量；5000 Hz 下约 1266 个）。

### 1.4 验证

用同一被试磁体外的 EEG 作参考，比较多种校正方法：
- 平均功率谱的**相对谱距离**（论文式 4）；
- 睁眼/闭眼的 alpha 增幅，用 **Jaccard 重叠**（式 5）和**地形图距离**（式 6）；
- 在 He 泵开 / 关两种条件下，与 AAS、OBS、OBS-ICA 比较；
- 第三部分是 EEG–BOLD 分析（本流程不包含）。

本仓库文档记载：论文报告 1–20 Hz 附近约 10 dB 的降幅（我们在 sub-BOR5A 上复现到 8–14 dB）。

---

## 2. 本流程保持不变的部分

| 部分 | 状态 |
|---|---|
| CWL 回归核心 | **直接调用原作者的 `run_cwregression` 代码**，4 s Hann 窗、±21 ms、taper 1、重叠 50%，不做修改（MATLAB 版）。Python 版是按论文和原代码移植的 |
| 梯度校正算法 | 与 Bergen 的模板选择同一逻辑（`m_rp_info` / `m_moving_average` 直接使用原文件；模板由 1 Hz 高通副本平均，相关性差的模板按 2.5×中位数剔除） |
| 滤波 | 0.33–125 Hz 带通，EEG 与 CWL 同一滤波器 |
| 采样率 | 500 Hz |
| 校正对象 | 梯度校正作用于所有通道（含 CWL、ECG）；CWL 回归只作用于 EEG 通道 |
| 评估指标 | 论文式 4、5、6 已实现（`metrics.py`、`cwl_eoec_metrics.m`） |

---

## 3. 区别与改进总览

图例：**A** = 实现方式改进（算法不变）；**B** = 新增方法（论文没有，需在方法部分注明）；**C** = 验证和工程改进；**D** = 未做。

| 步骤 | 原论文 | 本流程 | 类别 | 默认 |
|---|---|---|---|---|
| 1. 梯度校正 | Bergen GUI | 同一算法的无界面实现，可批处理；支持 fMRIPrep 头动文件 | A | 开 |
| 1a. 头动阈值 | 0.3 mm | 阈值可调（`threshold`）；本数据最大帧间位移 0.293 mm，0.3 下无事件，等同 moving average | A | 0.3 |
| 1b. 梯度残差 | 无 | **OBS 残差 PCA**（Niazy 2005），4 个主成分，60 Hz 以上 | B | 关（推荐开） |
| 2. 滤波与降采样 | 0.33–125 Hz；文中以 500 Hz 计 | 相同，并明确降到 500 Hz | A | 开 |
| 3. CWL 回归 | 原代码 | 原代码不变；自动裁掉未校正首尾 | A | 开 |
| 3b. 估计伪影低通 | 无 | **对回归估计出的伪影做 90 Hz 低通**后再减 | B | 关（推荐开） |
| 3c. BCG AAS | 仅作对比方法 | CWL 后可选再做心跳锁定 AAS | B | 关 |
| 4. 评估 | 与磁体外记录比较 | 已实现论文指标；新增随机相位对照、残峰分析、逐 volume 检查、**分阶段谱和去噪报告** | C | – |
| 5. 运行安全 | – | 运行清单、标记间隔检查、单位换算、复现信息 | C | 开 |
| 6. OBS / OBS-ICA 对比、EEG–BOLD | 有 | 未实现 | D | – |

---

## 4. 改进细节

### 4.1 梯度校正：无界面实现与头动输入（A）

- **原论文**：Bergen EEGLAB 插件的图形界面，每次手工操作，不便批量和复现。
- **本流程**：`cwl_gradient_correct.m` 与 Bergen `CorrectionMatrix` 同一逻辑，无界面；直接读取 fMRIPrep 的 `*_desc-confounds_timeseries.tsv`（`trans_x/y/z`、`rot_x/y/z`，单位 mm 与弧度，与 SPM `rp` 相同）。
- **自动裁掉未校正的首尾**：校正范围（第一个标记到最后一个标记 + TR）以外没有模板，仍含梯度伪影，输出直接删去。
- **省略的部分**：Bergen 的基线校正和削波处理。检查本数据无削波（最大 10.5 mV，量程 16.4 mV）。
- **效果（sub-BOR5A）**：层频率谐波（16.84 Hz 及其倍数）在 63 个 EEG 通道平均谱上降低约 4–5 个数量级。

### 4.2 OBS 残差去除（B，新增，不在论文中）

- **问题**：模板减法后，84.5、101、117.75 Hz 处仍有残峰，高出背景 10–14 dB；CWL 回归只能降 3–6 dB。
- **诊断**：
  - 峰间距与 `SliceTiming` 对应的层频率一致（16.86 Hz vs 16.84 Hz），精细结构落在 1/TR 网格上 → 是梯度伪影的层频率谐波；
  - 逐 volume 时间错位标准差 ≤ 0.005 个采样点（约 1 µs）→ 不是时间抖动；
  - 0 个样本到达量程 → 不是削波；
  - 与头动（FD）相关 r = 0.05 → 不是头动；
  - 少数 volume 的残差高出 10–60 倍，所有通道（含 CWL）同时出现 → **伪影形状在部分 volume 中变化，平均模板无法跟踪**。
- **方法**：Niazy et al. (2005) 的最优基集（FASTR 的 OBS 步骤）。对每个通道，把模板减法后的残差高通到 60 Hz，按 volume 做 PCA，去掉前 4 个主成分；残差异常大的 volume（> 100× 中位数）不参与估计。
- **效果**：84.5 / 101 / 117.75 Hz 处回归后的峰高从 6.9 / 8.5 / 6.7 dB 降到 2.5 / 1.7 / 2.3 dB；0–125 Hz 内不再有高于背景 3 dB 的峰；1–30 Hz 功率变化 −0.01 dB，30–60 Hz −0.03 dB。
- **风险**：与 TR 锁相的诱发活动（刺激严格按 TR 呈现）可能被部分去除，任务态实验需注意。
- **开关**：`obs_npc`（0 = 关）、`obs_hp`（60 或 70 Hz）；Python：`--obs-npc`、`--obs-hp`。使用时请在方法部分注明 Niazy et al. (2005)。

### 4.3 对估计伪影低通（B，新增，不在论文中）

- **问题**：CWL 回归后，90 Hz 以上的功率反而比回归前高 2–3 dB，101、118 Hz 又出现小峰。
- **原因**：回归系数由能量占主导的低频（BCG）决定，同一组系数也作用于碳线环自身的高频噪声；EEG 中没有这些成分，减去它们等于把它们加进 EEG。随机相位对照中同样可见。本被试氦泵伪影很弱，高频几乎没有可去的伪影。
- **第一次尝试失败**：对碳线环信号本身低通，会使带延迟的回归变得病态，在过渡带（80–105 Hz）放大噪声，已弃用。
- **最终方法**：回归照常进行，只对**估计出的伪影**做 8 阶零相位低通（`cwl_lowpass = 90`），再从 EEG 减去。
- **效果**：90–125 Hz 与回归前重合；alpha 功率比 0.0463 → 0.0462，1–7 Hz 功率比 0.1061 → 0.1061（基本不变）。
- **限制**：90 Hz 以上的伪影（氦泵开启时的 60–120 Hz 氦泵伪影）不会被去除。**氦泵开着的数据应关闭此选项。**
- **开关**：`cwl_lowpass`（0 = 关）；Python：`--artifact-lowpass`。

### 4.4 CWL 后的 BCG AAS（B，可选，默认关闭）

- **原因**：残余 BCG 集中在碳线环覆盖不到的电极（左后颞边缘最明显）。我们在一次测试运行上看到：心跳锁定残差 13.7 → 2.68 µV（磁体外约 1.04 µV），其中 TP9、TP7、P7、PO7、P5 占全部超出量的约 46%。
- **方法**：Allen et al. (1998) 平均伪影减法，R 波从 ECG 检测（10–35 Hz 带通、差分平方、80 ms 滑动平均，最小峰间距 0.7×心动周期，对扫描仪内高大 T 波稳健）；每个心跳减去最近 21 个心跳的平均。
- **代价**：只去除与心跳严格锁定的部分，**心跳锁相的脑活动（如心跳诱发电位）也会被去除**；真实数据上的效果尚需验证。
- **开关**：`bcg_aas = 1`（可加 `bcg_hr`）；Python：`--bcg-aas`、`--bcg-hr`。输出 `05_bcg_aas.set`。

### 4.5 CWL 回归的有效区间（A）

- **原工具行为**：回归不校正开头约半个窗和结尾约 1–3 个窗（实测不校正开头约 2 s、结尾约 5 s）。
- **本流程**：输出只保留原工具实际校正过的区间，区间由原工具返回的**窗口权重**（权重 > 0）确定，而不是"数值是否发生变化"——后者在零 CWL、零拟合或短记录时会得到空区间，且边缘恰好拟合为零会让裁剪位置漂移。

### 4.6 验证和诊断（C，新增）

| 工具 | 内容 | 论文中 |
|---|---|---|
| 随机相位 CWL 对照（`cwl_diagnose`） | 用相位随机化的 CWL 重做回归，作为过拟合下限 | 无 |
| 残峰分析（`analyze_residual_peaks.py`） | 残峰是否体积锁定、是否在 1/TR 网格上、逐 volume 和逐通道分布、与头动相关 | 无 |
| 逐 volume 检查（`check_volume_alignment.py`） | 原始数据上的时间错位和削波 | 无 |
| 心跳锁定残差（`check_bcg_residual.py`） | R 波锁定平均的残差及其头皮图，以磁体外为下限 | 无 |
| 去除成分的邻近电极相干（`check_removed_signal.py`） | CWL 是否误去脑电 | 无 |
| 论文指标 | 谱距离（式 4）、EO/EC alpha 的 Jaccard（式 5）、地形图距离（式 6） | 有 |
| **去噪报告（`denoise_report.py`）** | 每一步的平均谱、相对磁体外的 dB 差、指标表、每一步的头皮分布图、后部 alpha 峰（扣 1/f） | 无 |

去噪报告的指标：1–30 Hz RMS、全频带和 1–30 Hz 谱距离、7 个频段功率及相对磁体外 dB、层频率谐波残峰、心跳锁定残差、每步 1–30 Hz 功率变化、CWL 去除方差比例、对数谱相关、alpha 地形相关，以及后部 alpha 峰频率、峰高和相对磁体外的功率。

**论文式的 EO/EC Jaccard（式 5）和地形距离（式 6）** 需要同一被试磁体内外都有睁眼和闭眼的段。我们的两段记录都是闭眼静息，所以只能直接比较 alpha 峰，做不了反应性比较。

### 4.7 运行安全与可复现（C，新增）

- **运行清单** `run_manifest.json`：每次运行先删除输出目录里旧的阶段文件，结束时写入本次实际生成的阶段；评估脚本只采用清单里列出的阶段，避免旧的 `02b` / `05` 被当成本次结果。仍建议每组参数使用独立输出目录。
- **标记检查**：体积标记间隔偏离 TR 超过 1%（漏标、重复、多余）直接报错；否则漏掉的那一段没有模板校正却会被保留。
- **BrainVision 单位**：按通道单位（V / mV / nV / µV）换算为 µV。
- **复现信息**：`EEG.etc.cwl_pipeline` 记录输入 `.vhdr`、全部原始标记位置、裁剪起点（`first_kept_raw_sample`、`output_start_raw_s`）、头动文件/阈值/事件数、模板数、回归参数、OBS 参数、创建时间等，便于与原始 EEG/fMRI 对齐。
- **Python 版输出**保留事件、通道类型（EEG/ECG/其余）和电极位置。
- **MATLAB 版独立化**：只需 MATLAB + Signal Processing Toolbox，不需要 EEGLAB。在合成 BrainVision 记录上与旧版（EEGLAB + bva-io + firfilt）逐阶段比较，02、02b、03、04、05 五个输出的数据差值均为 0。
- **跨语言检查**：`cwl_crosscheck` 覆盖模板矩阵、梯度校正、OBS（去除比例与输出波形）、R 波检测、BCG AAS、CWL 回归内部区域、降采样共 8 项；不覆盖 I/O、事件、带通和整条流程端到端。

---

## 5. 一次测试运行的结果（仅供参考，单个被试）

| 指标 | 03（回归前） | 04（CWL 后） | 磁体外 |
|---|---|---|---|
| 1–30 Hz RMS（µV） | 20.5 | 8.27 | 8.42 |
| 1–30 Hz 谱距离（%） | 436 | 69 | 0 |
| 心跳锁定 RMS（µV） | 13.7 | 2.68（占 20.7%） | 1.04 |
| CWL 去除方差比例 | – | 77% | – |
| 对数谱相关 | 0.905 | 0.939 | 1 |

仍存在的差异：60–125 Hz 比磁体外高约 4.5–6.7 dB（各步骤都没改善，更可能是磁体内的本底噪声）；8–13 Hz 比磁体外低约 2.6 dB（状态不同或被过度去除，尚需用 alpha 峰指标和 `check_removed_signal` 区分）；磁体外参考本身有 50 Hz 工频峰和 <1 Hz 漂移。这些是单个被试的结果，不能直接推广。

---

## 6. 尚未完成 / 与论文仍有差距

1. 论文作为对比的 BCG 方法：AAS 已作为可选步骤实现；**OBS 和 OBS-ICA 未实现**。
2. **EEG–BOLD 分析**（论文第三部分）未做。
3. 没有睁眼段时无法计算论文的 Jaccard（式 5）和地形距离（式 6）。
4. Python 流程尚未在真实数据上运行。
5. BCG AAS 在真实数据上的效果尚需验证（可能连心跳锁相的脑活动一起去掉）。
6. 90 Hz 以上的伪影（氦泵开启）不会被低通选项去除。
7. 50 Hz 工频峰仍在，可在后续分析中陷波。

---

## 7. 方法部分的建议写法

> EEG–fMRI data were preprocessed following van der Meer et al. (2016) [原方法]: gradient artifacts were removed from all channels (EEG, ECG and carbon-wire loops) by head-motion-informed average-artifact subtraction (25 templates, 0.3 mm displacement threshold; Moosmann et al., 2009; Bergen EEG-fMRI toolbox), followed by identical zero-phase band-pass filtering (0.33–125 Hz) of EEG and loop signals, downsampling to 500 Hz, and windowed carbon-wire-loop regression (4 s Hann windows, 50% overlap, ±21 ms delays; original CWRegrTool code). Departures from the published pipeline (all parameter values and the choice to enable each step reported in the Supplement): (i) the residual gradient artifact (slice-frequency harmonics) was further reduced by an optimal-basis-set residual PCA (4 components, >60 Hz; Niazy et al., 2005); (ii) the artifact estimated by the regression was low-pass filtered at 90 Hz before subtraction, so that loop noise was not added to the EEG above this frequency; (iii) [optional] remaining heartbeat-locked residuals were removed by average artifact subtraction (Allen et al., 1998). Data quality was assessed against an eyes-closed resting recording from the same participant outside the scanner (relative spectral distance, band power, heartbeat-locked residual, posterior alpha peak).

（把未启用的步骤从上面的句子里删掉；括号里的参考文献请核对后再引用。）

---

## 8. 参考文献

- Allen, P. J., Polizzi, G., Krakow, K., Fish, D. R., & Lemieux, L. (1998). Identification of EEG events in the MR scanner: the problem of pulse artifact and a method for its subtraction. *NeuroImage*, 8, 229–239.
- Allen, P. J., Josephs, O., & Turner, R. (2000). A method for removing imaging artifact from continuous EEG recorded during functional MRI. *NeuroImage*, 12(2), 230–239.
- Moosmann, M., et al. (2009). Realignment parameter-informed artefact correction for simultaneous EEG–fMRI recordings. *NeuroImage*, 45(4), 1144–1150.
- Niazy, R. K., et al. (2005). Removal of FMRI environment artifacts from EEG data using optimal basis sets. *NeuroImage*, 28(3), 720–737.
- van der Meer, J. N., et al. (2016). Carbon-wire loop based artifact correction outperforms post-processing EEG/fMRI corrections—A validation of a real-time simultaneous EEG/fMRI correction method. *NeuroImage*, 125, 880–894.

（以上条目为我根据记忆整理，页码和作者请在投稿前核对。）
