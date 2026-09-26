# ALOD 精简子项目迁移计划书


## Execution update (2026-09-22)

The user's current instruction supersedes earlier fixed physical-core/BLAS-one
requirements for production execution. Production launchers now use no CPU
binding, no hard 64-thread cap, and no forced BLAS thread count. Explicit test
thread settings remain reproducibility controls. Memory and disk supervision
remain in force.

Implement manual `ell_absolute_threshold=tau(k)` as an additional OR promotion
criterion on the existing lazy schedule. Integrate background audit workers
through the same executable, with immutable snapshots and resumable audit work.
Save training/ell/accepted boundaries automatically; recover interrupted output
transactions and reject duplicate writers. See `docs/async_execution.md` for the
exact contract and limitations. Full paper trajectories and remaining scientific
diagnostics are still deferred; bounded server tests are not full reproduction.


## Version 4 update (2026-09-22)

This update supersedes conflicting historical numerical settings below; see
`docs/paper_v4_review.md` for the complete English review and remaining work.
Current enriched tests use `J_ell = I - T_ell I_H`, with the full coupled PG
blocks. The current E1/E2 presets both have 25 cycles and two reference sweeps
(51 accepted states), normalized nominal balance threshold 0.3, and ell 2–4.
E3 adds k=128 and uses coarse levels 4/6/8/10/12 with reference gaps 3/4/5/6/7.
Raw-ratio/AOT policies and the 33-state E2 campaigns are historical controls.
Old AOT triangular optimizations do not apply to the kernel-lifted equations.

P7 must integrate the remaining generic pre-training defect controller,
exact-only/target-terminal audits, large-checkpoint and solver-memory recovery,
remaining endpoint diagnostics and portable result import. P8 must validate
current kernel-lifted E2, normalized E1, all E3 wavenumbers and complete audit
coverage on the server. Long experiments remain explicitly deferred; supplied
results must not be presented as reruns of this checkout. P9 release and the
current manuscript figure/table pipeline remain open. Preserve the original
PowerShell/SSH, tmux, physical-core and memory supervision requirements.


初稿：2026-09-09；本次核对更新：2026-09-11。保留原文件名，避免既有任务链接失效。历史状态：计划已按最新版论文更新。2026-09-11 已完成 P0/P1 验收；当前执行状态以 `docs/migration_status.md` 为准。

语言约定：本计划书按用户要求暂时保留中文；README、NOTICE 和 `docs/` 下的文档使用英文。所有提交内容（包括本计划书）仍不得包含本地绝对路径。

P0/P2 核对补充：论文规定准插值采用按单元面积加权的顶点平均，而归档源码采用相邻单元的算术平均。二者在均匀等面积网格上一致，在非均匀网格上存在已复现差异；反例见 `tests/quasi_average_probe.cpp`。P2 还确认归档 Riesz 插值支撑扩张规则在回归网格上生成 N³ 顶点 patch，而论文指定 N²。P2 默认采用论文的面积加权与 N² 定义，并分别保留显式归档兼容策略；旧版等价验收选用归档策略，论文策略另以独立核基与谱计算验收。P1 的 SLOD 保持原均匀网格行为。映射与数值证据见 `docs/paper_implementation_map.md`、`docs/provenance/p2_validation.json`，不得将历史轨迹改称论文策略的结果。


本次依据：`<PAPER_ROOT>/helmholtz_lod_certified_amsart_revised.tex`，修改时间 `2026-09-11 12:35:58`，SHA-256 `f714fb0ba66c2cc276437002aa2306d16271fd91fc88205ff8c3ca07ccdcbc37`。执行 agent 开始前须重新核对；若论文变化，先更新差异清单，不能沿用旧图号或旧数据表。

本版主要调整：将 E1 四组 nominal/family 标记对照及全状态载荷审计列入必需复现范围；将已从最新版删除的 E2 shared/fresh 图和五样本表降为可选诊断；补全局部能量算子、AOT、POD 和升层启发式的数学约束；加入最近发现的大字段 CSV 解析、审计内存和恢复数据来源问题。

目标：从现有 `LOD2d_C++` 中提取一个能够独立构建、运行、诊断和复现论文 E1/E2 的 `ALOD` 项目。包括 E1 主生产轨迹、E1 四组标记对照、E2 区域 AS 继承生产轨迹，以及 AFEM、UFEM、SLOD 基线。family-AFEM 是同一 AFEM 驱动的标记成员配置，不另造算法或问题；不建设通用 FEM/Helmholtz 软件平台。

建议目标目录为 `<ALOD_ROOT>`；实施时统一称 `<ALOD_ROOT>`。本计划及来源清单位于旧项目，迁移后应复制到新项目 `docs/`。新项目不得运行时引用旧项目的头文件、静态库、构建目录或绝对路径。

## 1. 给执行 agent 的任务说明

按本文 P0–P9 的依赖顺序实施。先冻结来源和得到可比较的旧版结果，再提取数值模块；先保持数学行为，再优化性能。每阶段提交代码、实际执行命令、验证结果和未完成项。跨会话用 `docs/migration_status.md` 续接，不凭上一条聊天摘要猜测实验参数。

执行约束：

- 保留旧项目、论文及所有原始实验数据；在新目录工作，不就地删除旧模块，不覆盖运行中的服务器目录。
- 不把 `results/`、旧可执行文件或旧静态库整体复制成新项目。仅引入必要源码、配置、明确列出的回归数据和来源清单。
- 不从目录名、`build_hash` 子串或遗留 `diagnostic_sc_lod_campaign` 推断新算法行为。旧版确实有这种分支选择方式，新版必须改为有类型、可校验的显式配置。
- 不擅自改局部算子、AOT、归一化、标记排序、训练/压缩阈值或积分规则来获得加速。改数学方法须另列对照，不能混入迁移等价性验证。
- 当前已完成 P0–P6 实现及有界验收；实际范围见 `docs/migration_status.md`、`docs/p3_p4.md` 和 `docs/p5_p6.md`。P5 已通过三处暂停恢复、state 16 独立 fresh 审计及全样本短程审计。P6 已完成独立服务器 16 物理核同状态三次短基准。完整生产轨迹及长程对照仍按用户授权延期到 P8，不能记为已运行；P7–P9 尚待完成。
- 每阶段若遇到源码/论文不一致，先形成具体差异和最小重现；不能静默选一个、宣布完全复现。

完成标准：从干净的新项目检出，在 WSL/Linux 上独立编译；E1/E2、三类基线及 E1 四组标记对照可运行；能从保存状态进行独立诊断；第 11.3 节列出的当前论文图表可再生；支持真正的断点续跑；保留必要资源保护与可核查计时。在线维数比较、同预算误差比较和墙钟比较分别验收，不统称为计算成本上的公平优势。

## 2. 权威来源及必须解决的版本问题

以下路径相对旧项目根目录，论文路径除外。本次论文、图表输入、对照配置和代码入口的 SHA-256 见同目录 `ALOD_SOURCE_EVIDENCE_20260911.json`。`ALOD_SOURCE_EVIDENCE_20260909.json` 保留为历史清单，其中旧论文哈希不再代表最新版。两者用于定位及防止来源变化，**不是已经完成的完整依赖清单**。

| 来源 | 用途 | 注意事项 |
|---|---|---|
| `<PAPER_ROOT>/helmholtz_lod_certified_amsart_revised.tex` 及其当前 `tables/`、实际引用的 `figures/` | 第 4 节 residual/localization、第 5 节区域增广、第 6 节算法、第 7 节实验和附录 | 以 LaTeX label/实际 include 为准；目录中仍存在但正文不再引用的旧图表不是必交成果 |
| `results/E1-periodic-lazy12-51-20260909/data/config.json` 和该目录 `runs/` | E1 最终 51 状态配置、网格、lazy 检查与误差基准 | 有部分配置藏在当时驱动代码中，JSON 不足以独立还原 |
| `results/E1-production51-paper-20260909/` | E1 论文数值、图表和基线核对 | 拟合使用完整更新周期终点 |
| `results/E1-family-marking-controls-20260910/CONFIGURATION.md`、`configs/`、`data/`、`paper-update/summary.json` | E1 固定 ell=3 的四组标记对照，Figure 3 和同预算表 | 51/51/97/97 状态、每状态 48 RHS；不替换 Figure 1 的变 ell 生产轨迹或其 87 点 AFEM |
| 同目录 `reproducibility/production-source/`、`build-provenance.json`、`reproducibility/recovery/` | 标记对照真实生产源码、构建及恢复校验 | 普通 AFEM 用原生产程序完整审计；family 恢复程序仅用于归档重放，不作为新项目默认算法 |
| `results/E2-regional-adaptiveell-inherit-controls-20260909/deployment/source/`、`deployment/configs/inheritell.json` | E2 生产版完整源快照与配置 | 数学实现优先从这里提取；不是主工作区所有文件的相同版本 |
| `results/E2-regional-adaptiveell-inherit-controls-20260909/outputs/inheritell/` | E2 正式接受状态 | 使用 `accepted_control_*.csv`，根目录 `control_*.csv` 可能是升层前解 |
| `results/E2-AOT-factor-reuse-20260909/` | 最新已验证的 AOT 分解缓存、runner 和 `.inc` | 曾链接生产快照静态库，不是独立项目；必须提取其依赖并重新全量构建 |
| `results/E2-production-inherit-paper-20260909/` | E2 生产、三类基线和冻结空间审计；另存 75 次历史独立 RHS 审计 | 当前稿不再引用 shared/fresh 图及五样本表；历史数据保留，不能自动重新插入或当作必跑项；旧计时不含后来 AOT 缓存 |
| `results/E2-paper-revision-20260909/CHANGELOG.md` | 历史修订记录 | 不覆盖 2026-09-11 最新正文；精简项目仍不承担全部严格区间认证功能 |
| `results/E2-regional-budget-inherit-controls-20260909/analysis/REPORT.md` | 预算/继承成本与局限 | 旧全域预算和新区域预算不能混淆 |
| `results/E2-regional-adjacent-fast-20260909/README.md` | 已实现的并行、缓存、相邻层诊断优化 | 昂贵全域 corrector 探针已由相邻层诊断替换 |

已核查的具体版本差异：主工作区 `include/helmholtz/adaptive/kernel_residual.h` 的 `AdditiveKernelRieszContext` 只有基础接口；E2 生产快照同名头文件含 `apply_selected` 和 `selected_values`。仅复制主工作区 runner 加当前头文件，不能视作区域 AS 已完整迁移。

P0 须完成 `source_manifest.json`：每个引入文件的原路径、哈希、来源快照、许可证/版权说明、修改理由；包含 compiler、CMake、Eigen、SuiteSparse、BLAS、OpenMP 版本及编译/链接命令。禁止继续以 `git_commit=WORKTREE` 作为唯一出处。E1 如需服务器源码，从其已归档部署目录只读取回，并核对可执行文件及构建来源，不猜测缺失补丁。

旧 `docs/E2_regional_AS_AOT_control_plan_20260909.md` 是早期方案，含已废弃的全域 corrector 探针和早期半径。只能当历史材料，不能覆盖本计划和最终生产配置。

### 2.1 论文—实现对应表

P0 将下列 label、实现函数、配置字段和回归输入写入 `docs/paper_implementation_map.md`。理论定义、历史执行方式和新接口分别记录，不因接口整理暗改数学对象。

| 当前论文 label | 必须实现/核对的内容 |
|---|---|
| `eq:local-kernel-space`、`eq:as-local-riesz` | 全体几何粗顶点的 `W_h(N^2(z))` 和正定局部 Riesz；边界约束、patch 重数、粗质量重分配 |
| `eq:certifier`、`eq:localized-coarse-smallness`、`eq:reference-interval` | Theta 的算子范数定义及 Ritz 近似；可解性/可靠性所需常数与实用升层阈值区分 |
| `eq:regional-as`、`eq:mode-deflation` | 完整 patch 区域选择、无权 AS 和两遍能量消去；原始 kernel 代表与投影后 trial 模式区分 |
| `eq:aot-matrix`、`eq:coupled-block-system` | 全参考域逆伴随响应、保留基础 test、块系统恒等式及数值残差 |
| `eq:regional-training-target`、`eq:pod-dictionary`、`eq:distance-gate`、`eq:residual-gate` | 冻结 base 目标、归一化 correction 坐标 POD、逐成员实际耦合解压缩检查 |
| `alg:singular`、`alg:smooth` | 固定网格 solve/训练/升层阶段，以及 `m_ref` 轮伴随参考细化 |
| `fig:e1-family-marking`、`tab:e1-family-control` | 新增四组固定 ell 标记对照、同预算逐样本插值与参考精度 |

若代码采用两侧最大值近似 Theta，或局部子空间/边界顶点集与论文写法不直接一致，必须先给出映射或固定网格等价验证；不把函数名一致当作定义一致。

### P5/P6 执行记录（2026-09-12）

- P5：便携二进制检查点、原子提交、日志前缀校验和恢复已实现；第一参考 sweep 后、训练前、ell 升层后三处恢复与不中断轨迹一致。只增加 horizon/资源的覆盖允许；数学字段变化禁止覆盖旧运行，新的运行按配置/样本内容产生独立 experiment ID（不创建 Git 分支）。
- 独立审计保存 e/f/g、E/E_ref/G、参考 FEM floor、rank-zero、fresh、pure 及参数/解身份；完整 48 RHS 和两个 pure 的短状态批大小 1/8 对比通过。真实小网格 state 16 的单文件 fresh 审计通过，无 H0–15 重放。
- 64 MiB 旧 CSV 字段上限、外置带版本/哈希的标记索引、canonical 表、独立完成状态和后处理失败隔离已实现。完整生产 `paper_complete` 仍须满足正式样本、目标 horizon 和校验；小状态验收不能标记论文复现完成。
- P6：生产仅训练族及名义 RHS；载荷并行、融合误差积分、基础 RHS/增长字典/POD 响应及小 PG 矩阵缓存、ell-only 参考/Riesz 复用、已接受训练解复用、NVB 暖启动注入已实现，实际 AOT/PG 残差门限保持。参考/FEM、load/solution 和缓存的有名稠密分配与进程峰值分别记录。
- 已按 PowerShell/SSH、独立目录、tmux、16 不同物理核、BLAS 单线程和内存/交换区监控完成 H6/h10、ell=2、16 RHS 的三次同状态短基准；中位数/范围和源码/二进制哈希见 `docs/provenance/p56_validation.json`。raw 上三角 PG 和统一对称 forward/adjoint 因子仅为实测候选，未替换生产公式；不宣称 32 线程或完整轨迹收益。
- P7 继续负责更完整的小/中集成、资源故障/并发作业监督、可移植结果包重建、剩余按需 candidate/adjacent/reference 深化/物理区域诊断，以及更细的 phase timer、闭包/晚期分配归因。P6 现有阶段计时不覆盖所有内部子阶段，禁止嵌套相加。
- P8 保留所有既列的 E1 51/97 状态、E2 33 状态及 every/lazy/inherit/reset 长程对照、完整审计和论文图表。大网格 Theta/AOT、晚期 checkpoint 容量/吞吐及可选 32 线程对照仍需后续服务器实验，不以本次小状态数值一致或低 RSS 代替。

## 3. 保留范围与删除范围

| 模块/功能 | 新项目处理 |
|---|---|
| 二维三角形网格、嵌套加细、闭包、边界标签、层级/注入 | 保留 E1/E2 实际使用部分，保留稳定元素 ID 和父子关系 |
| 共形 P1 Helmholtz 装配、正定能量矩阵、混合边界、制造解及积分 | 保留；特别保留凹角和局部 Gaussian 的积分路径 |
| 准插值、局部 kernel 约束、双边 LOD corrector、PG 解 | 保留；只保留已用 direct-Schur 路线及必要内部依赖 |
| 正定局部 Riesz、批量 RHS、全域粗估计子、区域 AS | 核心功能；共用 patch 数据和局部分解 |
| 区域预算、Krylov、压缩字典继承、POD、全域 AOT | E2 核心功能，不能保留为临时诊断注入代码 |
| matrix-free Ritz 的 `Theta_ell`、E1 lazy、E2 每状态升层 | 保留实际用到的计算路径；不迁移整个认证驱动 |
| 强残差 reference 标记、全域 Dörfler、粗更新伴随 `m_ref` 轮参考细化 | 核心调度，取代旧 reference-only/balance/epoch 策略组合 |
| AFEM、UFEM、SLOD | 保留两种论文问题的最小基线驱动；这里 SLOD 指现有标准固定层数 LOD 基线，不引入其他同名算法 |
| E1 nominal/family 标记对照与逐状态 48 RHS 审计 | 必需；复用 ALOD/AFEM 和标记器，以显式成员 ID 选择控制变量，不复制四份 runner |
| 网格/状态快照、检查点、独立 RHS 审计、区域/层间/floor/rank 诊断 | 保留，按需运行，不能改变生产状态 |
| balance gate 控制、`C_rel` 驱动的 reference 切换、旧 reference-only 控制器 | 删除；不留默认关闭但仍影响分支的生产路径 |
| 旧 balance 指标和 candidate dual | 仅可选诊断模块，默认不计算；不能读回主控制器 |
| 平衡通量估计器、全域 kernel 逆增广、显式奇异模式增广、旧固定 rank 训练 | 不迁移生产实现；除非闭包中的必要内部函数，亦不得自动注册旧模式 |
| hp/DG、污染/逆不等式研究、patch GMRES/多重网格/多种 Schwarz 预条件器、旧 E0/E3 等实验 | 不纳入功能范围；先做编译依赖闭包再剥离，不能按文件名盲删 |
| MPFR/MPFI 严格区间认证平台 | 不要求；保留实数谱诊断和必要数值失败检测，不把实用指标标为严格认证 |
| 旧部署补丁生成器、几十种日期脚本、重复巨型 runner、历史数据全集 | 不进入发布项目；以一个结构化驱动、一个诊断工具、一个部署入口替代 |

基线已有数据时，支持验证后导入，生成论文不强制重跑；但新项目必须具备重新生成这三类基线的能力，才算能独立实现论文两个实验。

当前 E2 的 75 次 fresh RHS 审计、双层 corrector 余项估计，以及未来“多个独立奇异分量”的研究均不列为生产复跑任务。保留经济的 fresh、rank-zero、固定 ell/半径等诊断接口；不得借精简迁移扩充第三个论文问题。

## 4. 不允许在迁移中改变的数学语义

### 4.1 解、残差及估计子

统一记 `A_h=K_h-kappa^2 M_h-i*kappa*M_R` 为 Helmholtz 矩阵，`E_h=K_h+kappa^2 M_h` 为正定能量矩阵。论文内积对第一个变量线性：`a(v,w)=w* A_h v`，`b_kappa(v,w)=w* E_h v`，星号表示共轭转置；不要把正定 `b_kappa` 误写为不定 Helmholtz 形式。

准插值按论文为 `I_H=E_H Pi_H^dg`：单元 L2 投影后按单元面积加权节点平均，Dirichlet 节点置零。P0/P2 核对源实现及 `I_H v_H=v_H`。估计子节点集为所有几何粗顶点，局部 kernel 为 `W_{h,z}={w in ker I_H: supp(w) subset closure(N^2(z))}`；不得将它与自由粗未知量集合、corrector 的 `N^ell(T)` patch 混用。受约束子空间可用等价隐式表示，不要求显式组装稠密 kernel basis。

- `U_base`：独立 rank-zero LOD PG 解。
- `U`：E1 为 rank-zero 解；E2 为最终采用的增广耦合解。
- `u_h`：当前参考网格 FEM 解，仅在需要误差诊断时求解。
- `u_c`：指定诊断 candidate 网格 FEM 解，不是强残差标记的前提。
- `r_mu=f_mu-A_h U_mu`。
- 每个全域粗节点 patch 在其受约束 kernel 上解正定 Riesz 问题；`eta_H,mu^2=sum_z ||xi_mu,z||_kappa^2`。它不是一次全域 kernel Riesz 范数，也不是 AS 向量和的能量范数。
- 节点平方指标按所属粗单元数均分，保持总和不变。训练族用 `d_mu=max(||U_mu||_kappa,1e-12)` 归一化后做平方指标均值标记，再补充当前最坏归一化样本的 Dörfler 条件。一次决定内冻结 `d_mu`；不得使用精确解范数进行生产归一化。
- E2 粗标记用 `eta_H(U_coupled)`。`eta_H(U_base)` 仅用于必要的基础计算/可选对照，不替换生产分母。
- 参考标记用嵌入当前解的强残差，按 `eq:strong-reference-indicator`：`h_K^2 ||f+kappa^2 U||_K^2`，内部边每个相邻单元分配 `0.5*h_e*||jump(partial_n U)||_e^2`，Neumann 为 `h_e*||partial_n U||_e^2`，Robin 为 `h_e*||partial_n U-i*kappa U||_e^2`。P1 单元内 Laplacian 为零；jump 为两侧外法向导数之和。核对旧实现的等价计数，不多计一次边。该量仅用于标记，论文未为嵌入的非 Galerkin 解调用强残差可靠性定理。两实验都在全域标记；AS 的区域掩码不作用于 reference Dörfler。
- mean 加 worst-member 补集只保证均值和所选最坏成员的 bulk 条件，不保证所有成员同时满足；nominal-only 是单成员特例。稳定并列排序及节点到单元重分配必须保存，不能用不同的 family 聚合替换。

### 4.2 区域 AS、训练和压缩

`D_R=Omega ∩ B_R(0)`，生产 `R=0.6`。只有完整 patch 支撑包含于区域内才选中；不能用 patch 中心代替。空掩码返回零增广和清晰状态，不自动扩大 R 或回退全域。

定义 `M_D g=sum_selected I_z B_z^{-1} I_z^* g`，其中 `B_z` 是正定能量限制。首种子 `M_D r` 复用同一残差、同一状态 estimator 的局部解；若残差已改变必须重新 apply，不能复用旧种子。后续方向按 `eq:mode-deflation`：对已接受、归一化的当前方向 `phi_j` 计算 `M_D A_h phi_j`，再投影消去基础 trial、保留字典和本块已有方向。两遍能量 Gram–Schmidt 后做相关性判定及归一化。旧代码若写成 `A_W`，须证明坐标变换等价；不能直接使用未消去的原始 `psi_j` 代替 `phi_j`。不得重新引入全域 `B_W^{-1}`。

恒等式 `Re(r* M_D r)=eta_D^2` 应可核查；重叠子空间求和不做 partition-of-unity 平均。若合并相同 patch 求解，应保留完整/选中节点的分别重数，不能把重复 patch 直接去重成权重一。

原始 kernel 代表应满足区域支撑和准插值约束；相对基础 LOD 能量投影后的 `phi` 可以含全域 LOD 分量。不要以 `phi` 域外不为零判错，也不要承诺区域 AS 完全不影响 localization 误差。

E2 训练规范：

1. 求当前 mesh/ell 的基础训练族，固定 `t_mu=max(1e-12,0.1*eta_D(U_base,mu))`。
2. 先检查继承字典；若 `max_mu eta_D(U_mu)/t_mu<=1`，无需增长。
3. 否则选最大预算比样本，每轮最多新增 2 个方向；工作维数上限 24。线性相关方向丢弃并记录。
4. 连续两次最坏预算比相对下降小于 1% 时按旧规则退出；保存达标、停滞、rank 上限、空区域和数值失败的不同状态。不得用退出码正常替代“训练达标”。
5. 能量 POD 使用 `Dtilde_mu=(U_full,mu-U_base,mu)/d_mu` 的 enrichment 坐标 `S=Phi_full* E_h [Dtilde_mu]`；对 `S` 做 SVD，取 `Phi_r=Phi_full U_r`。非空字典从 rank 1 起检查最小合格 rank，逐候选重新求实际耦合解。每个训练样本满足 `||U_full-U_r||_kappa<=0.2*t_mu` 且 `eta_D(U_r)<=max(1.1*t_mu,1.03*eta_D(U_full))` 才接受；若无压缩候选合格保留完整字典，完整字典为空则 rank=0。不要把 POD 截断能量当作 PG 解误差，也不要额外套入旧 5% snapshot 门槛。
6. 工作维数、独立原始方向数、压缩在线维数分别输出；不能写成固定 16 维，也不能用 `16×4` 推断维数。实际生产在线为 1–2 维，但不能硬编码这个结果。

P4 新核对差异：旧 E2 冷启动在 family 贪心前先从 nominal 构造两列原始 kernel Krylov 向量，再消去基础 trial；第二列使用原始 kernel 方向而非规定的 phi。新实现严格采用上述空字典最大预算比和归一化 phi 递推，不恢复该旧启动分支。固定状态 base targets/Theta 严格回归，增广解单列数值邻近回归，不能声称完整 E2 旧轨迹数值等价；后续 P8 图表须标注这一差异。

独立单 RHS 从空字典重训只用于共享空间审计，不参与预算判据。训练族/测试族/shift/pure 的样本身份在接口和结果中明确区分。

压缩保证是已测试训练成员的**绝对** reference-relative 误差扰动不超过 `epsilon_dist*t_mu`；不是相对误差上界，也不保证未见载荷。训练退出、预算达标、压缩接受为三个独立字段：训练停滞后保留完整字典不能被标成区域预算已达标。

### 4.3 继承和 AOT

- 只改变 ell：保留同网格的 kernel 代表，更新 LOD 投影、全域 AOT 和当前层训练目标，再检查预算；不默认清空，也不复用旧投影后的 trial/test 对。
- 改变 H 或 h：嵌套注入旧 kernel 代表，采用生产版的 `M_D,new E_new P psi_old` 修复并校验新约束/支撑。此操作会改变方向，不声称精确保留旧空间。
- 继承对象为压缩字典对应的 kernel 代表。可在构造/压缩时同步保存代表及坐标，避免后来额外解 `C T` 恢复；这是需等价性验证的优化，不得直接把 `phi` 当 kernel 代表。
- AOT 保持全参考域的 `A_h^* z=E_h phi`，随后按旧法对基础 test 空间做能量投影和正交化，构造实际耦合 PG 对。**局部 trial 种子不意味着可把 test 方程局部化。**
- AOT 分解可跨 RHS 和仅 ell 变化复用；同一个分解不意味着不同右端的解可复用。改变参考算子/自由节点时失效。保留单条目或明确内存预算的缓存。

保留论文的可核查块结构：令基础 trial/test 矩阵为 `B,C`，`Phi*E_h B=0`、`Phi*E_h Phi=I`，则精确 raw AOT 满足 `Zraw* A_h B=0`、`Zraw* A_h Phi=I`；使用 raw tests 时耦合系统为 `[C* A_h B, C* A_h Phi; 0, I]`。投影/正交化 tests 只是换基，必须保留 `Y^0`。数值验收记录这些偏差、base-test residual 与实际 PG residual；近似 AOT 不能假定块恒等式严格成立。

论文只由基础块可逆推出增广对可解，其准最优常数使用自身 `alpha_r`，不能宣称 `alpha_r>=gamma_h`。trial 包含关系不保证实际 PG 误差单调下降。精简项目不计算全部认证常数，但不能在日志或说明中输出未经验证的认证结论。

### 4.4 `Theta_ell` 与升层

论文定义 `G_ell(v_H)(w)=conj(a(w,T^*_(ell,h) v_H))`，`Theta_ell=sup_{v_H!=0} eta(G_ell(v_H))/||v_H||_kappa`，其中 eta 使用第 4.1 节同一组局部能量 Riesz。分母是未校正粗函数的能量范数，不是 LOD 函数、载荷或相邻层解差。实插值和复对称算子下 primal/adjoint 缺陷范数相等；若源程序取两侧最大值，先检查该对称性和固定网格数值等价，不能相加两侧值。

保留 `reference_defect_spectrum_matrix_free` / `build_reference_corrector_certificate` 实际用到的谱核，迁入 `localization`。记录 Ritz 近似、残差、收敛状态和暖启动身份；单个 Ritz 值及很小迭代残差不自动给出 Theta 的严格上界。严格区间认证是范围外功能。

`Theta_ell` 不是 `||U_(ell+1)-U_ell||`。控制分母用名义 RHS 的未归一化 `eta_H`，不乘 `C_rel`。Ritz 容差生产为 `1e-4`，最多 750 迭代；保留暖启动及状态身份，不能把未收敛结果当正常数值。

- E1：初始/规定终点/距上次检查满 `2^ell` 个已求解状态才计算 Theta；比值 `>1.2` 优先升层，否则 `Theta<0.2` 接受，其他情况保持。中间参考 sweep 计入状态。迁移测试必须核对升层后重解及下一次计数的旧版语义。
- E2：每个最终训练压缩状态检查 `Theta/eta_H(U)>0.1`，固定 H/h 升一层，更新空间、目标和字典，重训/压缩并重检，直到通过或 ell=4。仅最终采用状态进入主曲线。
- 每次检查后重置 lazy 计数，而不是只在升层后重置。固定 ell 的 E1 对照即便保留 Theta 观察，也不进入升层分支。
- 零分母：分母为零而 Theta 正时触发至上限；Theta 为零则无需升层，二者为零记录 `zero_defect`、ratio=NA。图表中未计算/无意义比值记 NA，不能拿上次 Theta 冒充本步值。
- `ell=4` 未通过必须记录 `ell_cap_reached`；两策略是实用启发式，不输出未经验证的“已认证精度达标”。

`q_ell=Theta_ell/eta_H(U_nominal)` 对载荷幅值敏感：保持空间及线性求解器固定，将载荷乘 `c` 后分母乘 `|c|`，q 除以 `|c|`。区域训练/压缩的绝对 floor 还可能改变整条自适应轨迹，不能由此声称轨迹严格缩放不变。保留论文原阈值，不擅自加归一化；可选固定空间幅值诊断验证这个限制。该比值判据不替代理论条件 `C_a*delta_ell^2 < alpha_hat_(H,infinity)`，也不直接控制 reference floor。

## 5. 生产、基线和标记对照 preset

下表冻结当前论文生产，不用早期诊断配置覆盖。旧 JSON 中无效遗留字段不等于新 preset 需要实现。

| 参数 | E1 | E2 |
|---|---|---|
| 问题 | 单位方形局部振荡光滑制造解，case R1 | L 形域凹角 + 波包，case S |
| kappa / 初始层级 | 16 / H6,h10 | 16 / H6,h10 |
| 边界/精确解 | 原论文及 `paper_cases.cpp` | 原论文 `boundary-weight-gaussian-alpha80`，名义波振幅 0.5；不加显式奇异模式 |
| 训练/test/shift | 16/24/8；训练 seed 20260901，test seed 20260902 | 16/24/8，seed 20260831；另 2 个 pure RHS 仅诊断 |
| 粗 Dörfler | 0.15，训练族 | 0.15，耦合解训练族 |
| 参考 sweeps | `m_ref=2`，theta 序列 `[0.3,0.2]` | `m_ref=1`，theta 序列 `[0.2]` |
| 参考标记 | 嵌入 LOD 强残差，全域 | 嵌入耦合解强残差，全域 |
| 层级闭包 | 保留 local reference gap=2，active patch layers=0 的实际闭包语义 | 不启用旧物理半径 matching；不额外强加该局部 reserve |
| 增广 | 无，rank=0 | 区域 AS R=0.6，预算及继承见第 4 节 |
| ell | 从 2 到最多 4；lazy 1.2，tau=0.2 | 从 2 到最多 4；每状态比值阈值 0.1 |
| 完整生产 | 25 粗更新周期，51 已解状态 | 32 粗更新周期，33 已解状态 |
| 历史额外检查 | 状态 35、44 保留检查；最终 50 检查 | 升层前/旁路不增加正式状态数 |
| 积分 | triangle 12 / Gaussian 16 / singular 24 / recursive 8 | 同前三项 / recursive 6 |
| 局部 / 参考求解 | direct Schur / UMFPACK，复用相同 patch 分解 | 相同；启用已验证 AOT 分解缓存 |
| reference / candidate 上限 | 500000 / 600000 自由未知量 | 300000 / 400000 自由未知量 |
| 原生产线程 | 16 | 16；所有 BLAS 线程 1 |

精确制造解、源项及边界必须从旧函数提取并数值核对，不只重写论文公式；尤其是 E2 角点分支、极坐标角度范围、Gaussian 归一化、Robin 符号以及奇异积分。

问题注册表须明确：E1 `z0=(0.75,0.5)`、中心盘半径 0.05，y=0/1 Dirichlet、x=0 Neumann、x=1 Robin；E2 `Omega=(-1,1)^2 minus ([0,1]x[-1,0])`，凹角角度 `0<theta<3*pi/2`，两条凹口边 Dirichlet、其余外边 Robin、无 Neumann。E2 用固定形状 `u_sing=(1-x^2)^2(1-y^2)^2 r^(2/3) sin(2 theta/3)`，族参数 `c in [0.75,1.25]`、`|z-(-0.5,0.5)|<=0.05`、`a in [0.25,0.75]`、`phase in [-pi/4,pi/4]`；包络及相位平移严格按 `eq:e2-family`，shift 样本使用名义振幅和相位。制造解形状只用于构造 load 和审计，不作为增广种子。

将完整样本参数导出为固定 `rhs_e1.json`、`rhs_e2.json`，保存样本 ID、角色及全部参数。seed 仅作来源说明：不同随机库或调用顺序可能让同 seed 生成不同实验。两实验训练族均已包含名义成员；48 个正式族样本=16 train+24 test+8 shift，不能再加一个 nominal 变成 49 个。E2 另有 2 个 pure 诊断样本，不并入正式族。E1 四组对照的 sample 0 为 nominal，0–15/16–39/40–47 分别为三种角色，逐状态使用相同身份和坐标。

基线 preset 也要保存原有效配置：AFEM 使用名义 RHS、theta=0.15，论文终止误差目标 E1 为 0.008、E2 为 0.003；UFEM 做均匀加细；SLOD 固定 ell=3、粗/参考四级间隔。E2 的完整 UFEM 为 H6–H19，SLOD 为 H6–H15。其他终止层数、积分和族评估参数从已核对的基线 manifest 导入，不根据图上截点重建实验 horizon。精确误差可用于这两个 AFEM 基线的指定终止目标，仍不得进入标记或 ALOD 决策。

最新论文的完整基线状态数为 E1 AFEM/UFEM/SLOD=87/15/11，E2=88/14/10。E1 四组对照中的两个 97 状态 AFEM 使用固定更新数，不沿用主基线的 exact-error stop；必须是不同 preset/run ID，不能把 97 点覆盖进原 87 点表格。

### 5.1 两种可核查运行模式

- `paper-replay`：保留旧 E1 35/44 检查及旧审计安排；E2 可重现旧升层双分支诊断，但主状态只采用 inherit。用于核对历史轨迹。
- `production`：相同数学更新与参数，仅运行采用的分支；审计从检查点独立进行。旧终点检查只在显式 `extra_check_states` 中保留，不能隐藏在代码。若删除这些检查，标记为新轨迹，不声称逐点等同论文。

观察性 balance 计算的开关不应改变主解/网格/ell。删除它之后如轨迹改变，应定位隐式状态副作用，不能解释成算法本来如此。

### 5.2 E1 四组受控实验（当前正文必需）

拟实现 `configs/e1-marking-{nominal-alod,family-alod,nominal-afem,family-afem}.json`，由一份公共配置和可验证覆盖生成；解析后必须导出四份完整有效配置。唯一主动对照变量是参与标记的载荷集合，ALOD 同时改变粗标记和参考标记的成员，AFEM 改变其强残差标记成员。标记集合自然形成不同网格，不冻结共同终态网格。

| 参数 | nominal-ALOD | family-ALOD | nominal-AFEM | family-AFEM |
|---|---|---|---|---|
| 标记成员 | `[0]` | `[0,...,15]` | `[0]` | `[0,...,15]` |
| 初始 H/h | 6/10 | 6/10 | 6/无参考 | 6/无参考 |
| ell / rank | 固定 3 / 0 | 固定 3 / 0 | 不适用 | 不适用 |
| theta_H | 0.15 | 0.15 | 0.15 | 0.15 |
| m_ref / theta_h 序列 | 2 / `[0.3,0.2]` | 同左 | 不适用 | 不适用 |
| reference closure | local reserve=2、active layers=0 | 同左 | 不适用 | 不适用 |
| 固定更新周期 / 状态数 | 25 / 51 | 25 / 51 | 96 / 97 | 96 / 97 |
| 每状态评估 | 相同 48 RHS | 相同 48 RHS | 相同 48 RHS | 相同 48 RHS |

共同条件：E1 制造解与固定样本表、kappa=16，积分 12/16/24、递归 8，16 线程、BLAS=1；ALOD direct-Schur 与参考 UMFPACK，AFEM 批量 UMFPACK。关闭 balance 控制及默认 observer，不求 candidate FEM 解。ALOD 保留初末/每 8 状态的 Theta 观察但不许升层；不能把主生产的额外 state 35/44 检查移入此 preset。参考解及 exact error 仅供审计，不影响任何一组的标记或固定 horizon。

资源起点：ALOD reference/work 上限 600000/700000 自由未知量，AFEM 500000 未知量及 1200000 粗单元，墙钟上限 14 天；已有服务器每组 RSS 限额 ALOD 125 GiB、AFEM 24 GiB。新部署仍按第 12 节实测留量，不能用主生产 E1 的 500000/600000 上限悄悄替代。JSON 中旧 gate/epoch/物理半径字段不进入新接口。

验收包含：四组样本完全一致；初态两组 ALOD 的解/误差/参考精度一致；单组所有状态恰有 48 条审计记录；16/1 标记成员数正确；mean/worst bulk 成立；至少一个标记集合确实不同；正式数据行总数 `48*(51+51+97+97)=14208`。其余 15 个训练载荷在 nominal 组只评估，不参与网格或停止选择。

既有数据已经完成，可验证后导入生成图表；新项目须通过名义/家族两种模式的短轨迹数值回归并支持完整 horizon 重跑。不要为“生成论文图”默认再占用服务器重跑四组。

### 5.3 配置与曲线用途分离

`experiment_id` 区分 E1 变 ell 生产、E1 固定 ell 标记对照和 E2 区域 AS 生产；`method` 与 `marking.member_ids` 正交。`stopping` 显式选择固定 cycles 或仅用于基线的 exact-error target；`audit` 选择应采集的状态/样本；`plot` 仅控制展示、拟合和插值，不反向影响求解。

所有比较从共同固定样本表取值。E1 主生产 51 状态和两个固定 ell 的 51 状态不是同一轨迹；旧文稿中的最终 2422 在线 DoF 不得填到新对照的 2522/2878 两行。标记收益应按共同在线预算评价，不能用不同完整终点的误差直接计算 family 收益。

### 5.4 本次 P3/P4 验收与延期清单（2026-09-11）

用户已确定：E2 lazy ell check 只复用 E1 的检查时机，保留 E2 阈值 `Theta/eta_H > 0.1`。独立保留 `every` 模式。lazy 在初始、终点、额外指定状态及距上次检查 `2^ell` 个已接受状态时检查；每次检查均更新时钟。触发 E2 升层后仍在同一状态内训练/压缩/复查到接受或 ell=4。未检查状态的 Theta 为 null，不能沿用旧值冒充新值。

P3/P4 本次完成实现和有界验收，包括 m_ref=1/2/3 完整/中断周期、E1 升层附近六状态旧版对照、四组 48 样本短程契约、区域恒等式/POD/继承/AOT 失效检查及 E2 同网格强制 2→3→4。完整轨迹不作为本次对话必须启动的长作业，以下项目明确延期，不能遗漏或以旧导入结果替代新结果：

| 延期项目 | 后续阶段及验收 |
|---|---|
| E1 主生产 25 周期/51 状态，额外检查 35/44 | P8：manuscript 模式完整曲线和事件 |
| E1 nominal/family ALOD 各 51 状态 | P8：每状态同 48 样本审计，共同在线预算比较 |
| E1 nominal/family AFEM 各 97 状态 | P8：每状态 48 样本，固定 horizon |
| E2 every/inherit 32 周期/33 状态 | P8：按新规范递推的完整生产轨迹，标注旧冷启动差异 |
| E2 lazy/inherit 33 状态 | P8：与 every 比较误差、ell 事件及实际耗时 |
| E2 lazy/every 的 reset 对照 | P8：继承收益及共同状态/预算比较 |
| reference FEM floor、fresh、pure RHS 的独立审计 | P5 已实现并完成短程验收；P8 保留完整轨迹审计 |
| 大网格 AOT/Theta、粗暖启动注入、内存及可选 32 线程 profile | P6/P7：实际测量，不提前宣称收益或完整可扩展性 |

所有新说明文档使用英文，本计划保留中文。full presets 仅保存有效配置，不表示已运行或当前存储上限可容纳所有后期状态。checkpoint/resume 已在 P5 完成；P6 完成有界性能测量，大规模资源监督及吞吐量仍需 P7/P8。后续长任务继续遵守本计划的 PowerShell/SSH、独立目录、tmux、16 物理核、BLAS 单线程及内存/交换区监控规则。

## 6. 新的统一参考细化调度

配置采用 `reference.sweeps_per_coarse_update=m_ref` 和长度恰为 `m_ref` 的 `theta` 数组。可接受单个 theta 并显式展开为数组；输出解析后的有效配置。`m_ref` 为正整数，不是加细层数，不是单元内部递归次数。

建议状态机：

```text
初始化 H、h -> solve_accept(初始状态) -> 保存状态
for coarse_cycle in 1..max_coarse_updates:
    用当前最终接受族解同时形成粗标记和第一轮参考标记
    对两套网格加细并完成相容/嵌套/指定 gap 闭包
    将临时 work mesh 提升为新参考 mesh
    solve_accept(cycle, sweep=1) -> 保存状态
    for sweep in 2..m_ref:
        保持当前粗网格
        用上一轮新接受解重新计算强残差及参考标记
        加细参考网格并完成闭包、提升
        solve_accept(cycle, sweep) -> 保存状态
```

闭包实现须保留新核对出的旧 E1 细节：reference-only sweep 仍由当前新解产生候选粗标记，在候选粗网格上完成相容及局部 gap 闭包后丢弃粗网格提案，仅保留新参考网格。活跃区包含 NVB 一致性补细的全部单元并扩展 ell+1 圈；不能简单省掉候选粗闭包，否则短轨迹会从 state 2 起改变。

`solve_accept` 内完成对应 E1/E2 的空间更新、训练、压缩和 ell 检查。E1 在第二 sweep 时可升 ell，但不增加粗网格。E2 在 ell 循环中不能提前提交升层前解。

独立计数：`coarse_cycle`、`reference_sweep`、`state_id`、`ell_event_id`、`mesh_revision`。初始状态为 state 0；完整 K 周期应有 `1+K*m_ref` 个接受状态；ell 子迭代/诊断不计额外正式状态。可在中间 sweep 达到资源/用户上限，记录未完成周期及下一阶段，不强行凑满。

验收必须包括 `m_ref=1,2,3`；`m_ref=3` 只是功能测试，不能写成已验证其科学效果优于论文配置。E1 的 `[0.3,0.2]` 不得简化成一次 theta=0.5 或连续两次标记旧解。

## 7. 项目结构与提取边界

建议保持少量库和两个 C++ 入口，不按旧 benchmark 数量建立几十个可执行文件：

```text
ALOD/
  CMakeLists.txt
  README.md
  include/alod/                 # 小而有类型的公共接口
  src/mesh/                    # 网格、嵌套、闭包、准插值
  src/fem/                     # P1、边界、能量、积分
  src/lod/                     # 双边 corrector、PG、localization
  src/as/                      # patch Riesz、区域算子、训练、POD、AOT
  src/adapt/                   # family 标记、m_ref、ell 状态机
  src/problems/                # 仅 E1/E2 及 pure 分量
  src/io/                      # 检查点、结果、网格导出、provenance
  src/diagnostics/             # 按需启用，不能控制 adapt
  apps/alod_run.cpp            # E1/E2/AFEM/UFEM/SLOD
  apps/alod_audit.cpp          # 指定 checkpoint 的独立诊断
  configs/                     # E1/E2、三类基线、四组标记对照、smoke、诊断模板
  data/rhs/                    # 固定样本参数
  tests/fixtures/              # 小规模原版数值和矩阵/网格
  scripts/run_server.py        # 多作业资源和核心分配
  tools/analyze.py             # 数据校验、曲线、拟合、表格、区域标注
  docs/                        # 算法映射、实验协议、迁移记录
```

源码提取索引：

| 新模块 | 旧源码起点 | 拆分要求 |
|---|---|---|
| mesh / hierarchy | `src/mesh/{refine,edges}.cpp`、`src/lod/quasi_interp.cpp`、`src/helmholtz/adaptive/hierarchy.cpp` | 保存层级、节点顺序和闭包；检查真实传递依赖 |
| fem / problems | `src/helmholtz/{operators,boundary,quadrature,manufactured}.cpp`、`src/helmholtz/benchmarks/paper_cases.cpp` | 只保留两问题；源项/精确解/积分分别可复用 |
| lod | `model.cpp`、`corrector.cpp`、`patch_system.cpp`、`patch_solver.cpp`、`corrector_pipeline.cpp`、`src/lod/patches.cpp` | 不将所有 hp/Schwarz 路线拖入最终构建；必要通用函数单独提取 |
| estimator / regional AS | E2 快照 `adaptive/kernel_residual.{h,cpp}` | 分离 patch 几何、局部分解、apply/estimate、区域 mask 和聚合 |
| localization | `adaptive/certificates.cpp` 的 reference defect/Ritz 路径 | 保留实际用到的数值谱核；隔离证书注册、区间数、全认证驱动 |
| training / AOT | 优化目录 `reference_epoch_paper_runner.cpp`、`E2_region_controls_helper.inc`、相关 `.inc` | 把成员注入和日志副作用重构为正常类；提取所有传递 include |
| scheduler / marking | `reference_epoch_paper_runner.cpp`、`adaptive/reference_epoch_driver.cpp` | 重新组织状态机，保留标记数学；不整体移入旧 epoch 驱动 |
| baselines | `bench_helmholtz_adaptive_paper.cpp` 及其实际求解路径 | 共用问题、FEM、网格、误差和输出模块 |
| diagnostics / plots | 上述生产快照及 `tools/analysis/prepare_E2_inherit_paper.py`、`prepare_E1_production45_paper.py`、`prepare_E1_family_marking_paper.py` | 去除日期路径和硬编码数据根；按当前 LaTeX 引用筛选成果，去掉旧 E2 fresh 默认叠图；保留真实数据语义和拟合规则 |

接口原则：`Problem` 提供 load 与可选 exact evaluator；自适应代码不获取 exact error、`u_h` 或 `u_c`。共享 `MeshState/OperatorContext` 只读对象；每个求解/缓存对象明确所属 mesh、算子和空间版本。诊断通过不可变 `AcceptedState` 输入工作。

其中 AFEM 的论文终止检查放在基线驱动的外层，不把 exact evaluator 下传标记器。E2 所需 Base LOD 解属于训练算法，不能为减少“审计成本”而一并移除。

缓存生命周期至少覆盖下表；缓存条目要有统计和内存上限：

| 变化 | 可复用 | 必须更新/重检 |
|---|---|---|
| 仅 RHS | 网格、算子、corrector、局部 Riesz 和 AOT 分解 | load、解、残差、目标；exact/floor 缓存按 RHS 区分 |
| 仅 ell | reference 算子、参考 FEM 解、AOT 分解、同 H/h 的 Riesz patch context | LOD trial/test、投影、Theta、Base 解、预算、耦合小系统；kernel 代表重投影 |
| 仅 reference h | 在不变局部几何上经身份验证的局部对象 | reference 算子/自由节点、Riesz context、AOT、corrector、注入/字典修复、所有网格相关指标 |
| 粗 H 改变 | 仅经过依赖检查确认未受影响的对象 | 准插值/kernel、patch、corrector、字典修复与投影、所有控制量；不能只比较 Nref |
| kappa、边界或矩阵系数改变 | 纯几何对象在条件满足时 | 所有相关数值分解；同维数、同 mesh ID 也不代表同算子 |

开始用已验证的实际稀疏矩阵比对；未来要以版本 ID 降低 O(nnz) 检查成本，必须先证明所有算子修改入口都会生成新身份，不能只把哈希检查删掉。

构建依赖以 C++20、CMake、Eigen、OpenMP、SuiteSparse 为起点；CHOLMOD/UMFPACK 按实际闭包保留，别把 Eigen fallback 默认为生产等价后端。消除未用 TBB、MPFR/MPFI 等依赖。Python 分析只需明确列出的 numpy/matplotlib 等小集合。不要因拆文件引入额外 web/UI/服务框架。

## 8. 结果、诊断和检查点协议

### 8.1 总是记录的轻量数据

`manifest.json` 保存解析后配置、来源/新源码 commit、脏补丁哈希、二进制哈希、后端/依赖版本、线程/CPU 集、样本哈希、运行状态及结束原因。

每接受状态 `states.csv` 保存计数、mesh/字典哈希、ell、`N_H/N_reference/N_work/N_online`、标记数和闭包增加量、工作/在线 rank、训练停止原因、预算最大比、压缩检查、PG/AS/AOT 残差摘要。`N_online=N_H+r_online`；E1 rank=0。DoF 用处理边界后的自由未知量，另存节点/单元数。

另设 `samples.csv`、`marking.csv` 和 `plots.json`：逐状态/样本保存角色、nominal/marking 标志、参数身份、解/网格身份及精确范数；保存所选 worst member、归一化尺度、mean/worst 总量和捕获量；保存每幅图实际选用的状态、阈值、拟合窗口、插值端点/权重及源哈希。E1 四组完整输出需覆盖每状态 48 RHS，可离线批量补齐，但在补齐前不可标为 `paper_complete`。

误差字段明确对应：原始数据 `E=e/||u_exact||`，`F=f/||u_exact||`，`G=g/||u_exact||`；最新论文把 `F` 写作 `E_ref`。导入时保留原列并记录映射，不能把 `F` 当载荷泛函或更换分母。AFEM 的独立 reference error/gap 为 NA；遗留 `N_reference=0` 表示不适用，不表示 reference floor 为零。

`ell_checks.csv` 独立存 Theta、当前分母、比值、计算的状态身份、检查原因、触发/达到上限及谱收敛证据。非检查状态 Theta 为 NA；若展示缓存值，必须附 `computed_state_id`、网格/ell 身份并明确不是当前检查。

`timing.csv` 用唯一 event ID、parent ID 及 phase 记录。至少拆 mesh/closure、装配、corrector、Ritz、base solve、AS prepare/apply、训练评估、POD、AOT factor/solve/project、耦合求解、载荷、误差积分、reference solve、I/O。不能把嵌套 timer 相加作为总墙钟，也不能把一次共享构建在两分支重复计费。

### 8.2 按需诊断菜单

| 诊断 | 保存的量 | 默认策略 |
|---|---|---|
| exact / floor / gap | `e=\|\|u_exact-U\|\|`、`f=\|\|u_exact-u_h\|\|`、`g=\|\|u_h-U\|\|`，绝对量及同 RHS 精确范数归一化 | 生产不依赖；论文模式名义逐状态、全族指定状态或离线补齐 |
| E1 marking control | 四组每状态 48 RHS 的 E/E_ref/gap、标记成员与集合、mean/worst bulk、同预算配对误差 | 当前论文必需数据；审计可离线分批，不能用只有终态的旧 AFEM family CSV 补造逐状态曲线 |
| reference/candidate | `\|\|u_h-U\|\|/\|\|u_c-U\|\|`、candidate 网格构造身份及 DoF | 关闭；显式选择 checkpoint 和工作网格定义后计算 |
| 旧 balance 对照 | 逐 RHS `eta_H/eta_dual^c`、训练族 max/p90/p50/min，必要时与真实 gap 比相除 | 关闭；沿用旧局部 candidate Riesz 定义及数据口径，不乘 C_rel，不加入控制分支 |
| localization | `Theta/eta_H(coupled)`、`Theta/eta_H(base)`、相邻层 `\|\|U_next-U\|\|/g`、next/current error | 相邻层仅指定状态；生产 Theta 仍按各自必需频率 |
| 区域误差 | D 内/外 base 与 coupled gap，pure-corner/pure-wave，区域支撑、kernel 约束、AS 恒等式 | 生产保留廉价约束；昂贵积分指定状态；圆穿单元的积分不能用重心归类替代 |
| 训练/压缩/rank | 每增长步最大预算比、最坏样本、独立方向数、工作/在线 rank、压缩前后距离与区域残差 | 主训练已有量直接记录；固定 rank/压缩关闭只作独立对照 |
| 字典继承 | incoming/repaired/dropped rank、kernel/支撑偏差、修复耗时，ell-only 与 mesh-change 分开 | 轻量必记；从空重训分支按需 |
| shared/fresh RHS | 相同 H/h/ell 下单 RHS 空字典训练；shared/fresh rank、预算、误差、signed percent | 可选，默认关闭；历史 H=0/16/32 × (名义+24 test)=75 次已完成，可验证导入；不要求再跑或重新插入论文 |
| E2 same-state rank-zero/AS | 同一 H/h/ell 和 RHS 下的 U_base/U_coupled、区域内外误差、eta_D/eta_H、成本 | 保留经济的可选诊断，优先复用训练已有 base 解；不能用不同自适应网格的 rank-zero 轨迹宣称隔离了增广收益 |
| 载荷缩放 / 双层 corrector | 固定空间下 q 的缩放；指定 L>ell 的 Delta_(ell,L) 及可选余项 | 仅诊断模板；双层算子范数不同于单 RHS 解差，最新版明确未评估该双层诊断；不要求为复现生产求全域 ideal corrector |
| reference 深化诊断 | 固定 H/ell、逐轮 theta、实际 solve 后 f/g/e、DoF/时间/内存 | 独立 audit 模式，最多轮数显式给定；不返回生产调度 |
| fixed-ell / regional 半径敏感性 | ell=2/3/4、R、各自轨迹及同 checkpoint 对照 | 小型配置模板，不保留旧全套活动脚本 |
| 内存/性能 | 进程组 RSS、MemAvailable、线程 CPU、分解命中/重建次数、scratch/cache 估算 | 资源监测必需；详细 profile 可选 |

全域 AS 可作为“全选 patch”的低维护诊断特例；不需保留旧全域预算训练框架或全域 kernel 逆。不实施昂贵全域 corrector 作为常规探针。相邻层冻结 kernel 代表、重新投影并构造 AOT，不重训的分支与“升层后重新达标”的分支分开命名。已有逐成员残差、区域和 base 解数据应先复用，不为每个诊断再建一个大求解器。

任一诊断都不得修改生产字典、计数、Ritz 暖启动或标记缓存。最好在独立进程从 checkpoint 运行；同进程审计需证明结束后主状态哈希和后续输出一致。

比值分母接近数值精度时 NA 并说明原因。floor/e 与 gap/e 是范数比，不可相加；`eta_D^2/eta_H^2` 是重叠 patch 指标比，不是物理区域误差占比。相邻层差是敏感性，不是真正 localization 分量，也不保证升层后的误差下降。

### 8.3 真正的断点续跑

之前 E1 45→51 状态因无可恢复状态而重放旧轨迹；审计也通过回放到 H0/H16/H32 开展。新项目将此作为明确需要解决的成本问题。

检查点至少包含：

- 粗/参考网格的节点、单元、边界、父子树、稳定 ID、加细层、闭包/注入所需状态；只存末点 VTU 不足以可靠续跑。
- 接受状态的完整计数、下一 sweep/阶段、ell、上次 lazy 检查位置、额外检查表和 Ritz 暖启动信息。
- 最终原始 kernel 代表、压缩坐标/在线维数、当前接受解或可确定性重建数据、样本表和必要随机状态。
- 配置、算子、样本、网格、字典及格式版本哈希，结束/暂停原因、日志已提交行号。

数值 LU 因子可以不落盘，恢复时重建；第一版优先可移植和正确性。可选空间快照用于减少离线审计的 corrector 重建；明确存储成本和版本校验。

采用临时文件写完、校验、原子提交 checkpoint，再提交接受状态日志；异常写入不得覆盖上一个有效状态。恢复不能重复写一行接受状态或重做已经提交的参考 sweep。支持只增加 horizon/资源上限的受控配置覆盖，数学参数变化产生新分支 ID。

验收：不中断轨迹 vs 在第一参考 sweep 后、ell 升层后、训练边界暂停恢复的最终网格/ell/rank/误差一致；至少能从 H16 checkpoint 独立跑 fresh RHS，而不重放 H0–15。

### 8.4 完整性、旧数据导入和失败隔离

E1 四组最近暴露的失败发生在求解结束后的解析阶段：`element_ids` 超出 Python CSV 默认 131072 字符，验证器报错后队列中断另外两组。这不是 OOM。新格式将长元素列表存入有版本/校验的独立索引文件，CSV 保存集合 ID、数量和哈希；旧 CSV 导入器明确设定至少已验证的 64 MiB 字段上限，并给出超限错误，不能截断读取。

独立记录 `solver_completed`、`validation_passed`、`audit_complete`、`paper_complete`，只有目标 horizon、记录完整、校验通过才汇总为成功。一个作业的后处理失败不默认杀死其他已授权、资源正常的作业；共享资源 watchdog 与数据校验失败采用不同退出原因。审计失败从已提交 checkpoint 补做，不重跑网格自适应前缀。

导入旧四组时，科学主数据只认 `data/<arm>/family_control_samples.csv` 和完成/哈希清单。family-ALOD 的 states 0–48 来自原运行，49–50 来自经验证恢复；family-AFEM 为原 0–95 加恢复 96。其恢复 `runs/*` 内旧格式审计前缀可能含跳过记录，不能覆盖合并主表。原标记集合与前缀误差一致的校验必须保留；新项目不继承这种跳过审计的生产分支。普通 AFEM 的 97 点为一次完整新审计，前 87 点与旧主基线 DoF 一致、名义误差最大相对差约 `5.06e-13`。

可移植结果包必须包含有效配置、样本表、完成标记、canonical 数据、实际图表输入和脚本依赖清单；在一个新的本地目录解包再跑验证/绘图，防止缺 `all-controls-completed.json` 等入口文件。论文标签到输出文件由单独清单管理，不复用只能匹配旧正文段落的日期式论文替换脚本。

## 9. 现有性能证据和优化优先级

### 9.1 已实现，必须带入而不是重复立项

| 优化 | 已有证据 | 迁移要求 |
|---|---|---|
| 局部 AS 分解复用、批量 RHS、种子复用 | `AdditiveKernelRieszContext` 及区域生产快照 | estimator/AS 使用同一 context；保留重叠权重和局部残差校验 |
| RHS 载荷、精确误差/区域积分 OpenMP 并行，积分缓存 thread_local | `E2-regional-adjacent-fast` 说明及源码 | 不退回串行；输出按样本 ID 固定顺序；检查首次初始化竞态 |
| 同状态 reference 解/精确范数/floor 缓存；相邻层构建在诊断间复用 | 同上 | 当前层/下一层、压缩前/后不重复构造；缓存键区分算子和 RHS |
| AOT 分解跨 RHS、仅 ell 改变复用 | `E2-AOT-factor-reuse/validation.json`：26 次调用、1 次分解、25 次命中，数值差约 2.45e-14；ell/mesh 对照约 4.21e-14 | 保留 exact free-node/稀疏结构/系数检查及变化时失效；单条目，控制驻留内存 |
| 核心分配和 BLAS 单线程 | 多次部署已避免进程只继承单核 affinity | 在 OpenMP 初始化前分配完整多核 CPU 集；检测每个工作线程 affinity |
| E1 四组全 RHS 审计的批量参考/AFEM 求解和 RHS 并行积分 | `E1-family-marking-controls-20260910/reproducibility/production-source/`；普通 AFEM 全 97 状态完成 | 提取共享算子/因子与样本顺序保证；不要为新增 family-AFEM 再复制整套基线驱动 |

AOT 缓存仅做过数值/小规模命中验证；小网格分解约 0.005 秒、命中检查约 0.00034 秒，不是大规模总训练加速比。不能把已实现的优化重新宣传为本次新结果。

### 9.2 待优化，按实际收益排序

E2 inherit 33 状态历史运行：总墙钟 58 分 08.5 秒；所有 family 训练约 105.14 秒（约 3%），误差审计约 2076.22 秒（约 59%），载荷装配约 346.94 秒（约 10%），峰值 RSS 31.14 GiB、平均 CPU 1076%（约 10.76 个逻辑核的 CPU 时间占用）。该运行包含旁路诊断和旧 AOT 计时；这些比例有重叠阶段口径限制，不能直接相加拆成完整账本。

新增 E1 对照的完整 nominal-ALOD 全审计墙钟为 31 分 50.33 秒、平均 CPU 1227%、峰值 RSS 81.52 GiB；nominal-AFEM 为 23 分 38.68 秒、1217%、1.425 GiB。它们是不同方法/网格/工作量，不能据此计算加速比。两个 family 恢复作业跳过已有前缀的重复审计，其墙钟不能与上述完整全审计计时直接比较。ALOD 81.52 GiB 的具体分配来源还需测量，不在未 profile 前归因给某个单独缓存。

| 优先级 | 实施内容 | 怎样验证；不能承诺什么 |
|---|---|---|
| P1 | 算法与审计解耦；训练只处理必要训练族，test/shift/pure exact 和 FEM 审计离线；默认不算 candidate 解/dual | 同配置审计开关不改变主轨迹；单独报告 production wall 和 audit wall。不以减少论文数据为加速手段，所需图表数据离线补齐 |
| P1 | 检查点直达审计和续跑，去掉到指定 H 的重复重放 | 从同状态重建/恢复与旧回放解一致；同时报告磁盘容量及恢复时间 |
| P1 | E1 四组逐状态审计有界批处理，识别 `Nref×RHS` 的 load/solution/residual 副本和 patch 因子峰值；长标记集合外置，流式校验 | 保存完整 14208 行语义及确定性标记；报告审计批大小、实际分配峰值和 wall，不能只凭 RSS 猜内存来源或删掉 held-out 点 |
| P1 | 载荷、积分几何、exact norm/floor/reference solve 按 mesh/算子/RHS 复用；RHS 分批处理 | 变化 RHS 不得误用旧解；从快照实际 profile 识别残留重复工作；线程/批大小受内存约束 |
| P1 | 生产升层只执行采用的继承分支；双分支比较迁入 audit | accepted 数值及后续轨迹等价；旧双分支事件成本不重复计费 |
| P2 | 固定 AOT/PG 矩阵做多 RHS block solve；缓存基础族解、E*Phi/A*Phi 和小矩阵，降低每个压缩 rank 的全量重装配 | 保留真正的每 RHS 解和压缩后距离/区域预算检查；确认所用 solver 的 solve 线程安全，不并发共享不安全 factor 对象 |
| P2 | 评估用 raw AOT 的块上三角系统复用基础分解：先求增广系数，再解基础块；POD 线性组合复用响应 | 候选优化，先验证 `eq:coupled-block-system` 的数值偏差、投影 test 换基和完整 PG 解等价；近似 AOT 不可强行将非零块置零，未经验证不替换生产路线 |
| P2 | 直接保存 kernel 代表和压缩坐标，避免每次从 phi 解 `CT` 恢复 | 小/中状态验证与旧代表张成空间、PG 解、继承后预算一致；不得跳过新 mesh 修复 |
| P2 | 对训练迭代只计算所需区域残差/种子；最终接受时计算全域 eta_H 供粗标记 | 目标及最坏样本序列与旧实现一致；最终全域粗标记不可被区域指标替代 |
| P2 | patch 缓存与每线程临时数组设内存预算，限制 Nref×RHS×多分支的稠密副本 | 测峰值 RSS，不能只测最终 RSS；分批 RHS 时保持相同数学聚合和稳定排序 |
| P3 | 测量 Theta 的 patch-Riesz apply、Ritz matvec、暖启动，减少不必要重复算子构造 | 保留生产阈值、频率和 1e-4 容差；保留 E2 every-state 独立模式；新增用户指定的 lazy 模式仅改变检查时机，阈值仍为 Theta/eta_H>0.1；不得把旧值当新值 |
| P3 | 评估统一 forward/adjoint LU 服务、符号分解复用、内存受限因子缓存 | 须先核实后端对共轭转置求解的支持及数值等价；作为候选，不保证当前接口已实现 |

全域 AOT 的大稀疏分解仍可能在更细参考网格上成为耗时和内存热点；缓存只消除重复分解，不消除一次分解本身。局部 AOT、迭代近似 AOT、改变区域 R 或预算不作为本次等价迁移优化。

性能验收采用同 H/h/ell/字典输入/诊断强度/线程和同后端的固定状态，在独立资源窗口比较；报告冷缓存/热缓存、至少 3 次短基准的中位数和范围。不同自适应终点不能当严格加速比。16→32 线程是否更快需单独测试，不能从服务器核多直接推断。

## 10. 分阶段任务与验收

下列命令/文件名是**拟实现接口**，不是当前已有可调用程序。先实现命令再在状态记录中填写实测命令和结果。

| 阶段 | 依赖 | 具体交付物 | 验收后才能继续 |
|---|---|---|---|
| P0 来源冻结 | 无 | source manifest、当前论文 label/include 清单、E1/E2 生产及四组对照有效配置/样本、源快照差异、最小 golden 输出、依赖闭包 | 能指出每个算法来自哪份代码；列出 JSON 外的隐藏分支；主生产与标记对照不混用；旧版小状态可重现 |
| P1 独立最小核 | P0 | CMake、mesh/FEM/边界/积分、两问题、三基线最小驱动、单元回归 | 新目录全量构建且不链接旧库；制造解/边界/网格及基线 smoke 通过 |
| P2 LOD 与 estimator | P1 | 双边 LOD/PG、批量 Riesz、显式成员 family/nominal 标记、强残差、Theta | 固定小网格与旧版解/指标一致；核约束、边界几何顶点和逐单元标记质量守恒；AFEM 与 ALOD 指标不混用 |
| P3 E1 调度与对照 | P2 | `m_ref` 状态机、lazy/固定 ell 策略、E1 主生产及四个对照 preset、事件输出 | m_ref=1/2/3 的完整/中断周期正确；E1 升层附近短轨迹与旧版对齐；四组相同样本、不同标记选择及 fixed horizon 正确 |
| P4 E2 区域方法 | P2、P3 | 区域 AS、预算、POD、继承、AOT 缓存、E2 ell 策略 | 空/全/部分 patch、AS 恒等式、强制 2→3→4、继承/reset 对照及 E2 固定状态数值回归 |
| P5 状态与独立审计 | P3、P4 | checkpoint/resume、audit 入口、全状态载荷审计、长字段兼容导入、可选 fresh/rank-zero 诊断 | 至少 3 个中断位置一致；审计开关不改主轨迹；direct H16 audit 不重放；独立校验失败不误杀其他作业 |
| P6 性能整理 | P4、P5 | 已实现优化完整迁入、P1/P2 待优化项的实际 profile 与实现 | 同网格数值等价；核心/内存正确；报告测得收益，无收益的项如实记录 |
| P7 小/中集成 | P6 | E1/E2 integration、基线及所有必需诊断回归、部署 smoke | 无旧库/旧绝对路径依赖；缺失数据、资源结束、失败退出均能正确解析 |
| P8 论文复现 | P7 | E1/E2 主生产 51/33 状态复现或明确差异，四组 51/51/97/97 及三基线验证导入/按需重跑，当前 7 图/8 表 | 四组 14208 行完整；共同预算与展示截点可追溯；无默认 E2 fresh 重跑或旧图叠加；不可用预填旧表冒充新结果 |
| P9 发布整理 | P8 | README、算法/配置/诊断/资源说明、来源与性能报告、可移植包 | 换一个干净路径构建与 smoke；仅含本范围功能；所有未达事项明确列出 |

每阶段更新 `migration_status.md`，记录完成/待办、命令、测试输入与哈希、数值最大差、资源和下一步。可自然切成逐阶段 agent 任务；不要并行修改同一 giant runner，模块边界稳定后再考虑独立任务。

建议用户执行入口：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DALOD_BUILD_TESTS=ON
cmake --build build -j <build_jobs>
ctest --test-dir build --output-on-failure -L smoke
build/alod_run --config configs/e1-paper.json --output results/e1
build/alod_run --config configs/e2-paper.json --output results/e2
build/alod_run --config configs/e1-marking-family-alod.json --output results/e1-family-alod
build/alod_run --config configs/e1-marking-nominal-afem.json --output results/e1-nominal-afem
build/alod_run --resume results/e2/checkpoints/latest --max-coarse-updates 40
build/alod_audit --checkpoint <saved-state> --config configs/audit-frozen-family.json
# 可选历史诊断，不是当前论文必跑项：
build/alod_audit --checkpoint <H16-checkpoint> --config configs/audit-e2-fresh.json
python3 tools/analyze.py --run results/e2 --paper-output artifacts/e2
python3 tools/analyze.py --campaign results/e1-marking-controls --paper-output artifacts/e1-controls
```

未知字段、theta 数组长度错误、训练/test 混用、错误 checkpoint/算子组合应报错；缺 SuiteSparse 不得静默切换生产后端。资源上限不得通过 `0` 的歧义值解除，使用显式可校验表达。

## 11. 数值及论文验收协议

### 11.1 必要测试

仅写验证数学/状态行为的测试，不写大量镜像实现细节的测试。

1. E1/E2 制造解源项、边界残差、积分收敛；基础 FEM 的离散残差和精确误差；加细嵌套注入保持 P1 函数。
2. 准插值投影及几何顶点/自由节点区别；局部 Riesz 的能量/核约束；全/空/部分 mask 的 AS 恒等式；相同 patch 合并后的正确重数；global eta 与区域 eta 的关系。
3. 归一化 family 均值质量、最坏样本补标条件、稳定并列排序；强残差边界/跳跃项与旧版一致。
4. AOT 方程、raw 块恒等式、base-test 与实际 PG 残差、缓存命中/失效：RHS 变、ell 变、同维数系数变、稀疏位置变、自由边界节点变均覆盖；POD 采用 correction 坐标而非只做绝对解快照截断。
5. 训练达标/停滞/上限、压缩失败回退、空区域；继承后新 kernel 和区域支撑；接受状态与升层前状态分离。
6. m_ref 和 lazy 的状态计数、强制升层/不触发/到上限、未完成周期、恢复/离线诊断隔离。
7. 线程数 1 与多线程数值一致到约定容差；线程 affinity 与实际参与线程数可观测；不能只看 OMP_NUM_THREADS。
8. 四组标记对照的成员选择、固定 ell、逐状态 48 RHS、原 AFEM exact-stop 与 97 状态 fixed-horizon 的分离；每样本先插值再分位数、不外推、参考维数使用真实端点，main target crossing 不插值。
9. 大于 131072 字符的旧标记字段完整导入；缺 RHS/重复行/半写 checkpoint/解析失败明确报错；不把 solver 正常完成等同于审计完整，不自动取消无关作业；便携包在独立目录可复现。

建议固定网格 oracle 比较容差起点为 `abs<=1e-12 + 1e-8*scale`，PG/AOT 等残差参照原验证标准（旧审计使用 `<1e-7`）；不同类型量采用有意义 scale，近零不得放大成假错误。这是新回归的起点，不更改原求解容差。旧 AOT 优化的 1e-14 量级差是实测证据，不要求所有平台逐比特相同。

自适应等价性不能只检查终点容差：相同配置/后端先逐状态检查 mesh/hash、ell、rank、worst-member、标记集；发现分歧停在首个状态，区分浮点并列排序和数学差异。确由浮点归约造成分支改变时，补做共同冻结状态比较，标记新轨迹，不能修改阈值强行对齐。

### 11.2 当前论文 oracle（用于发现错版本，不能作为硬编码答案）

| 量 | E1 | E2 inherit |
|---|---:|---:|
| 接受状态数 | 51 | 33 |
| 最终在线 DoF | 2422 | 8911 = 8910+1 |
| 最终参考 DoF | 316767 | 108169 |
| 最终 ell | 3；state 12 从 2 升 3 | 3；H_step 8 从 2 升 3 |
| 名义相对精确能量误差 | 0.009205652787925123 | 0.008294536369324071 |
| 相对 reference floor | 约 0.0039472696 | 0.004130545463867472 |
| 相对 reference gap | 约 0.0083161506 | 0.007192366436660725 |
| 论文尾部 DoF 阶 | 0.644557；最后 5 个完整周期 | 0.625629；最后 10 状态 |

当前 E2 冻结终态的 24 test 误差中位数/p90/max 为约 `0.008450/0.009000/0.009205`，max/nominal 约 `1.110`。这支持固定单一奇异形状加附近波包的低秩复用，不推广到多个独立奇异形状，也不以不同网格对比断言区域 AS 单独贡献。

历史可选 oracle：75 次 E2 单 RHS 审计已完成；24 test 的 `100*(E_shared/E_fresh-1)` 中位数在 H0/H16/H32 约 `-0.0217%、3.0310%、-1.4051%`。这些数据已从当前论文图表移除，只在显式启用旧审计复现时使用，不计入本版论文验收，不保证共享空间逐样本优于 fresh。

E1 tail 用 state 42/44/46/48/50，附 4/6 周期窗口；重复在线 DoF 的中间参考 sweep 保留曲线但不当作独立粗加细点混入拟合。E2 报 5/8/10 状态窗口和 R²。统一拟合 `log(error)=a-p log(N_online)`，p 是 DoF 阶，不直接称 h 阶。

E1 四组标记对照另外验收，不能替代上表主生产 oracle：

| 组别 | 完整状态 | 最终 N_on | 最终 N_ref | 最终 held-out max E | Figure 3 显示点数 / 末 state |
|---|---:|---:|---:|---:|---|
| nominal-ALOD | 51 | 2522 | 266247 | 0.01184514730 | 26 / 50 |
| family-ALOD | 51 | 2878 | 320194 | 0.01016350548 | 25 / 48 |
| nominal-AFEM | 97 | 192960 | 不适用 | 0.00606796517 | 81 / 80 |
| family-AFEM | 97 | 272164 | 不适用 | 0.00489073616 | 77 / 76 |

以上显示点数来自“每条保留轨迹首次测试 max<=0.012”的统一展示规则，不能作为运行停止条件。共同预算 500/1000/2000 的 family 相对 nominal 最大测试误差降幅分别为 ALOD `+7.51%/-2.85%/+2.50%`、AFEM `+1.03%/-0.24%/+2.94%`，负值表示变差。2000 时 ALOD 有 20/24 个测试成员改善、两组 reference-error p90 均约 0.620%；必须如实复现这种小且依赖预算的收益，不把“family 必须更优”写成测试断言。

同预算表的 E 均为百分数：2000 时 nominal/family-ALOD 测试 max 为 `1.352/1.318`，nominal/family-AFEM 为 `5.981/5.805`；两组 ALOD 的实际参考维数插值端点为 `149.1–198.3` 与 `132.9–179.3` 千。`E<=0.015` 的 family-ALOD/family-AFEM 首次保留状态分别用 1816/29700 在线 DoF，实际 max 为 0.01457/0.01497，ALOD 同时用 132888 参考 DoF。这是在线压缩证据，不是离线时间或峰值内存优势。

### 11.3 当前论文图表交付清单

本次核对的正文和附录共引用 7 幅图、8 张表。生成器按以下 LaTeX label/文件映射选择内容；编号只帮助核对，未来论文改动以新 label 清单为准。

| 位置 / label | 成果与数据口径 |
|---|---|
| Figure 1 `fig:e1-error-dof` / `e1_error_dof.pdf` | 变 ell E1 四方法名义 error–N_on；右图 ALOD exact/reference/gap；保留全部 51 状态，尾拟合用 5 个完整周期，不混入固定 ell 对照 |
| Figure 2 `fig:e1-final-meshes` / `e1_final_meshes.pdf` | 主生产 ALOD 粗/参考及原 87 点 AFEM 最终网格，区分节点数与自由 DoF |
| Figure 3 `fig:e1-family-marking` / `e1_family_marking_control.pdf` | 左 nominal/family-ALOD、右 nominal/family-AFEM；24 held-out 的中位数/p90/max；共同 1.2% 展示范围，原始完整轨迹保留 |
| Figure 4 `fig:e2-error-dof` / `e2_error_dof.pdf` | E2 四方法名义误差及 ALOD 三种误差分量；不默认画已删除的 fresh nominal 交叉标记 |
| Figure 5 `fig:e2-final-meshes` / `e2_final_meshes.pdf` | ALOD 粗/参考及完整 AFEM 最终网格；前两幅必须标注 `D=Omega intersect B_0.6(0)` 的阴影/虚线；R 不代表 LOD patch 物理半径 |
| 附录 Figures 6/7 `fig:e1-rhs-centres`、`fig:e2-rhs-centres` | E1 test/shift 中心及终态误差，E2 全部 48 样本中心投影及误差；E2 仍有振幅/相位等变化，不能把中心投影视为全部参数 |
| Table 1 `tab:adaptive-config` | 两主生产的有效配置；注明四组 E1 对照固定 ell=3 |
| Tables 2/3 `tab:e1-crossings`、`tab:e1-terminal-results` | 主生产 E1 离散首次目标越界及完整终态/尾拟合；普通 AFEM 为 87 点 |
| Table 4 `tab:e1-family-control` / `e1_family_control.tex` | 四组 N_on=2000 的名义/test 分位误差百分数、真实参考维数端点范围及 `E_ref` p90；不能漏普通 AFEM |
| Table 5 `tab:enrichment-config` | 区域 R、目标、块大小/工作 rank、停滞及两个压缩测试参数 |
| Tables 6/7 `tab:e2-terminal-results`、`tab:e2-rhs-methods` | E2 四方法完整终态/拟合；冻结最终空间 24 test 分位数及各自 max/nominal，非效率排名 |
| 附录 Table 8 `tab:e1-rhs-methods` | 原四种最终 E1 空间的绝对误差及自身 nominal 增长比；不作为受控 family marking 收益证据 |

共同统计规则：

- 主生产的等精度比较取全部已记录状态中首次 `E<=target`，不插值；E1 可以在中间参考 sweep 首次跨阈值。AFEM 在主图中是 nominal 自适应，ALOD 为 family 自适应，不把不同完整终点称为严格等精度计时比较。
- Figure 3 的 ALOD 只取初态和完成两轮参考细化的 state `0,2,...,50`；AFEM 取全状态。每条保留轨迹截到 max test error 首次不超过 0.012 的点，完整末段另存补充图，不删除原始记录。不能拿原 87 点 AFEM 的单次终态 family 审计拼造曲线。
- 四组共同预算固定为 500/1000/2000：每个同一 RHS 在相邻保留点作 log N–log E 插值，再对 24 test 作线性 type-7 分位数；不先插值分位曲线，不外推。参考相对误差同样逐 RHS 插值，N_ref 报实际两个端点范围。逐样本改善数由同 ID 配对计算。
- E2 主图截到各曲线首次不高于 0.008294536 的采样点；ALOD cutoff 来源为实际终态误差。完整基线终点表及 AFEM 完整末点网格另有口径，图注明示。新运行有数值变化时记录新 cutoff 和已画点表，不硬塞旧值。
- 当前不交付 `e1_rhs_growth.pdf` 作为 Figure 3，也不把未引用的 `e2_rhs_reuse.pdf`、`e2_rhs_samples.tex` 或 fresh nominal crosses 自动加入论文。这些旧成果仅放可选诊断输出目录。
- 可选图含 Theta、层间敏感性、区域内外 gap、训练预算/rank、同状态 rank-zero/AS、reference 深化及旧 balance。75 次 fresh 不列为重绘当前论文的前置条件。
- 原始 CSV 只读；派生数据含 source hash、实际状态/样本选择、拟合窗口、插值权重、参考端点和缺失值统计。NA 不填零，旁路状态不混入主曲线。生成表头用当前 `E_ref`，可保留原始 `F` 列映射。
- `tools/analyze.py` 产生数据、PDF/PNG、TeX table 和 artifact manifest，默认不改论文正文；需要装入稿件时另做受源文件哈希保护的复制/编译步骤，不运行旧的精确字符串替换脚本覆盖用户新稿。

## 12. 服务器运行与资源要求

实现阶段小测试在本地 WSL Ubuntu-22.04的`<WSL_ROOT>`；完整生产可通过powershell在 `ssh shuihan` 对应服务器隔离目录 `<REMOTE_ALOD_ROOT>`，用 tmux 管理。不得复用旧实验输出目录或已有 session 名。

部署脚本发现 CPU 拓扑与当前允许 CPU 集，优先分配互不重叠的物理核心；默认每作业 16 物理核心（不足则显式减少），BLAS=1。在创建 OpenMP worker 前设置进程完整 affinity，保留已验证 `OMP_PROC_BIND=false`/移除不匹配 `OMP_PLACES` 的起点；验证 worker 的实际集合。不要在一条已绑定单核的 shell 下仅增大线程数。

并发数由实测峰值和可用内存决定，不由 tmux 数量或空闲核数决定：

- 启动前读取 CPU 拓扑、MemAvailable、已有任务、对应规模峰值；给每作业估计峰值乘至少 1.5 安全系数。新作业预留加上活跃作业预计尚未增长的内存，必须仍留主机余量。
- 同一既有大内存服务器初始建议保留至少 80 GiB 或物理内存 20% 中的较大值；是配置起点，不是 OOM 保证。在小 WSL 用适合机器的更小作业和显式上限。
- 监测整个作业进程组 RSS、MemAvailable 和 swap，不仅主 PID。分别限制 reference/candidate 自由度、patch 工作量、稠密批大小及墙钟。
- E1 nominal-ALOD 全审计已测到约 81.52 GiB 峰值，不能用 E2 历史 31.14 GiB 替代它的预算。四组对照全审计与在线生产的内存模型分开；普通 AFEM 也只按其自身实测估算，不按终态 DoF 简单比例外推 LU 填充。
- 在昂贵 factor/mesh 分配前做估计和边界检查；进程组周期性 watchdog 作为补充。因子填充不可仅由 DoF 精确预测，留量不能保证所有输入绝不 OOM。
- 资源不足不启动下一批；支持在安全状态请求 checkpoint 后退出。失败不自动无休止重启或升高限制，保留明确停止原因。
- `completed` 必须同时满足目标状态数、数据完整和成功结束；资源/用户停止分别记 `paused_resource/paused_user`，不能当成数学收敛。
- 后处理失败保留 `validation_failed` 和已经提交的求解结果；只暂停依赖该验证的后续任务，不默认终止其他运行正常的独立实验。跨作业取消必须有明确资源或用户指令原因。

首次仅一组中规模性能/内存测量，验证后 E1/E2 可并行；fresh RHS 与误差审计可在不可变 checkpoint 上分批并行，但不要让每个 RHS 各复制一个大 AOT 分解而耗尽内存。优先同进程分解复用和有界多 RHS，是否进一步多进程以实测为准。

## 13. 发布前必须回答的核对问题

- [ ] 新项目是否在移开旧项目路径后仍能构建、smoke、导入数据和生成图？
- [ ] E1/E2 的实际算法是否由显式配置决定，而不是目录名/build hash？
- [ ] 三类基线是否可运行，已有论文基线能否验证导入？
- [ ] E1 四组对照是否共用显式成员标记接口，51/51/97/97 状态及 14208 行审计能否生成/验证导入？
- [ ] 是否区分 E1 变 ell 主生产、固定 ell 对照及 AFEM 的 87 点主基线/97 点对照？
- [ ] 全域 eta_H、区域 eta_D、Theta、相邻层差是否明确区分？
- [ ] Theta 是否保留正确粗能量分母，并明确 q 的载荷幅值依赖、Ritz 非严格上界及启发式与可解性条件的区别？
- [ ] reference sweep 是否真正用上一轮新解，并正确处理 m_ref=1/2/3？
- [ ] balance 已完全退出生产控制，reference/candidate FEM 解是否只在诊断出现？
- [ ] AOT 是否保持全域定义并正确复用分解，区域 kernel 是否完整迁移？
- [ ] E2 继承的是 kernel 代表，mesh-change 修复与 ell-only 重投影是否分开？
- [ ] 暂停恢复是否避免前缀重放，离线审计是否不污染主状态？
- [ ] 是否同时报告优化前/后数学一致性、实际 wall、峰值内存和诊断强度？
- [ ] E1/E2 图表是否由可追踪的数据生成，未计算项/截断/不同终点口径是否明确？
- [ ] 是否按当前 7 图/8 表交付，E1 使用四组绝对 held-out 误差，旧 E2 fresh 图/样本表不再默认生成到论文目录？
- [ ] 大字段标记、失败隔离、恢复主表来源、结果包独立解压复现是否验证？
- [ ] 是否记录仍未验证的事项：大规模 AOT 缓存总收益、32 线程收益、m_ref>2 科学效果、区域 AS 对 localization 的非严格分工？

最终交付应包含“已实现”“已复现”“仅有候选优化”三类清单，不能把计划中的性能改进写成已测得结论。
