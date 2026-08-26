# 精简 ALOD 数值实验项目开发计划书

版本：v0.1  
日期：2026-08-26  
项目目录：`D:\code\femcode\ALOD\ALOD_minimal`  
论文基线：`../ALOD_paper_sectioned/`  
代码参考：`D:\code\femcode\LOD2d_C++`

## 1. 计划结论

新项目只服务于论文第 6 节的二维 Helmholtz 数值实验，不建设通用有限元框架，也不迁移
`LOD2d_C++` 的全部功能。实现采用“先正确、再缓存、最后并行”的顺序，并设置两种明确的
结果等级：

- `practical`：全部公式和算法路径已实现，但几何常数、局部混合求解或谱界仍使用普通浮点
  数值，只能报告论文记号中的 `U_ex^pr`；
- `certified`：`mu`、局部混合重构、最小奇异值下界、最大广义特征值上界以及最终标量传播
  均有向外舍入或可验证包络，才允许报告 `U_ex` 和“严格认证”。

首个可用版本应先复现论文的两组制造解和四类方法轨迹；最终版本再通过严格证书门禁。
参考工程当前已有大量可复用的 Helmholtz、网格和 LOD 基础代码，但其 reference-epoch/
candidate 流程与本论文“细网格不是参考解、误差直接控制到精确变分解”的语义不同，不能整套
搬入。

## 2. 项目目标与完成定义

### 2.1 目标

1. 实现二维复值、混合 Dirichlet/Neumann/Robin 边界的 Helmholtz P1 有限元。
2. 实现嵌套网格上的 Petrov--Galerkin LOD：投影型准插值、原/伴随局部校正子、超采样 patch
   和粗空间求解。
3. 实现每步使用的离散核局部 Riesz 指标 `eta_H`，以及检查点使用的完整精确解证书。
4. 实现论文 Algorithm 1 的 `H-h-ell` 自适应控制、cache batch 和失败原因分流。
5. 复现方形光滑局部振荡算例与 L 形低正则算例，并生成论文所需误差、证书、网格和成本数据。
6. 保证制造精确解只用于事后误差验证，绝不进入 MARK、CERTIFY 或 STOP 决策。

### 2.2 完成定义

项目只有同时满足以下条件才视为完成：

- 一条命令可从干净构建生成两个主算例的全部原始 CSV/JSON/VTU 和最终图表；
- 每个 ALOD 检查点记录 `eta_eq, beta_l_lower, rho_p, rho_a, eps_p, eps_a, beta_inf, U_ex`；
- `U_ex >= ||u-U||_kappa` 在所有制造解回归点成立；正式 `certified` 结果还须通过验证算术门禁；
- `H <= h` 的嵌套关系、`mu <= mu_0 < 1` 和所有状态机前置条件均 fail-closed；
- 缓存版本与 full-rebuild 版本在小网格上逐步等价；
- 比较方法使用同一问题、边界、误差口径和输出合同；
- 任何资源上限均结构化退出，不产生伪完整结果。

## 3. 科学与算法基线

### 3.1 模型与离散

实现论文中的

`-Delta u - kappa^2 u = f`

及齐次混合边界。能量范数为

`||v||_kappa^2 = ||grad v||^2 + kappa^2 ||v||^2`。

网格只保留两层数学角色：自适应粗网格 `T_H` 与计算校正子的细网格 `T_h`。二者始终嵌套；
`T_h` 不得被称作全局参考解网格。准插值固定为论文的
`I_H = E_H o Pi_H^dg`，并显式验证投影性和局部支撑。

### 3.2 LOD 主链

1. 检查分辨率数 `mu` 和核强制性常数 `c_W`。
2. 在 `N^ell(T)` 上分别求原问题和伴随问题的局部核校正子。
3. 组装计算试验空间 `(I-Q_l,h)V_H` 和测试空间 `(I-Q_l,h^*)V_H`。
4. 求解非 Hermitian 的粗 Petrov--Galerkin 系统。
5. 对每个粗顶点求 `b_kappa` 局部 Riesz 表示，形成 `eta_H,z`、`eta_H,T` 和粗 Dörfler 标记。

### 3.3 精确解证书链

检查点必须按以下顺序执行，任何一步失败都不得继续宣称证书：

1. 用 compatibility correction 构造完整残量的 patchwise `RT2/P2`（或等价
   `BDM3/P2`）平衡通量，得到 `eta_eq` 和数据振荡项；
2. 对原/伴随 total corrector defect 构造 quotient residual majorant；
3. 由广义最大特征值问题得到有保证的 `rho_p`、`rho_a` 上界，并计算
   `eps_p = rho_p/c_W`、`eps_a = rho_a/c_W`；
4. 从粗 Petrov--Galerkin 矩阵和能量 Gram 矩阵得到 `beta_l` 的有保证下界；
5. 计算
   `beta_inf = pos(beta_l_lower - rho_p*rho_a/c_W) / ((1+eps_p)(1+eps_a))`；
6. 仅当 `eps_a < 1`、`beta_inf > 0` 且总缺陷接受条件通过时，输出
   `U_ex = (1/c_W + eps_a/((1-eps_a)*beta_inf))*eta_eq`。

普通双精度 Ritz 值不能直接作为最大特征值上界，普通 SVD 值也不能直接作为最小奇异值
下界。开发期允许输出 `practical`，正式证书必须使用残量包络、区间算术或等价的验证算法。

### 3.4 自适应状态机

每个迭代严格执行：

`ADMISSIBILITY -> BUILD/REUSE -> SOLVE -> ETA_H/MARK -> CHECKPOINT? -> CERTIFY -> REPAIR or STOP/REFINE_H`

- `eta_H` 的相对下降或最大延迟触发昂贵检查点；改变 `T_h` 或 `ell` 后强制检查；
- 证书失败时先计算离散 localization diagnostic；超过阈值则增加 `ell`，否则根据原/伴随
  主导模态的 fine-defect 指标局部加密 `T_h`；
- `T_h` 或 `ell` 改变即开始新的 cache batch；粗网格局部加密只失效受影响 patch；
- 尝试终止前必须执行一次新检查点，禁止沿用旧的 lazy 数值。

## 4. 数值实验矩阵

### 4.1 E0：小网格代数校准

- 方形与 L 形网格的边界标签、NVB 闭包和嵌套检查；
- P1 Helmholtz 装配与直接 FEM 制造解收敛；
- `I_H` 投影性、核约束、原/伴随校正子恒等式；
- RT2/P2 patch compatibility、散度、法向通量和全局残量审计；
- 小矩阵上将矩阵自由谱界与显式 dense 结果比较；
- full-rebuild 与 cache 路径逐迭代比较。

E0 通过后冻结 `theta_H, theta_h, mu_0, tau_loc, tau_tot, q_cert, m_cert`，正式运行期间不得
按图形效果单独调参。

### 4.2 E1：方形局部光滑振荡解

问题完全采用论文第 6.1 节：`Omega=(0,1)^2`，上下 Dirichlet、左侧 Neumann、右侧
Robin，局部中心 `(3/4,1/2)`，`alpha_loc=80`，相位 `exp(i*kappa*x)`。

比较：

1. exact-certified adaptive LOD；
2. 固定 `h`、固定 `ell` 的 LOD；
3. uniform P1 FEM；
4. residual-adaptive P1 FEM。

输出误差/粗自由度、误差/总局部工作量、证书 effectivity、粗/细网格、cache 命中和检查点成本。

### 4.3 E2：L 形低正则解

问题完全采用论文第 6.2 节：凹角在原点，角点项为
`b_boundary*r^(2/3)*sin(2*theta/3)`，另加中心 `(-1/2,1/2)`、`alpha_osc=25`、振幅
`0.25` 的局部振荡包。禁止恢复参考工程中的旧 radial-cutoff 版本。

比较：

1. exact-certified adaptive LOD；
2. 固定超采样、均匀细网格的 LOD；
3. adaptive P1 FEM；
4. 成本允许时增加 uniform P1 FEM。

重点检查粗网格是否同时追踪凹角与振荡包，以及 fine-defect 网格是否保持其独立的证书角色。

### 4.4 E3：波数稳健性

v0.1 建议使用 `kappa in {8,16,32}`；最终数值在 E0 后冻结。固定制造解空间位置参数和统一
的 `kappa*H` 初始分辨率上限，报告证书前因子
`C_cert = 1/c_W + eps_a/((1-eps_a)*beta_inf)`。只论证可靠性常数的波数稳健性，不声称
总计算量与波数无关。

## 5. 最小工程结构

```text
ALOD_minimal/
  CMakeLists.txt
  README.md
  DEVELOPMENT_PLAN.md
  configs/                 # 冻结的 E0/E1/E2/E3 配置
  include/alod/
    mesh/                  # 三角网格、NVB、patch、嵌套映射
    fem/                   # P1 与 RT2/P2 局部混合空间
    helmholtz/             # 边界、积分、算子和制造解
    lod/                   # I_H、原/伴随校正子、PG 系统
    estimator/             # eta_H、localization diagnostic
    certificate/           # eta_eq、defect、谱界、U_ex
    adaptive/              # Algorithm 1 与 cache batch
    experiment/            # 配置、结果等级、运行状态
    io/                    # CSV/JSON/VTU
  src/                     # 与 include 对应的实现
  apps/alod_run.cpp        # 唯一正式运行入口
  tests/                   # 单元、恒等式、回归和小型端到端测试
  scripts/                 # 构建、运行矩阵、校验、绘图
  results/                 # 默认不提交大体积中间文件
```

依赖保持为：C++20、CMake、Eigen3；OpenMP 与 SuiteSparse 仅作为可选性能后端；严格验证构建
再启用 MPFR/MPFI/GMP。绘图脚本只依赖 Python、NumPy、pandas 和 Matplotlib。

## 6. 参考工程取舍

| 参考模块 | 决策 | 用途或原因 |
|---|---|---|
| `mesh/types`, `mesh/refine`, `mesh/edges` | 精简迁移 | 网格、NVB、父子关系和嵌套基础 |
| `helmholtz/boundary`, `quadrature`, `operators` | 精简迁移并复核 | 混合边界、复值 P1 装配、积分 |
| `lod/quasi_interp`, `lod/patches` | 迁移并按论文测试 | `E_H o Pi_H^dg` 与粗 patch |
| `helmholtz/patch_system`, `patch_solver`, `corrector`, `model` | 提取 DirectSaddle/DirectSchur 主路径 | 原/伴随局部校正子和 PG LOD |
| `helmholtz/benchmarks/paper_cases` | 只迁移论文当前的 R1/S 公式 | 制造解、梯度、forcing、边界验证 |
| `adaptive/estimator`, `kernel_residual` | 提取局部 Riesz 原语 | 每步廉价 `eta_H` 与诊断 |
| `adaptive/candidate_flux`, `certificates`, `verified_spectrum` | 逐公式审计后重组 | 可复用 RT/谱计算原语，但不能默认视为完整严格证书 |
| `reference_epoch_*`, `practical_driver`, `reference_retraction` | 不迁移 | 数学目标和状态机与本论文不同 |
| `singularity_hybrid` | 不迁移 | 论文当前 Algorithm 1 不需要 |
| `hp_*`, `schwarz_*`, `two_level_schwarz`, `shifted_laplacian` | 不迁移 | 不属于第 6 节最小实验闭包 |
| 椭圆 LOD、历史 benchmark/results/server runbook | 不迁移 | 降低构建、测试和语义负担 |

迁移策略不是复制整个静态库后删文件，而是每个工作包只引入已通过该工作包测试的最小源码
闭包。保留来源文件与参考 commit 记录，方便审计差异。

## 7. 工作包与里程碑

| 工作包 | 主要交付 | 退出条件 | 估算 |
|---|---|---|---:|
| WP0 工程与合同 | CMake、配置 schema、状态/结果等级、输出目录 | clean build；未知配置字段拒绝 | 2--3 人日 |
| WP1 网格/FEM/算例 | NVB、P1 装配、混合边界、R1/S 制造解 | G0--G2 通过；FEM 收敛正确 | 5--7 人日 |
| WP2 Petrov--Galerkin LOD | `I_H`、patch、原/伴随校正子、PG 解 | manufactured LOD 与全局校正子回归通过 | 7--10 人日 |
| WP3 廉价指标 | 局部核 Riesz、`eta_H`、Dörfler、localization diagnostic | 残量恒等式和局部效率回归通过 | 5--7 人日 |
| WP4 平衡通量 | RT2/P2、compatibility correction、`eta_eq` | patch 与全局残量审计通过 | 8--12 人日 |
| WP5 total defect/稳定性 | quotient residual、`rho_p/a`、`beta_l/beta_inf`、`U_ex^pr` | dense 小问题交叉验证；无假证书 | 8--12 人日 |
| WP6 自适应与缓存 | Algorithm 1、`H-h-ell` 修复、cache batch | full-rebuild 等价；结构化退出 | 7--10 人日 |
| WP7 比较方法与实验 | fixed LOD、UFEM、AFEM、E1/E2/E3 runner 与绘图 | 一条命令复现 practical 矩阵 | 6--8 人日 |
| WP8 严格数值验证 | 验证谱界、向外舍入、证书证据链 | 正式结果可标记 `certified` | 10--20 人日 |

单人开发的 practical 版本约 8--11 周；严格 certified 版本预计再增加 2--4 周。估算不含大型
服务器主实验排队时间。

## 8. 验收门禁

| Gate | 验收内容 | 未通过时禁止 |
|---|---|---|
| G0 | clean configure/build/test；依赖可发现 | 迁移更多代码 |
| G1 | NVB 闭包、边界标签继承、`T_H <= T_h` 全程成立 | LOD 求解 |
| G2 | R1/S 的值、梯度、forcing、混合边界和积分收敛 | 主实验 |
| G3 | `I_H^2=I_H`、核约束、原/伴随校正子残量 | `eta_H` |
| G4 | `eta_H` 与直接离散核双范数一致；标记质量守恒 | 自适应粗加密 |
| G5 | RT2/P2 compatibility、散度、边界通量、全局残量恒等式 | `eta_eq` |
| G6 | `rho_p/a` 和 `beta_l` 与显式 dense 小问题一致且方向正确 | `U_ex` |
| G7 | `U_ex^pr` 覆盖制造解真实误差；失败状态 fail-closed | practical 图表 |
| G8 | full-rebuild/cache 逐步等价；命中与失效原因可审计 | 启用缓存计时 |
| G9 | 所有常数和代数误差有验证包络、最终区间向外舍入 | `certified` 标签 |
| G10 | 干净环境从配置重建 CSV/JSON/VTU/图，run ID 可追溯 | 论文最终结果 |

## 9. 输出与复现合同

每次运行写入 `results/<run_id>/`：

- `config.json`：规范化后的完整参数；
- `metadata.json`：论文哈希、参考代码 commit、当前 commit、编译器、依赖、硬件和结果等级；
- `iterations.csv`：每步 `H/h/ell`、DoF、标记、状态、动作、指标、检查点量、缓存与计时；
- `checkpoints.csv`：证书各分量、验证状态和失败原因；
- `summary.csv`：目标误差首次命中点和终止状态；
- `mesh_manifest.csv` 与选定 `.vtu`：initial、checkpoint、fine-repair 前后和 final 网格；
- `run.log`：完整诊断日志。

未计算的 lazy 字段写 `NA`，不得沿用上一检查点。计时分为 method time、事后 exact-error
evaluation time 和 artifact time；精确解不得出现在算法对象接口中。

## 10. 主要风险与控制

1. **参考代码语义污染**：旧 reference-epoch/candidate 变量容易误入新算法。控制方式是重新定义
   两网格状态对象，不迁移旧 driver。
2. **伪证书**：浮点 Ritz/SVD 数值方向不满足证明要求。控制方式是结果等级和 G9 硬门禁。
3. **RT2/P2 工作量被低估**：先在 WP4 做独立 patch 审计与内存测量，再接入自适应循环。
4. **L 形奇异积分不稳定**：使用论文给出的解析 forcing split，角点附近采用可验证的加密积分
   回归，禁止数值二阶差分生成 forcing。
5. **缓存破坏正确性**：所有 key 必须包含 mesh version、patch、`kappa`、`ell`、边界和
   interpolation 指纹；先 full-rebuild 对照再开启。
6. **算例调参造成选择性偏差**：所有控制参数在 E0 冻结，正式配置、代码和论文哈希写入 run ID。
7. **矩阵自由谱计算不收敛**：设置显式迭代上限和残量包络；无法给出正确方向的界时结构化退出。

## 11. 首轮开发冲刺

第一轮只执行 WP0--WP2：

1. 建立最小 CMake 目标 `alod_core`、`alod_run` 和 `alod_tests`；
2. 迁移网格/NVB、边界、积分、P1 Helmholtz 和 VTU 输出的最小闭包；
3. 固化 R1/S 当前论文公式并完成制造解测试；
4. 迁移 `I_H`、patch、DirectSaddle/DirectSchur 校正子与双侧 PG LOD；
5. 以小波数、小网格生成第一个端到端 `practical` LOD 解；
6. 产出“迁移源码清单 + 删除模块清单 + 基准测试报告”，评审通过后再进入证书实现。

首轮不接入 reference-epoch、自适应 candidate、hp-FEM、Schwarz、服务器脚本或历史结果。
