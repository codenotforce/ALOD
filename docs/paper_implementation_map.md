# 论文—实现对应表

以 `provenance/paper_inventory.json` 的论文哈希与标签为准。下列“原入口”相对相应来源根，未引入的功能不是当前可调用 API。

| 论文定义/标签 | 原入口 | 当前状态 / 后续字段与回归 |
|---|---|---|
| E1/E2 制造解与混合边界 | `src/helmholtz/benchmarks/paper_cases.cpp`，`make_shifted_r1_paper_case`、`make_parameterized_boundary_gaussian_s_paper_case` | 已提取；`alod::make_problem` 选择名义 E1/E2；`tests/fixtures/fem_E1.json`、`fem_E2.json` |
| P1 Helmholtz 与正定能量 | `src/helmholtz/operators.cpp`，`assemble_helmholtz_operators`、`solve_helmholtz_fem` | 已提取；UMFPACK 固定网格残差与旧版解对照 |
| `eq:local-kernel-space`、`eq:as-local-riesz` | `src/lod/quasi_interp.cpp`、生产快照 `adaptive/kernel_residual.*` | 未迁入；核约束、全几何粗顶点、patch 重数及质量守恒须单独验证 |
| `eq:certifier`、`eq:localized-coarse-smallness`、`eq:reference-interval` | `adaptive/certificates.cpp` reference defect/Ritz 路径 | 未迁入；保持粗能量分母，区分 Ritz 近似与严格界 |
| `eq:regional-as`、`eq:mode-deflation` | E2 生产 `AdditiveKernelRieszContext`、AOT 优化目录 helper | 未迁入；完整 patch 区域 mask、两遍消去及 AS 恒等式 |
| `eq:aot-matrix`、`eq:coupled-block-system` | `E2-AOT-factor-reuse-20260909` runner/helper | 未迁入；全参考域逆伴随、基础 test、块残差和缓存失效 |
| `eq:regional-training-target`、`eq:pod-dictionary`、`eq:distance-gate`、`eq:residual-gate` | 同上区域训练/POD 路径 | 未迁入；固定 base 目标、归一化 correction 坐标与逐成员重新耦合检查 |
| `alg:singular`、`alg:smooth` | `reference_epoch_paper_runner.cpp` 与 epoch driver | 未迁入；m_ref=1/2/3、接受状态/ell 事件独立计数 |
| `fig:e1-family-marking`、`tab:e1-family-control` | E1 对照生产源、canonical 四组 CSV | 已冻结 14208 行；当前只检查完整性/样本/残差，未实现同预算统计与新算法重跑 |

源快照不能混用：`provenance/source_differences.json` 记录工作区与生产头文件差异、AOT 优化入口哈希；`dependency_inventory.json` 提供保守依赖索引。构建日志、主机路径及旧二进制不进入版本控制。
