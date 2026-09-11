# ALOD

从 `LOD2d_C++` 提取的二维自适应局部正交分解研究代码，目标是独立复现论文 E1/E2、AFEM/UFEM/SLOD 基线及 E1 nominal/family 标记对照。

**当前为初始迁移检查点，不是论文复现发布版。** 已冻结来源与小规模回归数据，提取可独立构建的网格、边界、P1 Helmholtz 装配、积分和两种名义制造解。P0 完整验收及 P1 的基线驱动尚未完成；ALOD 调度、LOD corrector、区域 AS/AOT、恢复和论文图表生成仍在后续范围。详见 [迁移状态](docs/migration_status.md)。

## 构建与验证

支持 Linux / WSL Ubuntu 22.04，依赖 C++20、CMake 3.20+、Eigen 3.3+、OpenMP、SuiteSparse UMFPACK，以及 Python 3.10+。Python 工具只使用标准库。

```sh
sudo apt-get update
sudo apt-get install -y g++ cmake libeigen3-dev libsuitesparse-dev python3 git
cmake --preset release
cmake --build --preset release
ctest --preset release
```

默认使用 2 个构建任务；缺少 UMFPACK 会在配置时失败。新项目自行编译全部引入源码，不读取旧项目头文件或链接旧静态库。原生 Windows 构建尚未验证。

运行已实现的固定网格 FEM 检查：

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/alod_fem_smoke E1 6
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 build/alod_fem_smoke E2 6
python3 tools/check_repository.py
```

该入口固定使用 kappa=16、名义载荷与论文积分设置，输出 JSON。最后一个参数是 NVB 加细层级，限 0–8，用于小规模验证。它没有实现自适应生产运行或论文基线 horizon；粗网格误差较大不代表生产轨迹的终态精度。

## 当前内容

| 目录 | 用途 |
|---|---|
| `include/`、`src/` | 提取的网格/FEM 核和 E1/E2 名义问题入口 |
| `apps/fem_smoke.cpp` | 有限规模的固定网格数值检查 |
| `data/rhs/` | 两实验的 48 个正式样本参数；名义样本均为 0 |
| `tests/fixtures/` | 旧版 FEM oracle、E1 四组 canonical 审计、E2 accepted 审计 |
| `tools/` | 来源回归、数据完整性及本地路径检查 |
| `docs/provenance/` | 原始/引入哈希、来源清单、历史配置和版本核对 |
| `.github/workflows/` | Ubuntu 构建与 smoke CI |

仓库保留 14,208 条 E1 对照审计和 1,650 条 E2 接受状态审计作为明确列出的回归输入。E2 每状态的 50 条包含 48 个正式样本及 2 个 pure 诊断样本。它们是**导入的历史结果**，不能当作本项目新计算的生产结果。

旧 JSON 仅放在 `docs/provenance/configs/` 供审计，其中仍含历史 `build_hash`、`WORKTREE`、旧控制字段及已过期论文哈希；不作为新运行接口。有效语义与差异见 [实验说明](docs/experiments.md) 和 [论文实现映射](docs/paper_implementation_map.md)。

## 来源与后续迁移

[计划书](ALOD_SUBPROJECT_AGENT_PLAN_20260909.md) 保留既有文件名；执行记录以 [migration_status.md](docs/migration_status.md) 为准。已核对的当前论文含 7 幅图和 8 张表，旧 E2 fresh 图及五样本表不属于当前必需交付。

重建旧版小网格 oracle 时，显式提供本机归档源根目录变量；通常使用仓库内已冻结的 oracle 即可：

```sh
python3 tools/build_legacy_oracle.py --source "$E2_SOURCE_ROOT"
```

该开发命令在被忽略的 `_local/` 内编译旧核，仅在明确刷新回归输入时使用。生成器不会运行旧生产轨迹。不要用新输出直接覆盖 oracle 来消除数值差异。

所有提交文件必须使用仓库相对路径或符号根。构建目录、运行结果、日志、缓存、机器配置和凭据均已忽略；提交前检查暂存区实际内容：

```sh
git add .
python3 tools/check_repository.py --staged
git diff --cached --check
```

来源与许可情况见 [NOTICE.md](NOTICE.md)。
