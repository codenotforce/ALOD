# 实验与配置边界

当前可执行程序仅为 `alod_fem_smoke`，接受问题 ID `E1` 或 `E2` 和可选 NVB 层级。下表冻结后续迁移目标，**不是已经实现的生产配置 API**。旧配置归档在 `provenance/configs/`，不能直接传给 smoke 程序。

| 轨迹 | E1 主生产 | E1 标记对照 | E2 主生产 |
|---|---|---|---|
| 问题 | R1，混合 D/N/R 方形域 | 与主生产相同 | S，凹角与 Gaussian 波包 |
| kappa、初始 H/h | 16，6/10 | 16，ALOD 6/10，AFEM 6 | 16，6/10 |
| 正式样本 | 16 train / 24 test / 8 shift | 四组同一组 48 个参数 | 16/24/8；pure 另计 |
| ell | 初始 2，最大 4，lazy | ALOD 固定 3 | 初始 2，最大 4，每接受状态检查 |
| 参考 sweep | 2，theta=[0.3,0.2] | ALOD 同左 | 1，theta=[0.2] |
| 标记 | 16 个训练成员，theta_H=0.15 | nominal=[0] 或 family=[0..15] | 耦合解训练族，theta_H=0.15 |
| 完整状态 | 51 | ALOD 51/51、AFEM 97/97 | 33 |
| 增广 | rank=0 | ALOD rank=0 | 区域 R=0.6、预算、POD、继承、全域 AOT |

积分均为 triangle/Gaussian/singular = 12/16/24，E1/E2 最大递归分别 8/6。当前 smoke 也使用这些值。E2 名义参数为奇异系数 1、波振幅 0.5、相位 0、Gaussian alpha=80、中心 (-0.5,0.5)。E1 中心为 (0.75,0.5)。

## 旧 JSON 之外的分支

- 四组 E1 通过 `E1_FAMILY_CONTROL`、`E1_FAMILY_MARKING` 和 `E1_CONTROL_OUTPUT` 进入/选择对照。ALOD 和 AFEM 复用对应算法驱动；恢复补丁只用于原归档重放。新接口需显式成员 ID，不能继承这些环境变量作为数学开关。
- E1 主生产额外 Theta 检查状态为 35/44/50，不能从同为 51 状态的固定 ell 对照推断。
- E2 历史 `build_hash` 含区域、预算、继承和 ell 控制标记，`diagnostic_sc_lod_campaign` 也参与选路。当前保存原值用于定位；未来必须把真正有效行为变为有类型的配置。
- 旧 `git_commit=WORKTREE` 和旧稿 `manuscript_sha256` 不是当前源码/论文版本依据，应使用来源清单中的实际文件哈希。
- E1 原 AFEM 主基线有 87 个状态；固定 horizon 对照有 97 个状态。当前图的展示截点与求解终止条件不同，不能由图重建 horizon。
- E2 `control_samples.csv` 可能包含升层前解。本次仅导入 `accepted_control_samples.csv`，并保留其中的纯分量诊断身份。

## 数据语义

E1 原字段 `E=e/exact_norm`、`F=f/exact_norm`、`G=g/exact_norm` 保持原样。论文称 `F` 为 `E_ref`；AFEM 的 reference 数据不适用，不把 NA 改为零。E1 样本表的 `shift` 对应原审计的 `shift-audit`。

E2 正式参数表仅含样本 0–47；接受状态另含 48/49 的 pure-corner/pure-wave。它们不参与正式 48 样本的计数，不把它们混入训练或测试族。

引入文件和变换均在 `source_manifest.json` 逐项登记。原始实验目录保持只读。长标记字段的 64 MiB 读取测试仅验证解析能力；原始 marking CSV 尚未迁入，本次不声称验证了所有 mean/worst bulk 或轨迹标记集合。
