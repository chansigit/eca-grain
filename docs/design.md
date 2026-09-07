# ecagrain 设计记录：一条命令跑到底的 Build → Outlier → Diagnose → Recheck

> 2026-09-07 定名 **ecagrain / GRAIN**（**G**uarded **R**igor-**A**udited **I**n-lineage **N**eighborhoods，用户从两个候选中选定；grain = 一个聚合单元，γ 是 graining level），包 `ecagrain/`，CLI `python -m ecagrain`。用户要求不叫 metacell 以免与他人重名；本文档正文里的 metacell 一律指 grain。用法见仓库 `README.md`。

**日期** 2026-09-06 **分支** `metacell`（worktree `eca-rsi-metacell`，基于 main 4dff5cc）
**状态** 已实现并在四个单元上验证（本文档按时间顺序记录决策，最新决定在后面的小节）。本分支另立门户，不修改 `ecarsi/` 已有代码，不影响正在运行的批量作业。

## 目标与边界

- 输入：ecarsi 已发布的单元，`<unit>/release/final.h5ad`。用到：`layers['counts']`、`obs.eca_sample_id`、
  `obs.zmip_ann_coarse`（粗标签，分块边界）、`obs.zmip_ann_fine`（细标签，只作事后核对列）。
  截至 2026-09-06 Oak 上有 86 个已发布单元。
- 输出：每个单元一套 metacell 表，写到 rsi 运行目录**之外**的独立根目录（避免触发 organize 的未声明 H5AD 守门）。
- 下游主要用途是 GRN 推断：metacell 个数要**大致正比于细胞数**，不做类型平衡。
- "metacell"泛指"把细胞划分成同质群再聚合"，不特指某个算法。

## 总原则

1. **One-through。** 所有参数事先固定，一个单元一条命令跑到底，中间没有人工检查点，
   不根据结果回头调参重跑。诊断信息只进报告，不进控制流。失败条件只有两种：输入缺失、守恒检查不过。
2. **一个细胞都不丢账。** 细胞可以不进 metacell（outlier），但 membership 表里每个细胞都有一行、有去向、有原因；
   所有 metacell 的 size 加上 outlier 数必须等于 `final.h5ad` 的细胞数。
   这是本方案的约定，不是 metacell 算法的共性：MC2 会把 outlier 放进未分配桶，mcRigor 的默认用法是剔除可疑组，
   SuperCell 和 SEACells 全分配。
3. **只拆不合、不重洗牌。** 复查只在原分组内部细化，已通过的分组一个不动。
4. **确定性。** 固定随机种子，同样输入产出同样结果。

## 核心模型

metacell = 一个划分 + 一个聚合。划分的"同质"由块内 embedding 上的近邻关系定义，标签只画外圈。
聚合存 **counts 之和 + size**，不存归一化平均。

## 流程

### ① Build：SuperCell 算法，Python 实现

1. **分块** `eca_sample_id × zmip_ann_coarse`（2026-09-07 用户暂定：按粗 lineage；9 月 6 日晚曾试细标签，
   Fu2022 对照显示两者的 metacell 聚类和 marker 一致（ARI 0.80），而细标签分块把多样本数据切碎、比例失真，故回到粗标签；
   `--label-col` 可改）。样本分块保样本组成、避开 Harmony 过校正；标签是不可跨越的边界。
   细标签作为 `audit` 列记录每个 metacell 的多数票和纯度。
2. **按粗标签重嵌入，跨样本合并算坐标**：坐标的来源和分组的范围分开。每个粗标签把所有样本的细胞合在一起
   选 HVG（`batch_key=eca_sample_id`，只留样本内也高变的基因）、算 PCA；分组仍只在样本 × 粗标签的块内做。
   块只借用坐标系，只和自己块内的细胞连边；批次偏移在块内是常数，不改变块内近邻关系。
   全局 PCA 看不见类型内的亚状态方向，这一步比标签细不细重要得多。这和 ZMIP 按 lineage 跨样本重嵌入是同一做法。
   合并后该粗标签仍 < 100 细胞时，不重算，直接用 `final.h5ad` 的全局 `X_pca_harmony` 行。
3. **walktrap 切分**：块内 kNN 图上 walktrap，切成 max(1, round(N/γ)) 个社区。γ 是目标平均值不是固定大小：
   70 个细胞切 4 组，大小由树决定（如 25/22/15/8）。比例只依赖组数 ∝ 细胞数，round(N/γ) 在块级别保证了这一点；
   大小不均只影响每组噪声，报告给大小分布。保留整棵层次树，复查时直接往下切。
   **上限**：walktrap 实测切出 1 到 50 个细胞不等（γ=20），大的那些几乎必被判 dubious；超过 1.5γ 的社区在自己的
   子图上继续切到不超过为止（确定性、只拆不合）。
4. **小块规则**：块内细胞数 < γ 时整块为一个 metacell，不需要坐标。PCA 维数取 min(30, 合并数/5)，k 取 min(15, 块内数 − 1)。
5. 每个细胞恰好一个 metacell。

用 scanpy + igraph 在 dl2025 venv 里写，零新依赖；不用 SuperCell R 包本身。整条流程纯 Python，不引入 R。

### ② Outlier：MC2 的细胞级 deviant 规则，固定阈值

规则改自 MC2 的 `find_deviant_cells`（gaps 策略）：
- **检验基因** = 该粗标签重嵌入用的那约 2000 个 HVG（已去 lateral 名单）。不检验全部基因：状态差异由这些基因定义，
  非 HVG 的偶发爆发不该把细胞踢出去。
- 对每个 metacell 的每个检验基因：算每个细胞的 log2(该基因 counts + 1) − log2(细胞总 counts)，把组内细胞按这个值排序，
  相邻值之间出现 ≥ 3（8 倍）的断层，断层外侧的细胞在这个基因上偏离；一个断层最多切出 3 个细胞且 ≤ 组大小的 10%。
- 某个基因在块内标出超过 3% 的细胞，视为 noisy 基因（爆发式表达），整块忽略它。这是 MC2 那个 3% 的真实含义：基因级保护，不是细胞级门槛。
- 每块上限 25%，超出时按偏离基因数取前 25%。一次通过，不像 MC2 那样剔除后再迭代。
- **实测后加的三道保护**（MC2 的规则是给 UMI 数据和百细胞级 metacell 定的，直接搬到 γ≈20 会过敏）：
  (1) 计数先按块内中位深度（上限 1 万）缩放再加 1，否则 Smart-seq2 每个 dropout 都是 8 倍断层；
  (2) 断层阈值取 max(3, 块内全部组内相邻断层的 99.9% 分位)，让噪声大的技术自动抬高门槛；
  (3) 至少在 3 个检验基因上偏离才算 outlier，单基因爆发是噪声，外来细胞在很多基因上偏离。
  UMI 数据上 outlier 率降到 0.2% 到 2.4%；Smart-seq2（tabula-muris-facs）仍有 12%，**未解决**，见实测一节。
outlier 不进任何 metacell，留在 membership 表里带原因。它保护两类细胞：残余 doublet，以及块里只有三五个、
否则会被硬塞进邻居组的稀有状态。
去掉 outlier 后 size 低于 γ/2 的 metacell 不解散（避免 MC2 式的连锁 dissolve），只在报告里标出。

MC2 需要人整理的三份名单在这里的处置：**屏蔽基因**（2026-09-07 定）= 细胞周期（Tirosh）+ 解离应激（van den Brink，osp 的名单）
+ 热休克（HSPA/B/D/E/H、HSP90、DNAJA/B）+ 线粒体 + 核糖体 + Malat1，只从 HVG 选择里剔除、不从输出矩阵里删；
**血红蛋白不屏蔽**（图谱里有红细胞数据集，Hba/Hbb 是它们的身份基因）；noisy 基因用表达量下限和 3% 占比规则替代。

### ③ Diagnose：mcRigor 组级检验，Python 重写

对 Build 产出（去掉 outlier 后）做 mcRigor `DETECT` 的检验，判每个 metacell trustworthy / dubious。
它检验组内基因相关结构，不看标签，粗标签下也能抓出亚状态混合。
不引入 R：算法读自本地包源码（`mcRigor_function.R` + `mc_test_stats.cpp`），核心不到百行，用 numpy 重写：

1. 全单元 log 归一化（CP10k + log1p）、seurat_v3 选 2000 个 HVG（= Seurat `NormalizeData` + `FindVariableFeatures` vst）。
2. 每个 size ≥ 2 的 metacell：取成员细胞 × HVG，去掉在 < 10% 成员里表达的基因，按基因 z-score 得 `dat`（细胞 × 基因）。
3. 统计量 T(dat) = ‖corr(dat) − I‖_F / √(p(p − 0.5))，即基因间相关矩阵偏离单位阵的程度。
   `TT_div = T(dat) / T(按基因列各自打乱的 dat)`：分母打乱破坏基因间相关、保留边际分布。
4. 零假设：每个细胞的值在基因间打乱（row permutation）得到无基因结构的 `dat0`，`null = T(dat0) / T(再列打乱的 dat0)`，重复 Nrep 次。
5. 阈值：按 size 取 null 的 1 − test_cutoff 分位数，跨 size 做 lowess 平滑（f = 1/6）；`TT_div > 阈值(size)` 判 dubious。

验证（已做，`python -m ecagrain validate-rigor`）：六月 `supercell2.0/outputs/Ear_g20/`（3055 细胞、152 个 metacell）
的 R 版结果对比：用 R 自己选的 HVG 时 T_org 逐个一致到 1e-9；用本包的 HVG 选择（与 Seurat vst 重合 1949/2000）时
T_org 相关 1.0000、TT_div 相关 0.99、阈值最大差 0.005、判定一致 91%（两边 Nrep=1 的置换噪声）。
细节：size < 5 的 metacell 不检验（统计量退化），状态记 `untested`；Nrep 默认 20（Nrep=1 时每个 size 只有几个 null 值，
分位数噪声让 dubious 率虚高 8 个百分点）。
阈值曲线（size → 阈值）在**整个单元**上估一次并保存，复查复用；按块估会因小块 metacell 太少而不稳。
dubious 不等于标错：分块下标签纯度恒为 1，"标签纯度"在粗标签下没有信息量。

### ④ Recheck：在 walktrap 树上往下切一层

1. 每个 dubious metacell 在自己的 walktrap 子树上切成两半（不重新聚类）。
2. 子块用保存的单元级阈值曲线按 size 插值分类，不重估 null。
3. 只拆一层，地板 γ/2。到地板仍不通过的标 `residual_dubious`，保留。
4. 膨胀上限：被拒部分最多变成两倍，整体比例漂移压在百分之十几以内。GRNBoost2 / pySCENIC 不接受样本权重，漂移必须在这里控制。

### ⑤ 交付

守恒检查通过才写输出：`Σ size + n_outlier == n_cells`，逐块同样成立。

## 比例保证

| 层面 | 保证方式 |
|---|---|
| 类型之间 | 固定 γ → 个数 ∝ 细胞数；Recheck 膨胀 ≤ 2× 且仅限被拒部分；outlier 按块上限 25%、按块报告 |
| 类型内部亚状态 | 近邻分组 → 亚状态比例从"细胞比例"变成"metacell 比例" |
| 精确还原 | 每个 metacell 带 `size`，按 size 加权；outlier 计数在表里 |
| 样本组成 | 按样本分块，metacell 不跨样本 |

metacell 不是生物学重复：跨条件统计的单位仍是样本。

## 固定参数（默认值，命令行可覆盖，但不做逐数据集调参）

| 参数 | 默认 | 说明 |
|---|---|---|
| γ | 20 | 每组目标细胞数 |
| 组大小上限 | 1.5γ = 30 | 超过的在子图上继续切 |
| 复查地板 | γ/2 = 10 | 只拆一层；size < 20 的 dubious 直接 residual |
| 最小检验大小 | 5 | 更小的记 `untested` |
| 重嵌入 HVG / PCA / kNN | 2000（batch_key=样本）/ min(30, n/5) / k=min(15, n−1) | 坐标按粗标签跨样本算；该标签 < 100 细胞时用全局 X_pca_harmony |
| 屏蔽基因 | Tirosh 细胞周期 + van den Brink 解离应激 + 热休克 + mt + ribo + Malat1（不含血红蛋白） | 只从 HVG 剔除，`genes.blocked_mask` |
| outlier | 检验基因 = 该类 2000 HVG（去 lateral）；深度缩放到 min(块中位, 1e4)；断层 ≥ max(3, 块内断层 99.9% 分位)；每组每基因最多 3 个细胞且 ≤ 10%；基因标出 > 3% 细胞则忽略；≥ 3 个基因偏离才算；每块上限 25% | 改自 MC2 gaps 策略，加三道保护 |
| mcRigor | Nrep=20, feature_use=2000, test_cutoff=0.05, gene_filter=0.1 | Nrep 从六月的 1 上调 |
| 阈值范围 | 单元级 | |
| 随机种子 | 0 | |

## 输出契约（每个单元）

```
<metacell_root>/<unit>/
  membership.parquet    cell, metacell_id(outlier 为空), eca_sample_id, coarse, fine, block,
                        status(member|outlier), outlier_reason, level(0|1), mcRigor
  metacells.h5ad        X = counts 求和;obs: size, eca_sample_id, coarse, fine_majority, fine_purity,
                        block, level, mcRigor, gamma;var 与 final.h5ad 一致
  threshold.tsv         单元级 mcRigor 阈值曲线(size, thre)
  summary.json          细胞数、块数、metacell 数、outlier 数与率(按块)、dubious 率、拆分数、residual 数、守恒检查
  report.md             人读汇总
```

## 首轮实测（2026-09-06，`ecagrain` 0.0.1，默认参数，输出在 `$SCRATCH/eca-metacell-dev/`）

| 单元 | 细胞 | 样本 | 块 | metacell | outlier | build 时 dubious | 拆分 | residual 细胞占比 | 用时 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| mca1.1 Bladder | 2193 | 1 | 6 | 143 | 2.4% | 22.5% | 5 | 23.3% | 19 s |
| mca3.0 Prostate | 19771 | 15 | 95 | 1313 | 0.2% | 19.6% | 35 | 20.6% | 96 s |
| Fu2022 半月板（`--sample-col project`） | 25357 | 11 | 61 | 1566 | 1.5% | 26.6% | 47 | 27.8% | 193 s |
| tabula-muris-facs Heart（Smart-seq2，30 板） | 3729 | 30 | 203 | 322 | **11.6%** | 14.5% | 5 | 18.8% | 40 s |

守恒检查四个单元全部通过。发现：
- **dubious 不跟细标签混合走，跟大小走。** Bladder 里 residual 的细标签纯度（0.87）反而高于 trustworthy（0.79）；
  TT_div 与 size 的 Spearman 0.30。大 metacell 更异质、检验也更有力；小 metacell 混了细标签也检不出来。
  所以 mcRigor 在这里抓的是"同类型内的连续状态跨度"，不是"混了两种标签"。细标签纯度列是独立的诊断信息。
- 每个数据集稳定有两成左右的细胞落在 residual dubious 里。多数 dubious 的 metacell 小于 20，按地板规则不拆；
  把地板降到 5 只会把它们送进更没有检验力的复检，不做。residual 带标记交付，由下游决定。
- Smart-seq2 的 outlier 率仍 12%：块内 99.9% 断层分位数只抬到 3.5 到 5.3，扩增噪声下仍有大量"多基因偏离"。
  **未解决**；全长数据在图谱里是少数（facs 20 个器官），先记录。
- Fu2022 旧结构 release 没有 `eca_sample_id`，`orig.ident` 对 endo/immune 两个来源缺失，用 `project`。

## 第二轮：样本 × 细标签分块（2026-09-06 晚）

| 单元 | 标签数 | 块数 | metacell | 大小 10% / 50% 分位 | outlier | 用时 |
|---|---:|---:|---:|---:|---:|---:|
| mca1.1 Bladder（1 样本） | 18 | 18 | 134 | 7 / 15 | 2.7% | 20 s |
| mca3.0 Prostate（15 样本） | 58 | 645 | 1472 | 2 / 13 | 0.1% | 77 s |
| Fu2022（11 样本） | 89 | 677 | 1817 | 2 / 14 | 0.8% | 198 s |
| tabula-muris-facs Heart（30 板） | 80 | 1386 | 1389 | 1 / 1 | 0% | 24 s |

细标签 × 样本把多样本数据切碎：Prostate 和 Fu2022 有一成 metacell 只有 1 到 2 个细胞；
Smart-seq2 Heart 的 1386 个块里大多数不到 20 个细胞，metacell 中位大小 1，等于没做。单样本的 Bladder 不受影响。
结果页（`figures` 子命令）只展示最终 metacell，且**可视化不用细标签**（89 个细标签画不下，用户拍板）：
metacell 矩阵重新算 HVG → PCA → Leiden 1.0 得到 metacell 聚类；UMAP 按粗标签和按 metacell 聚类各一张；
marker 热图的行 = 每个 metacell 聚类前 5 个 Wilcoxon marker 去重，列按聚类分块、块内层次聚类排序、块间按块平均谱层次聚类排序，
列色条 = 粗标签（每个 metacell 的多数票）。

**页面规范（2026-09-07 用户定稿）**，每个单元两张图，不带面板字母，UMAP 一律正方形：
- 上图一个面板两张 UMAP、共用右侧 lineage 图例。左 `Metacells (opaque, dot size ∝ cell#) on raw cells (translucent)`：
  单细胞半透明小点、metacell 不透明大点，都按粗 lineage 着色。右 `Metacell UMAP · Leiden clusters`：metacell 自身 UMAP，
  Leiden 编号，点不透明、不描边、在最上层；底层是各粗 lineage 的**岛**：每个 lineage 的 metacell 按 6 倍中位近邻距离
  连边取连通分量成岛，≥ 4 点的岛用 shapely `concave_hull(ratio=0.2)`，外扩 1.2% 轴宽画轮廓线、外扩 2.2% 画半透明填充；
  岛的名字不上图。
- 下图分块 marker 热图（见上），基因名 8.5 pt。
- 单细胞 UMAP、metacell 按 lineage 着色的 UMAP、大小分布、比例图、dot plot 都试过并按用户要求移除（数据仍在 summary.json）。
dubious / 阈值 / 拆分是中间过程，留在 summary.json 和 obs 列里，不上页面。
图内文字一律英文（节点无中文字体）。

## 唯一待拍板

`<metacell_root>` 放哪。建议独立根 `$OAK/data/sc/_metacell/<batch>/<organ>/`，
不放进 `<organ>/` 里，否则 rsi 批量作业 rsync 组织目录时要再加一条 exclude。

## 讨论记录来源

- 六月 `supercell2.0/`：SuperCell 2.0 γ=20 + mcRigor 在 Tabula Sapiens 27 组织上的结果。dubious 率多数 10% 到 30%，
  偏向上皮 / 分泌 / 稀有类型；把可疑细胞倒在一起按 γ 10/5/3 重分后保留率从七成到 97% 以上。
  本设计吸收"变粒度必需、阈值复用"两条，改为树上原地切一层，以控制比例漂移。
- mcRigor：Liu & Li, Nat Commun 2025 (s41467-025-63626-5)。机制以本地包源码为准。
- MC2 outlier 规则：Ben-Kiki et al. 2022 Genome Biology，`metacells` 包 deviants 默认参数。
