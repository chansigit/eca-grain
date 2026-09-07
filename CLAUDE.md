# eca-grain: GRAIN（Guarded Rigor-Audited In-lineage Neighborhoods）

把已发布的 ECA-RSI 单元（`<unit>/release/final.h5ad`）一次性聚合成 grain（约 γ = 20 个表达相近的细胞，
counts 求和，**输出保留全部基因**），供 GRN 推断等下游使用。生态位置：eca-pp → eca-rsi（osp / msp / zmip）
→ release → **eca-grain**；只读 release，不改它。GitHub `chansigit/eca-grain`（公开），发行名 `eca-grain`，
导入名 `ecagrain`，命令 `eca-grain run|figures|validate-rigor`（= `python -m ecagrain`）。
名字里不用 metacell（用户要求，避免与他人重名）；代码注释和旧文档里的 metacell 一律指 grain。

## 运行（本机）

- 集群规则见组织级 CLAUDE.md；本会话通常已在 Slurm 计算节点 allocation 内（`$SLURM_JOB_ID` 已设），可直接跑 Python。
- **不要 pip install 进共享 venv** `/scratch/users/chensj16/venvs/dl2025/.venv`（正在跑批量作业的 ecarsi / osp / msp / zmip
  都 editable 装在里面）。本仓库只靠 PYTHONPATH：

  ```bash
  PY=/scratch/users/chensj16/venvs/dl2025/.venv/bin/python
  PYTHONPATH=/scratch/users/chensj16/projects/eca-grain $PY -m ecagrain run <final.h5ad> <outdir> [--sample-col project]
  ```

- 运行产物一律放仓库外。**正式输出**（用户 2026-09-07 拍板）：与 `rsi/` 并列的 `eca-pp/<Tissue>/grain/`，
  即 `$OAK` 下 `.../<dataset>/eca-pp/<Tissue>/grain/{report.html,metacells.h5ad,membership.parquet,summary.json,...}`。
  批量跑用 `$SCRATCH/eca-grain-jobs/`（`units.txt` 清单 + `grain_array.sbatch`，日志在 `logs/`）。
  开发用参考输出 `$SCRATCH/eca-metacell-dev/`：`bladder` / `fu2022` / `mca3_prostate` / `tmfacs_heart` 各有 `report.html`，
  合页 `results.html`；`*_v*`、`*_fine`、`*_coarse` 目录已过期。
- 输入是 eca-rsi 批量跑的 release（`$OAK` 下各数据集 `.../rsi/units/<unit>/release/final.h5ad`）。
  旧结构的 Fu2022 release 没有 `eca_sample_id`，用 `--sample-col project`。
- 给用户只报结果目录，不给投射 URL（用户要求）。节点无 CJK 字体，图内文字一律英文。
- venv 里 scikit-misc 与 numpy 不兼容，scanpy 的 seurat_v3 HVG 不可用，所以 `hvg.py` 自写 vst（与 Seurat 重合 1949/2000）。

## 已定设计（用户拍板；改动前先确认）

- one-through：参数固定（`run.py` 的 `DEFAULTS`），不逐数据集调参，不做人工循环。
- 格子 = sample × **coarse lineage**（`zmip_ann_coarse`）；fine 标签只记 `audit_majority` / `audit_purity`，不参与分组。
  坐标按 lineage 跨样本合并计算；lineage 不足 100 细胞回退到全局 `X_pca_harmony`。
- grain 数量与细胞数成比例，不做 1:1:1 平衡。
- outlier 按 MetaCells 2 gaps 规则丢，但每个细胞都在 `membership.parquet` 台账里；
  守恒 Σsize + outliers = n_cells，全局和逐格子都查，失败即报错。
- mcRigor 用 numpy 移植（与 R 包一致到 1e-9），不引入 R；dubious 只拆一层（地板 γ/2），拆不动的标
  `residual_dubious` 交付，不删。
- 屏蔽基因只影响 HVG 选择：细胞周期、解离应激（来自 osp）、热休克、线粒体、核糖体、Malat1；**血红蛋白不屏蔽**；
  输出 h5ad 的基因集与输入完全一致（测试锁定）。
- 结果页模板 `report.html` 已定稿：两张正方形 UMAP 共用右侧 lineage 图例（grain 叠在半透明原始细胞上；grain 自身
  UMAP 的 Leiden 1.0 聚类叠在 lineage 凹包岛上），下方按 grain 聚类 split 的 marker 热图；不显示 dubious / threshold /
  split 等中间结果；不加 panel 字母。改图先问用户。

## 开发约定

- ruff（行宽 120，规则见 pyproject）+ pytest；CI 在 GitHub Actions（lint + Python 3.10 / 3.12 测试 + wheel 导入）。
  提交前本地：

  ```bash
  ruff check ecagrain tests && ruff format --check ecagrain tests && PYTHONPATH=$PWD $PY -m pytest
  ```

- 版本号在 `ecagrain/__init__.py` 和 `pyproject.toml` 两处；改版本同时写 `CHANGELOG.md`（Keep a Changelog）。
- 文档分工：`README.md`（英文，当前真相）、`docs/design.md`（中文，按日期追加决策记录，不改写历史）、
  `docs/design-2026-09-06.html`（通俗版，已被前两者取代）。
- 依赖 `osp-sc` 只为 `DISSOCIATION_GENES_HS`；不引入 agent runtime。
- eca-rsi 仓库 `metacell` 分支的 worktree（`projects/eca-rsi-metacell`）里的副本已被本仓库取代。

## 未解决

- Smart-seq2（tabula-muris-facs）outlier 约 12%，UMI 数据 0.1% 到 2.4%。
- 每单元约两成细胞在 `residual_dubious` grain 里；mcRigor 的 dubious 更多跟 grain 大小相关，而非标签混杂。
- 首批 102 个单元（mca1.1 / mca2.0 / mca3.0 / tabula-muris-facs / tabula-muris-drop）的批量运行 2026-09-07 提交；其他数据集未跑。
- `validate-rigor` 依赖六月 `supercell2.0/outputs/<tissue>` 的 R 结果目录（用户旧工作目录，不在本仓库）。
