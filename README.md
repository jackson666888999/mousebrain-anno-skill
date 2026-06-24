# 🧠 mbanno — 小鼠全脑单细胞注释.skill

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Allen/BICCN](https://img.shields.io/badge/reference-Allen%2FBICCN-green.svg)](https://portal.brain-map.org/)

> **不要靠一个工具赌你的细胞类型。用 Allen/BICCN 权威参考 + 7 种算法交叉验证，让每一条注释都有据可查、有分可依。**

---

## 🎯 你写论文时最怕什么？

```
审稿人："你怎么确定这个 cluster 真的是 L5 ET 而不是 L6 CT？"
你：     "因为 Scanpy 的 marker 图……看起来像？"
审稿人："有证据吗？"
你：     🤷
```

**mbanno 的回答：**
```
细胞 #48291:
  consensus: L5 ET (置信度 0.91, 高)
  MapMyCells → L5 ET    ✓ (官方图谱)
  scANVI     → L5 ET    ✓ (深度学习)
  CellTypist → L5 ET    ✓ (统计回归)
  SingleR    → L5 ET    ✓ (相关性)
  Marker     → L5 ET    ✓ (Fezf2+, Bcl11b+, Slc17a7+)
  脑区       → 皮层一致 ✓
  ────────────────────────
  5/5 工具一致, 结论可靠。
```

---

## 💡 为什么选 mbanno？

| 不用 mbanno 时 | 用了 mbanno 之后 |
|---------------|----------------|
| 一个工具定终身，错了也不知道 | 7 个工具投票，谁对谁错一目了然 |
| 被审稿人质疑时无证据 | 每个细胞附 7 维置信度矩阵，可追溯到原始工具 |
| 不同批次/平台切换流程 | 同一套权重自动适配 scRNA/snRNA/空间 |
| "我们准确率最高"——无法验证 | 公开 benchmark，跑完再说，不夸大 |
| 不知道参考数据能不能商用 | data_sources.yaml 记录每份数据的许可 |

---

## 🔬 基于真实数据，不做空中楼阁

所有参考数据来自 **Nature 2023** 顶刊论文，CC-BY-4.0 许可，可验证可下载：

| 数据集 | 规模 | 来源 | 用途 |
|--------|------|------|------|
| **Allen Consensus-WMB** | ~7M 细胞 | Yao et al. *Nature* 2023 | 🥇 主参考 taxonomy |
| **MERFISH WMB** | ~10M 细胞 | Zhang et al. *Nature* 2023 | 🗺️ 空间验证 |
| **STARmap PLUS** | ~1M 细胞 | Shi et al. *Nature* 2023 | 🧭 跨平台验证 |
| **Slide-seq** | 全脑结构 | Liu et al. *Nature* 2023 | 🏗️ 脑区结构参考 |

完整清单：`data_sources.yaml` | 许可说明：`docs/data_licenses.md`

---

## 🧰 集成 7 种最主流注释方法

```
MapMyCells (官方)  ───── 0.35 ─┐
scVI/scANVI (深度)  ───── 0.25 ─┤
CellTypist  (快速)  ───── 0.15 ─┼──→ Consensus Label + Confidence
SingleR    (稳健)  ───── 0.10 ─┤
ScType     (规则)  ───── 0.10 ─┤
scGPT      (基础)  ───── 0.05 ─┤
空间一致性  ───────────── 0.05 ─┘
```

> ⚠️ 权重不写死——通过 benchmark 自动学习，不同脑区、不同平台各有一套。

---

## 快速开始

### 安装

```bash
pip install git+https://github.com/jackson666888999/mousebrain-anno-skill.git
```

或开发模式：

```bash
git clone https://github.com/jackson666888999/mousebrain-anno-skill.git
cd mousebrain-anno-skill
pip install -e ".[dev,benchmark,report]"
```

### 下载参考数据

```bash
# 下载 Allen Consensus-WMB taxonomy 元数据
mbanno download-reference --ref allen-consensus-wmb --version 20251031

# 列出所有可用参考数据源
mbanno download-reference --list
```

### 注释您的数据

```bash
# 基础注释
mbanno annotate --input query.h5ad --species mouse --tissue brain \
    --reference allen-consensus-wmb

# 完整 ensemble 注释 + 脑区约束
mbanno annotate --input query.h5ad \
    --methods mapmycells scanvi celltypist singler marker \
    --region cortex \
    --output results/
```

### 输出文件

```
results/
├── query.annotated.h5ad       # 带注释的 AnnData
├── query.annotations.tsv      # 完整注释表（含多工具标签 + 共识 + 置信度）
├── query.confidence.tsv       # 置信度详情
├── query.marker_validation.tsv # Marker 验证结果
└── query.method_report.html   # HTML 报告
```

### Benchmark

```bash
# 在公开数据集上评测工具组合
mbanno benchmark --config benchmark.yaml --output benchmark_results/
```

---

## 注释层级

输出 Allen/BICCN 标准层级：

```
neighborhood → class → subclass → supertype → cluster
```

- **推荐发表层级：subclass**（~100 种亚类，与 CCF 解剖结构对应）
- **高置信数据可用：supertype**（~300 种超型）
- **探索性分析：cluster**（~500+ 簇，需额外验证）
- **置信度不足时输出 `unassigned` 或 `ambiguous`**

---

## 置信度体系

每个细胞输出以下指标：

| 字段 | 说明 |
|------|------|
| `confidence` | 综合置信度 [0, 1] |
| `agreement_rate` | 多工具一致率 |
| `margin` | 最高分与次高分之差 |
| `entropy` | 预测熵（不确定性） |
| `marker_support` | Marker 基因支持度 |
| `region_consistency` | 脑区一致性 |
| `reference_version` | 参考数据版本 |

置信度等级：

- `high_confidence`：主工具一致，marker 支持，脑区一致
- `medium_confidence`：多数工具一致，marker 或脑区证据弱
- `low_confidence`：工具分歧大
- `ambiguous`：同类内亚型无法区分
- `novel_or_state`：偏离参考但 marker 支持已有大类，可能是状态变化

---

## 项目结构

```
mousebrain-anno-skill/
├── README.md
├── LICENSE
├── CITATION.cff
├── pyproject.toml
├── data_sources.yaml              # 数据源清单（DOI、许可、下载脚本）
├── model_cards/                   # 各参考模型说明
├── src/mbanno/                    # 核心 Python 包
│   ├── cli.py                     # CLI 入口
│   ├── io.py                      # 数据读写
│   ├── qc.py                      # 质量控制和预处理
│   ├── references.py              # 参考数据管理
│   ├── mapmycells.py              # Allen MapMyCells 接口
│   ├── scanvi.py                  # scVI/scANVI 封装
│   ├── celltypist_runner.py       # CellTypist 封装
│   ├── singler_runner.py          # SingleR 封装
│   ├── marker_rules.py            # Marker 规则引擎
│   ├── consensus.py               # 多工具 consensus
│   ├── confidence.py              # 置信度评估
│   ├── benchmark.py               # Benchmark 框架
│   └── report.py                  # HTML 报告生成
├── workflows/                     # Snakemake 流程
├── markers/                       # Marker 基因表
├── docs/                          # 文档
│   ├── data_licenses.md
│   ├── methods.md
│   └── benchmark_protocol.md
├── notebooks/                     # 教程
├── tests/                         # 测试
└── benchmark_results/             # Benchmark 结果（不提交到 git）
```

---

## Benchmark 框架

在声称性能之前，必须先跑公开 benchmark。支持的指标：

- Accuracy / macro-F1 / weighted-F1 / balanced accuracy
- Adjusted Rand index / hierarchical F1
- Calibration error / abstention-aware accuracy
- Cross-batch robustness / cross-platform transfer
- Rare cell recall / brain-region consistency

支持的评测拆分：

- leave-one-batch-out
- leave-one-brain-region-out
- leave-one-donor-out
- scRNA-train → snRNA-test
- snRNA-train → MERFISH/Slide-seq validation

---

## 开源合规

- **代码**：MIT License
- **不打包原始矩阵**：用户通过脚本从官方源下载
- **记录每个数据集的 license + citation + version**
- **模型权重单独发布**，注明训练数据来源和可商用限制
- **Marker 表**由脚本从公开数据重新计算，不复制论文大表

详见 `docs/data_licenses.md`。

---

## 引用

如使用本工具，请引用：

1. 本工具：见 `CITATION.cff`
2. **Allen WMB Atlas**：Yao et al. *Nature* 624, 317-332 (2023). doi:10.1038/s41586-023-06812-z
3. 各方法原文（CellTypist、scVI、SingleR 等），详见 `data_sources.yaml`

---

## 路线图

- [x] Phase 1: MapMyCells + marker rules + confidence report
- [x] Phase 2: scANVI + SingleR consensus annotation
- [ ] Phase 3: CellTypist + scGPT embedding（提升跨平台鲁棒性）
- [ ] Phase 4: 公开 benchmark，自动输出与现有工具对比
- [ ] Phase 5: 基于 benchmark 结果发布具体性能声明

---

**Slogan**: *A reproducible consensus annotation framework — not an unverifiable "best" claim.*
