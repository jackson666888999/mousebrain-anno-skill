# mbanno — Mouse Brain Consensus Annotator

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Allen/BICCN](https://img.shields.io/badge/reference-Allen%2FBICCN-green.svg)](https://portal.brain-map.org/)

> **A reproducible consensus annotation framework for mouse whole-brain single-cell and spatial transcriptomics using Allen/BICCN reference taxonomies.**

**mbanno** 不是一个声称"准确率超过所有工具"的黑箱，而是一个**基于权威参考 + 多模型交叉验证 + 可复现实验基准**的注释框架。工具的性能宣称只有在公开 benchmark 上跑出具体结果后才做出，且会明确注明"在哪个数据集/层级/指标上优于哪些工具"。

---

## 核心设计原则

### ✅ 科学合规
- 不声称"超过所有工具"——只在有 benchmark 数据支撑时才做出具体比较
- 所有参考数据来自公开发表、可验证的来源（Allen Institute / BICCN / Nature）
- 每个结果附带置信度、method agreement、参考版本号

### ✅ 多证据 consensus
- 多个工具各自给出 label，最终通过加权投票 + 置信度评估输出
- 权重**不写死**——通过 benchmark 自动学习不同脑区/平台的权重
- 置信度不足时输出 `unassigned` 或 `ambiguous`，不强制过细标签

### ✅ 数据合规
- **不打包任何原始矩阵或论文内容**——用户通过脚本从官方源下载
- 每个数据集记录 citation、DOI、license、version、download date
- 完整的数据使用说明见 `docs/data_licenses.md`

---

## 参考数据源

所有参考数据来自真实、可验证的公开来源：

| 数据集 | 论文 | 许可 | 用途 |
|--------|------|------|------|
| **Allen Consensus-WMB 2025-10-31** | Yao et al. *Nature* 2023 | CC-BY-4.0 | 主参考 taxonomy |
| **Allen WMB 10X v3 (Nature 2023)** | Yao et al. *Nature* 2023 | CC-BY-4.0 | scRNA-seq 参考 |
| **MERFISH WMB Atlas** | Zhang et al. *Nature* 2023 | CC-BY-4.0 | 空间验证 |
| **STARmap PLUS CNS** | Shi et al. *Nature* 2023 | CC-BY-NC-4.0 | 空间参考 |
| **Slide-seq cytoarchitecture** | Liu et al. *Nature* 2023 | CC-BY-4.0 | 脑区结构验证 |

完整数据源清单：`data_sources.yaml`

---

## 注释方法集成

| 模块 | 工具 | 用途 |
|------|------|------|
| 官方脑图谱映射 | **MapMyCells** | Allen/BICCN taxonomy label transfer（首选） |
| 深度生成模型 | **scVI/scANVI/scArches** | batch-aware integration + 半监督注释 |
| 快速监督分类 | **CellTypist** | 大规模初筛 |
| 参考相关性 | **SingleR** | 稳健 correlation-based 注释 |
| Marker 规则 | **ScType / 内置 WMB marker rules** | 可解释复核 |
| 基础模型嵌入 | **scGPT** | 辅助 embedding |
| 空间验证 | MERFISH / Slide-seq / STARmap | 脑区一致性检查 |

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
