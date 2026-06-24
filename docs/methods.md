# Methods — Consensus Annotation Framework

## Overview

mbanno implements a **layered consensus annotation** framework for mouse whole-brain
single-cell and spatial transcriptomics data, using Allen/BICCN reference taxonomies
as authoritative ground truth.

---

## Annotation Hierarchy

The tool outputs Allen/BICCN standard taxonomy levels:

```
neighborhood → class → subclass → supertype → cluster
```

| Level | Name | Count (~) | Usage |
|-------|------|-----------|-------|
| Level 0 | Class | 5-6 | Initial QC (neuron/glia/vascular/immune) |
| Level 1 | Neighborhood | ~15 | Functional grouping |
| Level 2 | Subclass | ~100 | **Recommended publication level** |
| Level 3 | Supertype | ~300 | Fine subtypes (high-quality data) |
| Level 4 | Cluster | ~500+ | Exploratory only |

Default publication level is `subclass`, which corresponds to CCF anatomical
brain regions and has the strongest cross-modality validation.

---

## Multi-Tool Ensemble

### Methods Used

| Priority | Method | Approach | Strength |
|----------|--------|----------|----------|
| 1 | MapMyCells | Allen official taxonomy mapping | Authoritative reference |
| 2 | scVI/scANVI | Deep generative model | Batch-aware integration |
| 3 | CellTypist | Logistic regression | Speed + interpretability |
| 4 | SingleR | Correlation-based | Robust to platform differences |
| 5 | Marker rules | ScType-like rule engine | Interpretable, biologically grounded |

### Consensus Formula

For each cell $i$ and annotation method $m$:

```
score_i = sum_{m} w_m * s_{i,m}  (for methods with agreement)
final_label_i = argmax(score_i)
agreement_rate_i = fraction of methods agreeing on final_label_i
```

Weights $w_m$ are:
- Initialized from benchmarks (default: MapMyCells 0.35, scANVI 0.25, etc.)
- Can be learned per brain region and per data modality
- NOT hardcoded — configurable and benchmark-driven

---

## Confidence Evaluation

Five confidence sources combined:

1. **Agreement rate** (0-1): Fraction of methods agreeing on the assigned label
2. **Prediction margin**: Score difference between top and second-best labels
3. **Marker support**: Expression level of canonical markers for the assigned type
4. **Region consistency**: Compatibility with known brain region-cell type mapping
5. **Entropy**: Prediction uncertainty across methods

### Confidence Levels

| Level | Criteria |
|-------|----------|
| `high_confidence` | Agreement > 0.8, margin > 0.3, marker support strong |
| `medium_confidence` | Agreement 0.5-0.8, reasonable margin |
| `low_confidence` | Agreement < 0.5, weak evidence |
| `ambiguous` | Multiple equally-scoring types at same level |
| `novel_or_state` | Deviates from reference but markers support basic class |
| `unassigned` | Insufficient evidence for any label |

---

## Region Constraints

Known brain region → cell type mappings from Allen CCF:

```
cortex → {L2/3 IT, L4 IT, L5 IT, L5 ET, L6 CT, L6b, Pvalb, Sst, Vip, ...}
hippocampus → {CA1, CA2, CA3, DG, Astro, Oligo, ...}
striatum → {D1 MSN, D2 MSN, Cholinergic, Astro, Oligo, ...}
...
```

Types not expected in a given region are flagged with lower confidence.

---

## Benchmark Protocol

### Metrics

- Accuracy (overall + per-class)
- Macro-F1 / Weighted-F1
- Balanced accuracy
- Adjusted Rand Index (ARI)
- Hierarchical F1 (penalizes distant misclassifications)
- Calibration error
- Abstention-aware accuracy
- Cross-batch robustness
- Rare cell recall

### Evaluation Splits

- Leave-one-batch-out
- Leave-one-brain-region-out
- Leave-one-donor-out
- scRNA → snRNA transfer
- snRNA → MERFISH/Slide-seq spatial validation

### Claims Policy

Performance claims are ONLY made when supported by reproducible benchmark results,
with explicit specification of:
- Which dataset
- Which taxonomy level
- Which metric
- Which baseline tools compared against
