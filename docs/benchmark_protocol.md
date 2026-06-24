# Benchmark Protocol

## Reproducible Evaluation of Mouse Brain Single-Cell Annotation

---

## Data Preparation

### Training Data

Allen/BICCN WMB 10X v3 reference:
- Source: s3://allen-brain-cell-atlas/releases/20231215/
- Split: 80% train / 20% test (stratified by subclass)
- Preprocessing: QC (min_genes=200, max_mt_pct=20), log-normalize, HVG 3000

### Validation Data

| Dataset | Modality | Cells | Use Case |
|---------|----------|-------|----------|
| WMB 10X test split | scRNA-seq | ~800K | Primary benchmark |
| MERFISH WMB | spatial | ~10M | Spatial consistency |
| Slide-seq | spatial | ~1M | Cross-platform |
| Hold-out donors | scRNA-seq | varies | Generalization |

---

## Task Definitions

### Task 1: subclass Classification
- **Input**: Normalized gene expression
- **Output**: Subclass label (~100 classes)
- **Primary metric**: Weighted-F1
- **Your data**: Allen WMB 10X v3

### Task 2: Cross-platform Transfer
- **Train**: scRNA-seq WMB 10X
- **Test**: snRNA-seq WMB Multiome
- **Metric**: Macro-F1, ARI

### Task 3: Spatial Validation
- **Train**: scRNA-seq WMB 10X
- **Test**: MERFISH cell-level predictions
- **Metric**: Region consistency score

### Task 4: Rare Cell Detection
- **Focus**: Cell types with < 0.1% abundance
- **Metric**: Recall@k

---

## Baselines

All baselines run with default parameters:

```yaml
baselines:
  - name: MapMyCells (official)
    method: allen_taxonomy_mapping
  - name: scANVI (vanilla)
    method: scanvi_label_transfer
  - name: CellTypist
    method: lr_classification
  - name: SingleR
    method: correlation
  - name: ScType
    method: marker_rules
```

---

## Running Benchmarks

```bash
# Prepare config
cat > benchmark.yaml << 'EOF'
tasks:
  - name: "subclass_classification_WMB_10X"
    dataset: "allen_wmb_2023"
    split: "leave_one_batch_out"
    level: "subclass"
  - name: "cross_platform_snRNA"
    dataset: "allen_wmb_multiome"
    split: "cross_modality"
    level: "subclass"
EOF

# Run
mbanno benchmark -c benchmark.yaml -o benchmark_results/
```

---

## Reporting Template

When reporting results, use this template:

```
On the Allen WMB 10X subclass classification task (held-out test split),
mbanno ensemble (MapMyCells + scANVI + CellTypist + marker rules) achieved:

- Weighted-F1: X.XX
- Macro-F1: X.XX
- Balanced accuracy: X.XX

Compared to individual baselines:
- MapMyCells alone: X.XX F1
- scANVI alone: X.XX F1
- CellTypist alone: X.XX F1
```

Do NOT write "exceeds all existing tools" without specifying exact
dataset/level/metric/numbers.

---

## Version Tracking

Every benchmark run records:
- Software versions (mbanno, scanpy, scvi-tools, celltypist, etc.)
- Reference data version (e.g., Consensus-WMB 20251031)
- Random seeds
- Hardware/environment
- Git commit hash
