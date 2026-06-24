# Model Card: Allen Consensus WMB Integrated Taxonomy (2025-10-31)

## Overview
- **Name**: Consensus Whole Mouse Brain Integrated Taxonomy
- **Version**: 2025-10-31 release
- **Institution**: Allen Institute for Brain Science
- **License**: CC-BY-4.0
- **URL**: https://alleninstitute.github.io/abc_atlas_access/descriptions/Consensus-WMB-taxonomy.html

## Data Description
Integrative consensus across scRNA-seq and snRNA-seq WMB taxonomies.
Provides unified cell type labels harmonized across sequencing modalities.

## Usage in mbanno
- Primary reference taxonomy for annotation
- Provides hierarchical labels (class → subclass → supertype → cluster)
- Used as ground truth in benchmark evaluation

## Citation
Yao, Z. et al. A high-resolution transcriptomic and spatial atlas of cell types
in the whole mouse brain. *Nature* 624, 317-332 (2023).
https://doi.org/10.1038/s41586-023-06812-z

## Download
```bash
mbanno download-reference --ref allen-consensus-wmb --version 20251031
```

## Limitations
- Based primarily on adult mouse brain (P56-P63)
- May not capture developmental or aged states
- Taxonomy is a consensus; individual samples may show deviations
