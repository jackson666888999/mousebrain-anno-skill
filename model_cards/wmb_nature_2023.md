# Model Card: Allen/BICCN WMB 10X v3 (Nature 2023)

## Overview
- **Name**: Whole Mouse Brain 10X v3 taxonomy
- **Version**: 20231215 release
- **Institution**: Allen Institute for Brain Science / BICCN
- **License**: CC-BY-4.0
- **DOI**: 10.1038/s41586-023-06812-z

## Data Description
~7M scRNA-seq profiles from adult mouse brain (P56-P63), ~4M after QC.
Integrated with MERFISH spatial data for spatial validation.
Primary paper for the WMB comprehensive cell type taxonomy.

## Usage in mbanno
- Training reference for scVI/scANVI models
- Benchmark ground truth
- Marker gene derivation

## Citation
Yao, Z. et al. *Nature* 624, 317-332 (2023).

## Download
Via Allen Brain Cell Atlas API:
https://alleninstitute.github.io/abc_atlas_access/

## Limitations
- 10X v3 chemistry only (may differ from v2, Multiome)
- Adult brain only (P56-P63)
- Marker genes are guide, not definitive for all conditions
