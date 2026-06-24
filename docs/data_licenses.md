# Data License Compliance

This document outlines the license terms for all reference datasets used by mbanno.

**IMPORTANT**: mbanno code is MIT-licensed. The reference datasets have their own
licenses. Users MUST comply with each dataset's terms.

---

## Allen Institute for Brain Science Datasets

### Allen WMB 10X v3 (Nature 2023)

- **DOI**: 10.1038/s41586-023-06812-z
- **License**: CC-BY-4.0
- **Terms**: Free to use, share, adapt for any purpose with attribution.
- **Attribution**: Yao, Z. et al. Nature 624, 317-332 (2023).

### Allen Consensus-WMB Integrated Taxonomy (2025-10-31)

- **Source**: alleninstitute.github.io/abc_atlas_access/
- **License**: CC-BY-4.0
- **Terms**: Same as above.

### MERFISH Whole Mouse Brain Atlas (Nature 2023)

- **DOI**: 10.1038/s41586-023-06808-9
- **License**: CC-BY-4.0
- **Terms**: Free to use with attribution.
- **Attribution**: Zhang, M. et al. Nature 624, 343-354 (2023).

### Slide-seq Cytoarchitecture (Nature 2023)

- **DOI**: 10.1038/s41586-023-06818-7
- **License**: CC-BY-4.0
- **Terms**: Free to use with attribution.
- **Attribution**: Liu, H. et al. Nature (2023).

---

## External Datasets

### STARmap PLUS CNS Atlas (Nature 2023)

- **DOI**: 10.1038/s41586-023-06569-5
- **License**: CC-BY-NC-4.0
- **Terms**: NON-COMMERCIAL USE ONLY. Cannot be used for commercial purposes.
- **Attribution**: Shi, H. et al. Nature (2023).

---

## Annotation Tools

| Tool | License | Commercial Use | Citation |
|------|---------|---------------|----------|
| MapMyCells | Allen TOS | Check terms | Allen Institute |
| scVI/scANVI | BSD-3 | Yes | Lopez et al. Nat Methods 2018 |
| CellTypist | MIT | Yes | Dominguez Conde et al. Science 2022 |
| SingleR | GPL-3 | Yes (with conditions) | Aran et al. Nat Immunol 2019 |
| scGPT | MIT | Yes | Cui et al. Nat Methods 2024 |
| ScType | GPL-3 | Yes (with conditions) | Ianevski et al. Nat Commun 2022 |

---

## mbanno Policy

1. **No raw data redistribution**: mbanno NEVER bundles raw matrices or paper
   content within the repository. Users download from official sources.

2. **License tracking**: Every downloaded dataset is accompanied by its license
   information in metadata files alongside the data.

3. **Marker genes**: Marker gene tables are script-derived from public data or
   curated from published literature (public domain when not copied from papers).

4. **Model weights**: Any pre-trained models will be published separately with
   clear training data provenance and license restrictions.

---

## Recommendations for Publications

When publishing results using mbanno, cite:

1. **mbanno**: See CITATION.cff
2. **Reference data**: Each dataset's original paper (see data_sources.yaml)
3. **Annotation tools**: Each method used (see table above)
