# RNA Velocity Benchmark
**Category:** Bioinformatics & Genomics | **Tags:** #single-cell #rna-velocity #scvelo #scanpy

### 📝 Description
This skill provides a **validated pipeline** for calculating and comparing RNA velocity trajectories using **scVelo**. It is specifically optimized to avoid common mathematical stability errors (like the scalar conversion error) by utilizing a curated, production-ready environment [6, 7].

### ⚙️ Prerequisites
To ensure stability and "stronger, more reliable performance," this skill requires the following environment [8]:
*   **Python**: 3.11+ (Required for distutils compatibility)
*   **Numpy**: `1.23.5` (Pinned to prevent internal scVelo solver incompatibilities)
*   **Packages**: `scvelo`, `scanpy`, `matplotlib`

### 🚀 Practical Example: Stochastic vs. Deterministic
This code performs a side-by-side benchmark using the built-in pancreas dataset.

```python
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

# 1. Load Data
adata = scv.datasets.pancreas()

# 2. Preprocessing (Standard Scanpy Workflow)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()

# 3. Computing Moments
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata)

# 4. Benchmark Different Velocity Modes
# Stochastic Mode
adata_stoc = adata.copy()
scv.tl.velocity(adata_stoc, mode='stochastic')
scv.tl.velocity_graph(adata_stoc)

# Deterministic Mode
adata_det = adata.copy()
scv.tl.velocity(adata_det, mode='deterministic')
scv.tl.velocity_graph(adata_det)

# 5. Visualization (Side-by-Side Comparison)
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
scv.pl.velocity_embedding_stream(adata_stoc, basis='umap', title='Stochastic Mode', ax=axes, show=False)
scv.pl.velocity_embedding_stream(adata_det, basis='umap', title='Deterministic Mode', ax=axes[5], show=False)
plt.tight_layout()
plt.show()
📖 References
scVelo Documentation: Optimized RNA Velocity Modeling
Scanpy: Single-cell analysis in Python
🛠 Installation
Create Directory: Create a folder named rna-velocity-benchmark in your scientific-skills/ directory.
Save File: Save this content as SKILL.md inside that folder.
Automatic Discovery: Restart your AI agent; it will automatically find and use this skill when relevant to your tasks
.

---