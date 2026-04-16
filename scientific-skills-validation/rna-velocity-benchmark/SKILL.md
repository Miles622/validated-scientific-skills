RNA Velocity Benchmark
Category: Bioinformatics & Genomics
Tags: single-cell, rna-velocity, scvelo, scanpy, benchmark
Description
This skill provides a validated pipeline for calculating and comparing RNA velocity trajectories using scVelo. It is specifically optimized to avoid common mathematical stability errors encountered in newer Python environments by utilizing a curated, production-ready environment. The skill includes a practical example that benchmarks the stochastic and deterministic velocity estimation modes, along with best practices for preprocessing and visualization. This serves as a foundational tool for researchers looking to analyze single-cell RNA-seq data with confidence in their velocity estimates.

--------------------------------------------------------------------------------
Prerequisites
This skill requires specific versioning:
Python: 3.11+ (Required for distutils compatibility).
Numpy: 1.23.5 (Pinned to prevent internal scVelo solver incompatibilities).
Key Packages: scvelo, scanpy, matplotlib.

--------------------------------------------------------------------------------
Practical Example
This code allows for a side-by-side benchmark of the stochastic and deterministic velocity estimation modes.
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

# 1. Load Data // Using the built-in pancreas dataset for benchmarking as an example
adata = scv.datasets.pancreas()

# 2. Preprocessing
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()

# 3. Computing Moments
# Calculate neighbors before moments to ensure stability in velocity estimation
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

# 5. Visualization
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
scv.pl.velocity_embedding_stream(adata_stoc, basis='umap', title='Stochastic Mode', ax=axes, show=False)
scv.pl.velocity_embedding_stream(adata_det, basis='umap', title='Deterministic Mode', ax=axes[11], show=False)
plt.tight_layout()
plt.show()

--------------------------------------------------------------------------------
References
scVelo Documentation: Optimized RNA velocity modeling
Scanpy

--------------------------------------------------------------------------------
How to Install this Skill
Create a folder named rna-velocity-benchmark in your scientific-skills/ directory.
Save the content above as SKILL.md inside that folder.
Restart your AI agent; it will automatically discover this skill and use it to help you with future single-cell tasks.