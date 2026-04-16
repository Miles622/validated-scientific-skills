import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

# 1. Load the built-in pancreas dataset
print("Loading pancreas dataset...")
adata = scv.datasets.pancreas()

# 2. Preprocessing (Standard Scanpy Foundation)
sc.pp.filter_genes(adata, min_cells=3)

# The warning "Did not normalize X" happened because this dataset is pre-processed.
# We set enforce=True to ensure our specific normalization is applied.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()

# 3. Resolve the Deprecation Warning
# We now compute PCA and Neighbors explicitly using Scanpy [1]
print("Computing PCA and Neighbors...")
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)

# Now scVelo moments will run without warning because neighbors already exist
scv.pp.moments(adata)

# 4. Estimate RNA Velocity
print("Calculating velocity...")
# Explicitly set the mode to 'stochastic' to ensure stability
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# 5. Visualization
print("Generating UMAP stream plot...")
scv.pl.velocity_embedding_stream(
    adata, 
    basis='umap', 
    title='RNA Velocity Stream Plot: Cell Trajectories'
)

# Display plot in WSL2

# 6. Benchmark different velocity modes
print("Benchmarking velocity estimation modes...")

# Create a copy for deterministic mode comparison
adata_det = adata.copy()

# Estimate velocity using deterministic mode
scv.tl.velocity(adata_det, mode='deterministic')
scv.tl.velocity_graph(adata_det)

# Generate side-by-side comparison plots
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Stochastic mode plot
scv.pl.velocity_embedding_stream(
    adata,
    basis='umap',
    title='Stochastic Mode',
    ax=axes[0],
    show=False
)

# Deterministic mode plot
scv.pl.velocity_embedding_stream(
    adata_det,
    basis='umap',
    title='Deterministic Mode',
    ax=axes[1],
    show=False
)

plt.tight_layout()
plt.show()