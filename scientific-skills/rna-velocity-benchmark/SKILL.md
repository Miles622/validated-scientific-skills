# RNA Velocity Benchmarking
**Category:** Bioinformatics & Genomics | **Tags:** #single-cell #rna-velocity #scvelo #scanpy

---

### 📝 Description
This skill provides a **validated pipeline** for calculating and comparing RNA velocity trajectories using **scVelo**. It is specifically optimized to avoid common mathematical stability errors by utilizing a curated, production-ready environment. 

The skill benchmarks the **Stochastic** and **Dynamical** models to evaluate their performance in inferring cellular trajectories. 

> **Note:** The dynamical model solves the full transcriptional kinetics for each gene—transcription ($\alpha$), splicing ($\beta$), and degradation ($\gamma$) rates—providing a more rigorous basis for inferring cellular latent time in complex lineages.

---

### ⚙️ Prerequisites
To maintain mathematical stability and avoid solver incompatibilities:

* **Environment**: Python 3.11+
* **Dependencies**: `scvelo>=0.3.0`, `scanpy`, `matplotlib`
* **System Requirements**: For WSL or headless Linux environments, the `matplotlib.use('Agg')` backend must be specified to prevent core dumps during image rendering.

---

### 🚀 Implementation: Benchmarking Pipeline
The following implementation utilizes the Pancreas dataset to compare steady-state stochastic assumptions against automated dynamical modeling.

```python
import matplotlib
matplotlib.use('Agg') 
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

def run_velocity_benchmark():
    # 1. Data Loading & Pre-processing
    adata = scv.datasets.pancreas()
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    
    # 2. Gene Selection & Stabilization
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    adata = adata[:, adata.var.highly_variable].copy()
    
    # 3. Compute First and Second Order Moments
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    # 4. Model Comparison
    # Stochastic Model (Probabilistic)
    adata_stoc = adata.copy()
    scv.tl.velocity(adata_stoc, mode='stochastic')
    scv.tl.velocity_graph(adata_stoc)

    # Dynamical Model (Differential Equations)
    adata_dyna = adata.copy()
    scv.tl.recover_dynamics(adata_dyna, n_jobs=1) 
    scv.tl.velocity(adata_dyna, mode='dynamical')
    scv.tl.velocity_graph(adata_dyna)

    # 5. Output Generation
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    scv.pl.velocity_embedding_stream(adata_stoc, basis='umap', title='Stochastic', ax=axes[0], show=False)
    scv.pl.velocity_embedding_stream(adata_dyna, basis='umap', title='Dynamical', ax=axes[1], show=False)
    
    plt.savefig('velocity_benchmark_results.png', bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    run_velocity_benchmark()

python3 benchmark_velocity.py
```

---

### 📊 Validation Metrics
The pipeline output is evaluated against biological ground truths using the following criteria:

* **Directional Correctness**: Vectors must align with the transition from *Ductal* progenitors to mature *Endocrine* clusters (*Alpha*, *Beta*, *Delta*).
* **Streamline Coherence**: Local neighborhood vectors should exhibit high cosine similarity, particularly in high-velocity "bottleneck" regions like `Ngn3-high EP`.
* **Kinetic Parameters**: The Dynamical model should identify gene-specific splicing rates ($\beta$) that produce recognizable **phase portraits**.
* **Latent Time**: Temporal ordering should correlate with established differentiation markers rather than simple distance-based pseudotime.

---

### 🛠 Installation & Usage
* **Directory Structure**: Store in `scientific-skills/rna-velocity/`.
* **Filename**: `SKILL.md`.
* **Execution**: Execute via `python3 benchmark_velocity.py`.

---

### 📖 References
* **Bergen et al. (2020)**, *Nature Biotechnology*: [Generalizing RNA velocity through dynamical modeling] (https://doi.org/10.1038/s41587-020-0591-3).
* **La Manno et al. (2018)**, *Nature*: [RNA velocity in single cells] (https://doi.org/10.1038/s41586-018-0414-6).

---