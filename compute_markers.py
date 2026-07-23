"""
Step 1: Load raw PBMC3k data, perform quality control, cluster from scratch
(Leiden algorithm), and compute marker genes per cluster.

Unlike the pre-processed/pre-labeled version of this dataset, this script
starts from raw counts with no cluster labels, so cluster identity is a
genuine open question for the RAG agent to answer later.

QC thresholds used here match the standard, published Scanpy PBMC3k
tutorial thresholds (not arbitrary).
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)
MARKERS_OUT = DATA_DIR / "marker_genes.csv"


def load_qc_and_cluster():
    print("Loading raw PBMC3k dataset...")
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()

    print("Calculating QC metrics...")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    print(f"Before QC filtering: {adata.n_obs} cells x {adata.n_vars} genes")

    # Standard PBMC QC thresholds:
    # - remove near-empty cells (likely broken droplets)
    # - remove likely doublets (unusually high gene counts)
    # - remove damaged/dying cells (high mitochondrial %)
    adata = adata[adata.obs.n_genes_by_counts > 200, :]
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.filter_genes(adata, min_cells=3)

    print(f"After QC filtering: {adata.n_obs} cells x {adata.n_vars} genes")

    # Normalize + log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save a copy BEFORE scaling — needed for correct log-fold-change math later
    # (scaled/centered values can be negative, which breaks logFC calculations)
    adata.raw = adata

    # Reduce to most informative genes, scale, and run PCA
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)

    # Build neighbor graph and cluster with Leiden (modern replacement for Louvain)
    print("Clustering with Leiden...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata)

    print("Cluster sizes:")
    print(adata.obs["leiden"].value_counts())

    # Compute marker genes per cluster (uses .raw for correct, unscaled logFC values)
    print("Computing marker genes (Wilcoxon rank-sum test)...")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)

    return adata


def save_markers(adata, top_n=20):
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    rows = []
    for group in groups:
        for i in range(top_n):
            rows.append({
                "cluster": group,
                "gene": result["names"][group][i],
                "score": result["scores"][group][i],
                "pval_adj": result["pvals_adj"][group][i],
                "log2fc": result["logfoldchanges"][group][i],
            })

    markers_df = pd.DataFrame(rows)
    markers_df.to_csv(MARKERS_OUT, index=False)
    print(f"Saved {len(markers_df)} rows -> {MARKERS_OUT}")
    return markers_df


if __name__ == "__main__":
    adata = load_qc_and_cluster()
    save_markers(adata)
