# Documentation — Concepts & Reference Papers

## 1. Core concepts

**Single-cell RNA-seq (scRNA-seq)** — measures which genes are active in
each individual cell from a tissue sample.

**Quality control (QC)** — before analysis, remove likely-broken cells
(too few genes detected), likely doublets (two cells captured as one, too
many genes detected), and damaged/dying cells (high mitochondrial gene %).

**Clustering (Leiden algorithm)** — groups cells with similar overall gene
expression into clusters, without needing to specify the number of
clusters in advance. Modern replacement for the older Louvain algorithm.

**Marker genes** — for each cluster, the genes significantly more active
there than in other clusters combined (computed via a Wilcoxon rank-sum
test here). These are the evidence for what a cluster biologically is.

**RAG (Retrieval-Augmented Generation)** — retrieve real facts from a
trusted source first, then have the LLM answer using only that retrieved
context, reducing hallucination and enabling honest "I don't know" answers.

## 2. How the 3 scripts connect

```
compute_markers.py  ->  marker_genes.csv (top genes per cluster, from-scratch clustering)
        |
gene_kb.py           ->  retrieval: gene symbol -> meaning (local dict + live API fallback)
        |
rag_agent.py         ->  combines both, asks Groq/Llama, grounded + cited answer
```

## 3. Reference papers (read in this order)

1. **CompBioAgent** — Wang et al., 2025 (preprint, bioRxiv). LLM agent for
   natural-language scRNA-seq querying and visualization.
   https://www.biorxiv.org/content/10.1101/2025.03.17.643771v1

2. **CellAgent** — 2024/2026, accepted as an ICLR 2026 poster (peer-reviewed).
   Multi-agent framework automating full scRNA-seq analysis workflows.
   https://www.researchgate.net/publication/380693960

3. **CellAtria** — Nouri, Artzi & Savova, 2026 (peer-reviewed, npj
   Artificial Intelligence, Nature). Agentic ingestion/standardization of
   scRNA-seq data via a chatbot interface — includes literature-driven
   dataset acquisition (URL/PDF parsing -> GEO retrieval).
   https://www.nature.com/articles/s44387-025-00064-0

4. **ELISA** — (preprint, arXiv). Survey/comparison of ~6+ agentic
   bioinformatics systems.
   https://arxiv.org/pdf/2603.11872


