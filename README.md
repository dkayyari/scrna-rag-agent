# scRNA-seq Marker Gene RAG Assistant

A retrieval-augmented Q&A agent that interprets single-cell RNA-seq
clustering results in plain English — built from raw, unlabeled data
through to a working, validated RAG pipeline.

Inspired by recently published agentic bioinformatics systems (CompBioAgent,
CellAgent, CellAtria), scoped down to one clear capability: **explaining
what a cluster of cells likely is, grounded in retrieved marker gene
annotations, with honest refusal when evidence is insufficient.**

## Why this project (and why it's honest about prior work)
LLM agents for scRNA-seq analysis are an active research area — this is
**not** a novel technique. This is a from-scratch reimplementation focused
on demonstrating the full pipeline and the RAG pattern, not a claim of new
science. See `DOCUMENTATION.md` for the papers this draws on.

## Pipeline
1. **`compute_markers.py`** — loads raw PBMC3k data (no pre-existing
   labels), runs full QC (removes near-empty cells, likely doublets, and
   high-mitochondrial/damaged cells using standard published thresholds),
   clusters from scratch with the Leiden algorithm, and computes marker
   genes per cluster via a Wilcoxon rank-sum test.
2. **`gene_kb.py`** — retrieval layer: looks up gene function from a small
   local curated KB (~17 canonical PBMC markers), falling back to a live
   mygene.info API call for the other ~13,000+ genes.
3. **`rag_agent.py`** — the RAG loop: retrieves marker genes + annotations
   for a cluster, grounds an LLM prompt (via Groq/Llama 3.3 70B) in that
   context, and answers questions with citations back to source genes.

## Validated results (real output from testing)

**Cluster 0 — correctly identified as T cells** (unlabeled, discovered from
scratch): cited CD3D, IL32, and LTB with correct log2FC/p-values, explicitly
noted that other retrieved genes (LDHB, TPT1) didn't provide cell-type
evidence rather than fabricating a role for them.

**Cluster 7 (36 cells, hardest case — no local KB genes matched)** —
correctly identified as antigen-presenting cells based on live-retrieved
annotations for HLA-DPA1, HLA-DRA, and CD74, and appropriately hedged
between B lymphocytes, dendritic cells, and macrophages rather than
overclaiming a single answer.

**Grounding/honesty test** — asked "what is the mutation rate of this
cluster's cells?" (information not present in any retrieved context). The
agent explicitly declined to answer rather than guessing, stating the
context didn't cover mutation rates. This is the core behavior that
distinguishes a grounded RAG system from a standard chatbot.

## Setup

```bash
pip install -r requirements.txt
```

Get a free Groq API key at console.groq.com (no credit card required), then:
```bash
export GROQ_API_KEY=your_key_here
```

## Run order

```bash
python src/compute_markers.py   # QC + cluster from scratch, ~2-3 min
python src/rag_agent.py          # asks a sample question
```

## Known limitations (worth stating explicitly, not hiding)
- Local gene KB covers only ~17 well-known marker genes out of 13,000+ in
  the dataset; everything else depends on a live API call (slower, needs
  internet, and mygene.info doesn't always have a summary for every gene).
- QC thresholds (200 < genes < 2500, mito% < 5) are the standard published
  values for this specific PBMC dataset — they would need to be re-tuned
  for a different tissue type or sequencing platform.
- Ribosomal genes (RPS*, RPL*) frequently appear as top-scoring "markers"
  for every cluster since they're highly expressed everywhere; they're
  statistically real but biologically uninteresting, and the local KB
  doesn't have entries for them.

## Reference papers
See `DOCUMENTATION.md` for the full reading list and how each system this
project draws inspiration from (CompBioAgent, CellAgent, CellAtria) differs
from what's built here.
