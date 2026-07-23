"""
Step 3: RAG-powered Q&A agent over scRNA-seq marker genes, using Groq (Llama).

Flow:
  1. Retrieve top marker genes for the cluster in question (from CSV).
  2. Retrieve gene annotations for those genes (local KB + live fallback).
  3. Feed the retrieved, grounded context to the LLM.
  4. LLM answers ONLY using retrieved context, citing sources, and explicitly
     declines to answer when the context doesn't cover the question.

Requires: GROQ_API_KEY environment variable set (free at console.groq.com).
"""

import os
import pandas as pd
from pathlib import Path
from groq import Groq
from gene_kb import retrieve_annotation

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
MARKERS_PATH = DATA_DIR / "marker_genes.csv"

client = Groq()  # reads GROQ_API_KEY from environment
MODEL = "llama-3.3-70b-versatile"


def get_top_markers(cluster: str, n: int = 5) -> pd.DataFrame:
    df = pd.read_csv(MARKERS_PATH)
    subset = df[df["cluster"].astype(str) == str(cluster)]
    return subset.sort_values("score", ascending=False).head(n)


def build_context(cluster: str, n: int = 5) -> str:
    top_markers = get_top_markers(cluster, n)
    context_lines = []
    for _, row in top_markers.iterrows():
        ann = retrieve_annotation(row["gene"], use_live_lookup=True)
        context_lines.append(
            f"- {ann['gene']} (log2FC={row['log2fc']:.2f}, adj p={row['pval_adj']:.2e}): "
            f"{ann['summary']} [source: {ann['source']}]"
        )
    return "\n".join(context_lines)


def ask(cluster: str, question: str) -> str:
    context = build_context(cluster)

    prompt = f"""You are a bioinformatics assistant helping interpret single-cell RNA-seq
clustering results. Answer the question using ONLY the retrieved marker gene
context below. Cite the gene(s) and source you used for each claim. If the
context doesn't contain enough information, say so explicitly instead of guessing.

Retrieved marker gene context for cluster {cluster}:
{context}

Question: {question}
"""

    response = client.chat.completions.create(
        model=MODEL,
        messages=[{"role": "user", "content": prompt}],
        max_tokens=500,
    )
    return response.choices[0].message.content


if __name__ == "__main__":
    cluster_id = "0"
    question = "What cell type does this cluster likely represent, and why?"
    print(f"Q: {question}\n")
    print(ask(cluster_id, question))
