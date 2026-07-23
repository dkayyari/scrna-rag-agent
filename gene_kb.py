"""
Step 2: Gene annotation knowledge base for RAG retrieval.

Two-tier lookup:
  1. Local curated dictionary — canonical PBMC marker genes, for fast/offline
     lookup. Compiled from standard immunology references (UniProt, NCBI Gene)
     and the marker gene set used in the Scanpy PBMC3k tutorial.
  2. Live fallback via the mygene.info API (NCBI-derived) for any gene not
     in the local dictionary — covers the much larger set of ~13,000+ genes
     that survive QC filtering, at the cost of needing internet access and
     being slower than the local lookup.

Known limitation: the local KB covers only ~17 well-known marker genes out
of thousands present in the dataset. This is a deliberate scope choice, not
an oversight — see README for discussion.
"""

import requests

LOCAL_GENE_KB = {
    "CD3D": "T cell receptor complex component; canonical marker for T cells.",
    "CD3E": "T cell receptor complex component; canonical marker for T cells.",
    "CD2": "T cell surface marker, involved in T cell adhesion and activation.",
    "IL7R": "Marker for naive/memory CD4+ T cells.",
    "IL32": "Cytokine produced by activated T cells.",
    "LTB": "Lymphotoxin beta; expressed in lymphocytes, involved in immune signaling.",
    "CD14": "Marker for classical monocytes; LPS co-receptor.",
    "LYZ": "Lysozyme; highly expressed in monocytes/macrophages.",
    "MS4A1": "Marker for B cells (CD20).",
    "CD79A": "B cell receptor complex component; B cell marker.",
    "NKG7": "Marker for NK cells and cytotoxic lymphocytes.",
    "GNLY": "Granulysin; expressed in NK cells and cytotoxic T cells.",
    "FCGR3A": "Marker for non-classical (CD16+) monocytes and NK cells.",
    "PPBP": "Platelet basic protein; marker for platelets/megakaryocytes.",
    "PF4": "Platelet factor 4; marker for platelets/megakaryocytes.",
    "FCER1A": "Marker for dendritic cells.",
    "CD8A": "Marker for cytotoxic (CD8+) T cells.",
}


def retrieve_annotation(gene_symbol: str, use_live_lookup: bool = True) -> dict:
    """Retrieve an annotation for a gene symbol. Returns gene/summary/source."""
    gene_symbol = gene_symbol.upper()

    if gene_symbol in LOCAL_GENE_KB:
        return {
            "gene": gene_symbol,
            "summary": LOCAL_GENE_KB[gene_symbol],
            "source": "local curated KB",
        }

    if use_live_lookup:
        try:
            resp = requests.get(
                "https://mygene.info/v3/query",
                params={"q": f"symbol:{gene_symbol}", "fields": "summary,name", "species": "human"},
                timeout=5,
            )
            hits = resp.json().get("hits", [])
            if hits and "summary" in hits[0]:
                return {
                    "gene": gene_symbol,
                    "summary": hits[0]["summary"][:500],
                    "source": "mygene.info (NCBI-derived)",
                }
        except Exception as e:
            print(f"Live lookup failed for {gene_symbol}: {e}")

    return {
        "gene": gene_symbol,
        "summary": "No annotation found in local KB or live lookup.",
        "source": "none",
    }


if __name__ == "__main__":
    print(retrieve_annotation("CD3D"))
    print(retrieve_annotation("LDHB"))
