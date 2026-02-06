import networkx as nx
import warnings
from .constants import STRING_API_URL, STRING_SPECIES, STRING_MIN_SCORE


def fetch_string_ppi(genes, species=STRING_SPECIES, required_score=STRING_MIN_SCORE):
    try:
        import requests
    except ImportError:
        raise ImportError("requests is required for STRING API. Install with: pip install requests")
    G = nx.Graph()
    G.add_nodes_from(genes)
    if not genes:
        return G
    request_url = f"{STRING_API_URL}/tsv/network"
    params = {
        "identifiers": "\r".join(genes),
        "species": species,
        "required_score": required_score,
        "caller_identity": "bioinformatics_pipeline_pmc11805404",
    }
    try:
        resp = requests.post(request_url, data=params, timeout=60)
        resp.raise_for_status()
    except Exception as e:
        warnings.warn(f"STRING API request failed: {e}. PPI network will have no edges.")
        return G
    lines = resp.text.strip().split("\n")
    if not lines or len(lines) < 2:
        return G
    header = lines[0].split("\t")
    if "preferredName_A" not in lines[0]:
        return G
    try:
        idx_a = header.index("preferredName_A")
        idx_b = header.index("preferredName_B")
        idx_score = header.index("score")
    except ValueError:
        idx_a, idx_b, idx_score = 2, 3, 5
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) <= max(idx_a, idx_b, idx_score):
            continue
        name_a = parts[idx_a].strip()
        name_b = parts[idx_b].strip()
        if name_a in genes and name_b in genes and name_a != name_b:
            try:
                score = float(parts[idx_score])
            except (ValueError, IndexError):
                score = 0.4
            G.add_edge(name_a, name_b, weight=score)
    return G


def construct_ppi_network(genes, k=10):
    print(f"\n{'='*70}")
    print("Constructing PPI Network (STRING database)")
    print(f"{'='*70}")
    print(f"  Species: Homo sapiens (9606), min confidence: 0.4 (highest)")
    G = fetch_string_ppi(genes)
    degree_centrality = nx.degree_centrality(G)
    sorted_genes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)
    core_genes = [gene for gene, _ in sorted_genes[:k]]
    if len(core_genes) < k:
        for g in genes:
            if g not in core_genes:
                core_genes.append(g)
                if len(core_genes) >= k:
                    break
        core_genes = core_genes[:k]
    print(f"✓ PPI network from STRING")
    print(f"  - Nodes: {G.number_of_nodes()}")
    print(f"  - Edges: {G.number_of_edges()}")
    print(f"  - Top {len(core_genes)} core genes (by degree centrality)")
    print(f"\nTop {len(core_genes)} Core Genes (by degree centrality):")
    for i, (gene, centrality) in enumerate(sorted_genes[:k], 1):
        print(f"  {i}. {gene}: {centrality:.3f}")
    if G.number_of_edges() == 0:
        print("  (No interactions returned by STRING for these genes; core list is first k genes.)")
    return G, core_genes
