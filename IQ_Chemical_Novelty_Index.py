import os
import logging
import pandas as pd
import numpy as np
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from typing import Dict, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Constants
WEIGHTS = (0.33, 0.33, 0.33)  # (Ans, Hcos, Hmqs)
MQS_STRATA = [(i/100, (i+5)/100) for i in range(0, 100, 5)]
COS_STRATA = [(i/100, (i+5)/100) for i in range(70, 100, 5)]


def load_data(graph_path: str, metadata_path: str) -> pd.DataFrame:
    """Load and merge graph data with metadata"""
    G = nx.read_graphml(graph_path, node_type=str)
    nodes = (
        pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient='index')
        .reset_index()
        .rename(columns={
            'index': 'id',
            'ATTRIBUTE_Strain': 'ATTRIBUTE_Species',
            'ATTRIBUTE_Bacteria': 'ATTRIBUTE_Genus'
        })
    )
    for col in ['MQScore', 'cosine_score']:
        if col in nodes.columns:
            nodes[col] = pd.to_numeric(nodes[col], errors='coerce')

    metadata = pd.read_csv(
        metadata_path,
        sep=';',
        dtype={'ATTRIBUTE_Species': 'category', 'ATTRIBUTE_Genus': 'category'}
    ).rename(columns={
        'ATTRIBUTE_Strain': 'ATTRIBUTE_Species',
        'ATTRIBUTE_Bacteria': 'ATTRIBUTE_Genus'
    })

    merged = nodes.merge(
        metadata,
        on=['ATTRIBUTE_Species', 'ATTRIBUTE_Genus'],
        how='left'
    )
    return merged


def calculate_strata_counts(
    df: pd.DataFrame,
    value_col: str,
    strata: list,
    unique_id: bool = True
) -> pd.DataFrame:
    """Count observations per species across strata"""
    df = df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors='coerce')
    exploded = df.assign(
        ATTRIBUTE_Species=df['ATTRIBUTE_Species'].str.split(',')
    ).explode('ATTRIBUTE_Species')

    bins = [s[0] for s in strata] + [strata[-1][1]]
    labels = [f"{s[0]:.2f}-{s[1]:.2f}" for s in strata]

    mask = exploded[value_col].notna()
    exploded.loc[mask, 'stratum'] = pd.cut(
        exploded.loc[mask, value_col], bins=bins, labels=labels, include_lowest=True
    )

    if unique_id:
        counts = (
            exploded.drop_duplicates('id')
            .groupby(['ATTRIBUTE_Species', 'stratum'], observed=True)
            .size()
        )
    else:
        counts = (
            exploded.groupby(['ATTRIBUTE_Species', 'stratum'], observed=True)
            .size()
        )
    return counts.unstack(fill_value=0)


def calculate_shannon_index(counts: pd.DataFrame) -> pd.Series:
    """Compute Shannon diversity index for each species"""
    total = counts.sum(axis=1)
    p = counts.div(total, axis=0)
    p_nonzero = p.where(p > 0)
    return -(p_nonzero * np.log(p_nonzero)).sum(axis=1)


def process_cosine_scores(edges: pd.DataFrame, nodes: pd.DataFrame) -> pd.DataFrame:
    """Extract and annotate cosine scores per species"""
    if 'source' in edges.columns and 'target' in edges.columns:
        src, tgt = 'source', 'target'
    else:
        src, tgt = 'node1', 'node2'

    adj = pd.concat([
        edges[[src, tgt, 'cosine_score']],
        edges.rename(columns={src: tgt, tgt: src})[[src, tgt, 'cosine_score']]
    ])

    merged = adj.merge(
        nodes[['id', 'ATTRIBUTE_Species']],
        left_on=src,
        right_on='id',
        how='left'
    )
    return merged[[src, 'cosine_score', 'ATTRIBUTE_Species']].rename(
        columns={src: 'node', 'cosine_score': 'Cosine_score'}
    )


def integrate_npatlas(
    score_df: pd.DataFrame,
    npatlas_url: str
) -> pd.DataFrame:
    """Merge NPAtlas compound info into score dataframe"""
    npa = pd.read_excel(npatlas_url, engine='openpyxl')[[
        'compound_name', 'compound_inchi', 'compound_smiles',
        'origin_type', 'genus', 'origin_species'
    ]].rename(columns={
        'genus': 'ATTRIBUTE_Genus',
        'origin_species': 'ATTRIBUTE_Species'
    })

    npa['Canonical_Smiles'] = npa['compound_smiles'].apply(
        lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)) if pd.notnull(s) else None
    )

    merged = score_df.merge(
        npa,
        left_on=['Smiles', 'INCHI', 'Compound_Name'],
        right_on=['compound_smiles', 'compound_inchi', 'compound_name'],
        how='inner'
    )
    return merged.dropna(subset=['origin_type']).drop_duplicates()


def calculate_iq_components(
    total_counts: pd.Series,
    annotated_counts: pd.Series,
    shannon_cos: pd.Series,
    shannon_mqs: pd.Series,
    weights: Tuple[float, float, float] = WEIGHTS
) -> pd.DataFrame:
    """Combine Ans, Hcos, Hmqs into IQ score"""
    ans = annotated_counts.div(total_counts).fillna(0)
    hcos = shannon_cos.div(shannon_cos.max()).fillna(0)
    hmqs = shannon_mqs.div(shannon_mqs.max()).fillna(0)

    w_ans, w_hcos, w_hmqs = weights
    iq = w_ans*ans + w_hcos*hcos + w_hmqs*hmqs

    return pd.DataFrame({
        'Ans': ans,
        'Hcos': hcos,
        'Hmqs': hmqs,
        'IQ': iq
    })


def generate_plots(iq_df: pd.DataFrame, output_dir: str):
    """Create and save plots"""
    os.makedirs(output_dir, exist_ok=True)

    # Bar plot of IQ by species
    plt.figure(figsize=(15, 8))
    sns.barplot(
        data=iq_df.sort_values('IQ', ascending=False),
        x='ATTRIBUTE_Species', y='IQ', hue='ATTRIBUTE_Genus', palette='viridis'
    )
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'IQ_barplot.png'))
    plt.close()

    # Distribution of IQ by genus
    plt.figure(figsize=(15, 8))
    sns.boxplot(data=iq_df, x='ATTRIBUTE_Genus', y='IQ')  # palette removed to avoid future warning
    sns.swarmplot(data=iq_df, x='ATTRIBUTE_Genus', y='IQ', color='black', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'IQ_distribution.png'))
    plt.close()

    # Stacked parameter contribution per species
    params = iq_df[['Ans', 'Hcos', 'Hmqs']].copy()
    # species are in the index already
    params_norm = params.div(params.sum(axis=1), axis=0)
    params_norm.plot.bar(stacked=True, figsize=(15, 8), colormap='Set2')
    plt.xlabel('ATTRIBUTE_Species')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'parameter_contribution.png'))
    plt.close()


def main():
    GRAPH_PATH = 'Actino_network.graphml'
    METADATA_PATH = 'Metadata_Actinomyces.csv'
    OUTPUT_DIR = 'Results'
    NPATLAS_URL = 'https://www.npatlas.org/static/downloads/NPAtlas_download.xlsx'

    logging.info("Loading data...")
    nodes = load_data(GRAPH_PATH, METADATA_PATH)

    G = nx.read_graphml(GRAPH_PATH, node_type=str)
    edges = nx.to_pandas_edgelist(G)

    total_counts = (
        nodes['ATTRIBUTE_Species']
        .str.split(',')
        .explode()
        .value_counts()
    )

    mqs_nodes = nodes.dropna(subset=['MQScore'])
    mqs_counts = calculate_strata_counts(mqs_nodes, 'MQScore', MQS_STRATA)
    shannon_mqs = calculate_shannon_index(mqs_counts)

    cos_df = process_cosine_scores(edges, nodes)
    cos_counts = calculate_strata_counts(cos_df, 'Cosine_score', COS_STRATA, unique_id=False)
    shannon_cos = calculate_shannon_index(cos_counts)

    annotated_counts = (
        mqs_nodes['ATTRIBUTE_Species'].str.split(',').explode().value_counts()
    )

    iq_df = calculate_iq_components(total_counts, annotated_counts, shannon_cos, shannon_mqs)

    species_genus = nodes[['ATTRIBUTE_Species', 'ATTRIBUTE_Genus']].drop_duplicates()
    plot_df = iq_df.join(species_genus.set_index('ATTRIBUTE_Species'))

    logging.info("Generating plots...")
    generate_plots(plot_df, OUTPUT_DIR)

    logging.info("Done.")


if __name__ == "__main__":
    main()