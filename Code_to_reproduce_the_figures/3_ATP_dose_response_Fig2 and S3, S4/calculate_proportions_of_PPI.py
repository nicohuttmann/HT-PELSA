# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import os
import json
from tqdm import tqdm
import scipy.stats as stats
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import networkx as nx
import itertools

# %%
def download_protein_info(protein_ids, uniprot_json = True, alphafold_pdb = True):

    for int_prot in tqdm(protein_ids):
        if uniprot_json:
            dest_json = f'data/protein_data/{int_prot}.json'
            if not os.path.exists(dest_json):
                url = f'https://www.uniprot.org/uniprotkb/{int_prot}.json'
                r = requests.get(url)
                with open(dest_json, 'w') as f:
                    f.write(r.text)

        if alphafold_pdb:
            dest_pdb = f'data/protein_data/{int_prot}.pdb'
            if not os.path.exists(dest_pdb):
                url = f'https://alphafold.ebi.ac.uk/files/AF-{int_prot}-F1-model_v4.pdb'
                r = requests.get(url)
                with open(dest_pdb, 'w') as f:
                    f.write(r.text)

def get_interpro_table(protein_ids):
    interpro_df = []
    for int_prot in protein_ids:
        input_json = f'data/protein_data/{int_prot}.json'
        with open(input_json, 'r') as f:
            data = json.load(f)
        for xref in data.get('uniProtKBCrossReferences', []):
            if xref.get('database') == 'InterPro':
                protdf = pd.DataFrame({
                    'uniprot_id': int_prot,
                    'interpro_id': xref.get('id', ''),
                    'description': xref.get('properties', [{}])[0].get('value', '')
                }, index=[0])
                interpro_df.append(protdf)
    interpro_df = pd.concat(interpro_df, ignore_index=True)
    return interpro_df

def analyze_interpro_enrichment_by_category(df_reg, df_annot, category, min_proteins=3, top_n=20, adjust_method='fdr_bh'):
    """
    Performs enrichment analysis for InterPro domains in a given regulation category.

    Parameters:
    - df_reg: DataFrame with 'uniprot_id' and 'regulation' columns.
    - df_annot: DataFrame with 'uniprot_id', 'interpro_id', 'description'.
    - category: Regulation category to analyze (e.g., 'stabilized' or 'destabilized').
    - min_proteins: Minimum number of category proteins required per domain.
    - top_n: Number of top enriched domains to plot.
    - adjust_method: Method for multiple testing correction (default 'fdr_bh').

    Returns:
    - DataFrame with enrichment stats.
    """
    # Handle multiple IDs separated by semicolon
    df_reg = df_reg.copy()
    df_reg['uniprot_id'] = df_reg['uniprot_id'].str.split(';')
    df_reg = df_reg.explode('uniprot_id')
    
    # Merge annotation with regulation table
    merged = df_annot.merge(df_reg, on='uniprot_id', how='inner')
    
    # Define sets
    all_proteins = set(df_reg['uniprot_id'])
    category_proteins = set(df_reg[df_reg['regulation'] == category]['uniprot_id'])

    # Group by InterPro domain and calculate stats
    results = []
    for interpro_id, group in merged.groupby('interpro_id'):
        proteins_with_domain = set(group['uniprot_id'])
        with_domain_and_category = len(proteins_with_domain & category_proteins)
        if with_domain_and_category < min_proteins:
            continue

        a = with_domain_and_category
        b = len(category_proteins - proteins_with_domain)
        c = len(proteins_with_domain - category_proteins)
        d = len(all_proteins - category_proteins - proteins_with_domain)

        contingency = [[a, b], [c, d]]
        odds_ratio, p_value = stats.fisher_exact(contingency, alternative='greater')
        description = group['description'].iloc[0]

        results.append({
            'interpro_id': interpro_id,
            'description': description,
            f'{category}_with_domain': a,
            'total_with_domain': len(proteins_with_domain),
            'odds_ratio': odds_ratio,
            'p_value': p_value
        })

    df_results = pd.DataFrame(results)

    if df_results.empty:
        print(f"No enriched domains found for category '{category}'.")
        return df_results

    # Adjust p-values
    df_results['p_adj'] = multipletests(df_results['p_value'], method=adjust_method)[1]
    df_results = df_results.sort_values('p_adj')

    return df_results

def plot_enrichment_significance(input_dict, top_n=10, figsize=(10, 4), pval_thresh=0.05):
    """
    Create subplots where bars represent enrichment significance (-log10 adjusted p-value).

    Parameters:
    - input_dict: dict of {label: enrichment dataframe}.
    - top_n: number of top terms per category to plot.
    - figsize: base figure size (per subplot).
    - pval_thresh: FDR threshold for significance marker.
    """

    n_panels = len(input_dict)
    fig, axes = plt.subplots(1, n_panels, figsize=figsize, constrained_layout=True)
    if n_panels == 1:
        axes = [axes]

    for ax, (label, df) in zip(axes, input_dict.items()):
        if df.empty:
            ax.axis('off')
            ax.text(0.5, 0.5, f"No enriched terms for '{label}'", ha='center', va='center')
            continue

        df = df.copy()
        df['neg_log10_padj'] = -np.log10(df['p_adj'].replace(0, np.nan))
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['neg_log10_padj'])

        # Select top N most significant terms
        df = df.sort_values('p_adj').head(top_n)

        # Plot
        sns.barplot(
            data=df,
            x='neg_log10_padj',
            y='description',
            color='crimson',
            ax=ax
        )

        # Annotate each bar with p_adj inside
        for bar, pval in zip(ax.patches, df['p_adj']):
            width = bar.get_width()
            y = bar.get_y() + bar.get_height() / 2
            text = f"{pval:.2e}"
            ax.text(
                width - 0.05 * width, y,
                text,
                ha='right',
                va='center',
                color='white' if width > 1.5 else 'black',
                fontsize=9
            )

        ax.axvline(-np.log10(pval_thresh), color='black', linestyle='--', linewidth=1)
        ax.set_title(f"{label.capitalize()} proteins", fontsize=13, weight='bold')
        ax.set_xlabel(r"$-\log_{10}$(adjusted p-value)")
        ax.set_ylabel("InterPro domain")

    plt.show()
    return fig, axes

def compute_category_distances(df_proteins, ppi_graph):
    """
    Computes shortest path distances between protein pairs grouped by unordered regulation categories.

    Parameters:
    - df_proteins: DataFrame with 'uniprot_id' and 'regulation'.
    - ppi_graph: networkx.Graph object (undirected, unweighted).

    Returns:
    - DataFrame with columns: protein_A, protein_B, category_pair (tuple), distance
    """
    df = df_proteins[['uniprot_id', 'regulation']].drop_duplicates()
    category_map = df.groupby('regulation')['uniprot_id'].apply(list).to_dict()

    results = []
    unique_category_pairs = list(itertools.combinations_with_replacement(category_map.keys(), 2))

    for catA, catB in tqdm(unique_category_pairs, desc="Computing category distances"):
        protA = category_map[catA]
        protB = category_map[catB]

        for a in protA:
            for b in protB:
                if a == b:
                    continue  # avoid self-pairs
                if not (a in ppi_graph and b in ppi_graph):
                    dist = None
                else:
                    try:
                        dist = nx.shortest_path_length(ppi_graph, source=a, target=b)
                    except nx.NetworkXNoPath:
                        dist = None

                # Enforce undirected category pair
                pair = tuple(sorted([catA, catB]))

                results.append({
                    'protein_A': a,
                    'protein_B': b,
                    'category_A': catA,
                    'category_B': catB,
                    'category_pair': pair,
                    'distance': dist
                })

    return pd.DataFrame(results)


    df_plot = df_distances.dropna(subset=['distance']).copy()
    df_plot['category_label'] = df_plot['category_pair'].apply(lambda x: f"{x[0]} ↔ {x[1]}")

    mean_distances = df_plot.groupby('category_label')['distance'].mean().reset_index()

    # arrangein decreasing order
    mean_distances = mean_distances.sort_values(by='distance', ascending=False).reset_index(drop=True)

    plt.figure(figsize=figsize)
    ax = sns.barplot(data=mean_distances, x='category_label', y='distance', color='steelblue')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Mean Shortest Path Distance")
    plt.xlabel("Category Pair")
    plt.title("Average Network Distance Between Protein Categories")
    plt.tight_layout()
    plt.show()

def compute_direct_interaction_proportions(df_proteins, ppi_graph):

    # Clean input data and group by category
    df = df_proteins[['uniprot_id', 'regulation']].drop_duplicates()
    category_map = df.groupby('regulation')['uniprot_id'].apply(list).to_dict()

    # All unique category pairs, including self-comparisons
    unique_category_pairs = list(itertools.combinations_with_replacement(category_map.keys(), 2))

    results = []

    for catA, catB in tqdm(unique_category_pairs, desc="Computing direct interaction proportions"):
        protA = category_map[catA]
        protB = category_map[catB]

        total = 0
        connected = 0

        if catA == catB:
            # Use combinations to avoid redundant and self-pairs
            for a, b in itertools.combinations(protA, 2):
                if a in ppi_graph and b in ppi_graph:
                    total += 1
                    if ppi_graph.has_edge(a, b):
                        connected += 1
        else:
            # Cross-category pairs — no redundancy
            for a in protA:
                for b in protB:
                    if a in ppi_graph and b in ppi_graph:
                        total += 1
                        if ppi_graph.has_edge(a, b):
                            connected += 1

        pair = tuple(sorted([catA, catB]))
        prop = connected / total if total > 0 else None

        results.append({
            'category_pair': pair,
            'n_total': total,
            'n_connected': connected,
            'proportion_connected': prop
        })

    return pd.DataFrame(results)

def fisher_comparisons(results_df, baseline_label='unresponsive ↔ unresponsive'):
    # Create a lookup from category_pair to stats
    stats = results_df.set_index('category_pair')[['n_total', 'n_connected']].to_dict('index')
    
    # Target comparisons
    test_labels = [
        'destabilized ↔ destabilized',
        'destabilized ↔ stabilized',
        'stabilized ↔ stabilized'
    ]
    
    # Baseline
    base = stats[baseline_label]
    base_connected = base['n_connected']
    base_not_connected = base['n_total'] - base_connected

    # Store results
    fisher_results = []

    for label in test_labels:
        row = stats[label]
        a = row['n_connected']
        b = row['n_total'] - row['n_connected']
        c = base_connected
        d = base_not_connected

        contingency = [[a, b], [c, d]]
        oddsratio, pvalue = fisher_exact(contingency, alternative='greater')  # One-sided

        fisher_results.append({
            'test_pair': label,
            'connected': a,
            'not_connected': b,
            'baseline_connected': c,
            'baseline_not_connected': d,
            'odds_ratio': oddsratio,
            'p_value': pvalue
        })

    return pd.DataFrame(fisher_results)


# %%
data = pd.read_csv('data/ATP_PELSA.csv')
filt_data = data[['Protein.Group', 'regulation', 'mean_EC50', 'CV_EC50']].copy()

# create output directotry in results
if not os.path.exists(f'results_r1/'):
    os.makedirs(f'results_r1/')

# here filtering of peptides
filt_data = filt_data[((filt_data['regulation'] == 'stabilized') & (filt_data['mean_EC50'] >= 3) & (filt_data['CV_EC50'] <= 0.25)) | ((filt_data['regulation'] == 'destabilized') & (filt_data['mean_EC50'] >= 3) & (filt_data['CV_EC50'] <= 0.25)) | (filt_data['regulation'] == 'unresponsive')]

filt_data = filt_data[filt_data['regulation'] != 'underquantified']
filt_data = filt_data.groupby('Protein.Group').agg({
    'regulation': lambda x: 'destabilized' if 'destabilized' in x.values else ('stabilized' if 'stabilized' in x.values else 'unresponsive'),
    'mean_EC50': 'mean',
    'CV_EC50': 'mean'
}).reset_index()
filt_data['uniprot_id'] = filt_data['Protein.Group'].str.split(';')
filt_data = filt_data.explode('uniprot_id')
filt_data = filt_data[['uniprot_id', 'regulation', 'mean_EC50', 'CV_EC50']].copy().drop_duplicates()
# for stabilized and destabilized, apply a filter which is mean_EC50 >= 3 and mean CV_EC50 <= 0.25

filt_data['regulation'].value_counts()

# %%
all_unique_proteins = filt_data['uniprot_id'].unique()
print('Analyzing {} proteins'.format(len(all_unique_proteins)))
download_protein_info(all_unique_proteins)
interpro_df = get_interpro_table(all_unique_proteins)

# %%
stabilized_results = analyze_interpro_enrichment_by_category(filt_data, interpro_df, category='stabilized')
# save to results 
stabilized_results.to_csv(f'results_r1//stabilized_interpro_enrichment.csv', index=False)
destabilized_results = analyze_interpro_enrichment_by_category(filt_data, interpro_df, category='destabilized')
# save to results
destabilized_results.to_csv(f'results_r1//destabilized_interpro_enrichment.csv', index=False)
res_dict = {'stabilized': stabilized_results, 'destabilized': destabilized_results}
fig,axs = plot_enrichment_significance(res_dict, top_n=10)
# save figures
fig.savefig(f'results_r1/interpro_enrichment.pdf', bbox_inches='tight')

# %% [markdown]
# ## Network analysis

# %%
biogrid = pd.read_csv('data/BIOGRID-ALL-4.4.247.tab3.txt', sep='\t', low_memory=False)
ecoli_biogrid = biogrid[(biogrid['Organism Name Interactor A'] == 'Escherichia coli (K12/W3110)') & ( biogrid['Organism Name Interactor B'] == 'Escherichia coli (K12/W3110)')]
ecoli_biogrid = ecoli_biogrid[['SWISS-PROT Accessions Interactor A', 'SWISS-PROT Accessions Interactor B']].rename(
    columns={
        'SWISS-PROT Accessions Interactor A': 'uniprot_id_A',
        'SWISS-PROT Accessions Interactor B': 'uniprot_id_B'
    })
ecoli_biogrid = ecoli_biogrid.drop_duplicates()
ecoli_biogrid = ecoli_biogrid[(ecoli_biogrid['uniprot_id_A'] != '-') & (ecoli_biogrid['uniprot_id_B'] != '-')]

# %%
ecoli_net = nx.from_pandas_edgelist(
    ecoli_biogrid,
    source='uniprot_id_A',
    target='uniprot_id_B',
    create_using=nx.Graph()
)

# %%
# count how many proteins are in the network
print(f"Number of proteins in the network: {ecoli_net.number_of_nodes()}")
print(f"Number of interactions in the network: {ecoli_net.number_of_edges()}")

# check hwo many proteins from filt_data are in the network
filt_proteins = set(filt_data['uniprot_id'].unique())
proteins_in_network = set(ecoli_net.nodes())
common_proteins = filt_proteins.intersection(proteins_in_network)
print(f"Number of proteins in filt_data: {len(filt_proteins)}")
print(f"Number of proteins from filt_data in the network: {len(common_proteins)}")

# %%
ppi_result = compute_category_distances(filt_data, ecoli_net)
ppi_summarized = ppi_result.dropna(subset=['distance']).groupby('category_pair').agg({
    'distance': 'mean'
}).reset_index()
ppi_summarized['category_pair'] = ppi_summarized['category_pair'].apply(lambda x: f"{x[0]} ↔ {x[1]}")
ppi_summarized.to_csv(f'results_r1//category_distances.csv', index=False)

# %%
sns.barplot(
    data=ppi_summarized,
    x='category_pair',
    y='distance',
    color='steelblue'
)
# save
plt.xticks(rotation=45, ha='right')
plt.ylabel("Mean Shortest Path Distance")
plt.xlabel("Category Pair")
plt.title("Average Network Distance Between Protein Categories")
plt.tight_layout()
plt.savefig(f'results_r1//category_distances.pdf', bbox_inches='tight')

# %%
ppi_proportions = compute_direct_interaction_proportions(filt_data, ecoli_net)
ppi_proportions['category_pair'] = ppi_proportions['category_pair'].apply(lambda x: f"{x[0]} ↔ {x[1]}")
ppi_proportions.to_csv(f'results_r1//direct_interaction_proportions.csv', index=False)
# create plot 
sns.barplot(
    data=ppi_proportions,
    x='category_pair',
    y='proportion_connected',
    color='mediumseagreen'
)
# modify pplot and save
plt.xticks(rotation=45, ha='right')
plt.ylabel("Proportion of Direct Interactions\nin BIOGRID")
plt.xlabel("Category Pair")
plt.title("Direct PPI Connectivity Between Protein Categories")
plt.ylim(0, ppi_proportions['proportion_connected'].max() * 1.2)
plt.tight_layout()
plt.savefig(f'results_r1//direct_interaction_proportions.pdf', bbox_inches='tight')

# %%
# compute fisher exact test for proportions
fisher_results = fisher_comparisons(ppi_proportions, baseline_label='unresponsive ↔ unresponsive')
fisher_results['p_adj'] = multipletests(fisher_results['p_value'], method='fdr_bh')[1]
fisher_results.to_csv(f'results_r1//interaction_fisher_results.csv', index=False)








