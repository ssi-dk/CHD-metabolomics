# filter_IIMN_adducts.py
#
# Filters adducts in a feature table identified using IIMN
# (Ion Identity Molecular Networking).
#
# For each feature (presumed metabolite) with multiple identified adducts,
# selects the adduct that:
#   1) Is observed in the most samples (signal > 0)
#   2) If tied, has the highest median signal

# Import necessary libraries
import pandas as pd
import networkx as nx

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths — edit these before running
# ---------------------------------------------------------------------------

feature_table_file        = Path("path/to/feature_table.csv")
ion_identity_network_file = Path("path/to/edges_msannotation.csv")
output_file               = Path("path/to/feature_table_IIMN_adduct_filtered.csv")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

ft_df    = pd.read_csv(feature_table_file, index_col=0)
edges_df = pd.read_csv(ion_identity_network_file)

# ---------------------------------------------------------------------------
# Save to current directory
# ---------------------------------------------------------------------------

ft_df.to_csv("original_feature_table.csv")
edges_df.to_csv("edges.csv")

# ---------------------------------------------------------------------------
# Get sample columns
# ---------------------------------------------------------------------------

sample_columns = [column for column in ft_df.columns if "Peak area" in column]
sample_names   = [sample_name.split(".")[0] for sample_name in sample_columns]

# ---------------------------------------------------------------------------
# Create graph object from edges
# ---------------------------------------------------------------------------

# Extracting the column names (excluding the first two) from edges_df to identify edge attributes.
edge_attributes = edges_df.columns[2:]

# Creating a graph from the edges DataFrame using NetworkX.
# The graph is created by treating 'ID1' and 'ID2' columns as source and target nodes respectively.
adduct_graph = nx.from_pandas_edgelist(edges_df, source="ID1", target="ID2")

# This loop iterates through each column in the feature table DataFrame (ft_df).
# - For each column, a dictionary mapping the node index to the column's value is created.
# - This dictionary is then used to set node attributes in the graph with the column name as the attribute key.
for column in ft_df.columns:
    nx.set_node_attributes(adduct_graph, pd.Series(ft_df[column], index=ft_df.index).to_dict(), name=column)

# ---------------------------------------------------------------------------
# Create list of connected adducts (sets)
# ---------------------------------------------------------------------------

# Identify connected components (adducts belonging to the same parent metabolite ) in the graph:
connected_adduct_sets = sorted(nx.connected_components(adduct_graph), key=len, reverse=True)

print(f"Number of connected adduct sets: {len(connected_adduct_sets)}")

avg_size = sum([len(s) for s in connected_adduct_sets]) / len(connected_adduct_sets)
print(f"Average size of connected adduct sets: {avg_size:.2f}")

# ---------------------------------------------------------------------------
# Adduct filtering
# ---------------------------------------------------------------------------

# Keep the adduct with the max number of samples with non-zero intensities. If several adducts have the same max, then choose the one with the highest median peak area. 
# In the original algorithm, if several adducts have the same number of max observationes, then the adduct with the first row label was selected. 

remove_feature_indices = []

for connected_adducts in connected_adduct_sets:
    # Select peak areas for connected adducts
    sel_mask = ft_df.index.isin(connected_adducts)
    ft_connected_adducts = ft_df.loc[sel_mask, sample_columns]

    # Find how many samples have non-zero values for each adduct in the set
    adducts_non_zero_samples_counts = ft_connected_adducts.gt(0).sum(axis=1)
    adducts_non_zero_samples_ids    = adducts_non_zero_samples_counts.index.to_list()

    # Find maximum number of samples with non-zero values
    max_non_zero_counts = adducts_non_zero_samples_counts.max()

    # Select the adduct with the maximum number of samples with non-zero values
    adduct_id_max = adducts_non_zero_samples_counts.idxmax()

    # If tied, select the adduct with the highest median value among candidates
    if sum(candidates := (adducts_non_zero_samples_counts == max_non_zero_counts)) > 1:
        adduct_id_max = ft_connected_adducts[candidates].median(axis=1).idxmax()

    # Remove the adduct to keep from the removal list
    remove_adducts_ids = adducts_non_zero_samples_ids
    remove_adducts_ids.remove(adduct_id_max)

    remove_feature_indices = remove_feature_indices + remove_adducts_ids

ft_df_filtered = ft_df.drop(remove_feature_indices)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

n_removed = ft_df.shape[0] - ft_df_filtered.shape[0]
print(
f"""
    The number of features in the original feature table: {ft_df.shape[0]}
    The number of features in the filtered feature table: {ft_df_filtered.shape[0]}
    The number of features removed: {n_removed} ({n_removed/ft_df.shape[0] * 100:.2f}%)
"""
)

# ---------------------------------------------------------------------------
# Export filtered feature table
# ---------------------------------------------------------------------------

ft_df_filtered.to_csv(output_file)
print(f"Filtered feature table saved to: {output_file}")