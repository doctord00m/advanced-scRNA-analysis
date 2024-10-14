import streamlit as st
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load preprocessed data
st.title("Advanced scRNA-seq Analysis App")

@st.cache_data
def load_data():
    adata = sc.read_h5ad("processed_pbmc_data.h5ad")
    return adata

adata = load_data()

# Sidebar: User options
st.sidebar.header("Clustering Options")
resolution = st.sidebar.slider("Select clustering resolution", 0.1, 2.0, 1.0, 0.1)
st.sidebar.header("Gene Expression Options")
gene = st.sidebar.text_input("Enter a gene to visualize", "CD3D")
cluster = st.sidebar.selectbox("Select a cluster for differential expression", adata.obs['leiden'].cat.categories)

# Apply clustering with selected resolution
sc.tl.leiden(adata, resolution=resolution, flavor="igraph", directed=False, n_iterations=2)
adata.obs['leiden_res'] = adata.obs['leiden']

# Main visualization
st.write("## UMAP Visualization with Selected Clustering Resolution")
fig, ax = plt.subplots()
sc.pl.umap(adata, color=['leiden_res'], ax=ax, show=False)
st.pyplot(fig)

# Gene Expression Visualization
st.write(f"## Expression of {gene} Across Clusters")
fig, ax = plt.subplots()
sc.pl.umap(adata, color=gene, ax=ax, show=False)
st.pyplot(fig)

# Differential Expression Analysis
st.write(f"## Top Differentially Expressed Genes in Cluster {cluster}")
sc.tl.rank_genes_groups(adata, 'leiden_res', groups=[cluster], reference='rest', method='wilcoxon')
de_df = sc.get.rank_genes_groups_df(adata, group=cluster)
st.dataframe(de_df.head(10))

# Heatmap of Marker Genes
st.write("## Heatmap of Top Marker Genes Across Clusters")
top_genes = de_df['names'].head(10).values
sc.pl.heatmap(adata, var_names=top_genes, groupby='leiden_res', show_gene_labels=True, cmap='viridis')
st.pyplot(plt)
