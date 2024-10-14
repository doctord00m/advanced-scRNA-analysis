import scanpy as sc

# Load and preprocess your data as before
adata = sc.read_10x_h5("5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, flavor="igraph", directed=False, n_iterations=2)

# Save the processed data
adata.write("processed_pbmc_data.h5ad")
