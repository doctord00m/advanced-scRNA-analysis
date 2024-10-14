---

# Advanced scRNA-seq Analysis App

Welcome to the Advanced scRNA-seq Analysis App! This app leverages Streamlit and Scanpy to provide an interactive and insightful analysis platform for single-cell RNA sequencing (scRNA-seq) data. 

## Features

- **Dynamic Clustering**: Adjust clustering resolution on the fly using the Leiden algorithm.
- **Gene Expression Visualization**: Input any gene to visualize its expression across clusters on a UMAP plot.
- **Differential Expression Analysis**: Identify top marker genes for specific clusters.
- **Customizable Heatmap**: Generate a heatmap of top marker genes to explore expression patterns across clusters.

## Visualizations

### 1. UMAP Visualization with Selected Clustering Resolution
Adjust the clustering resolution to explore how clusters change. This UMAP plot visualizes clusters based on the selected resolution.

![UMAP Visualization](./path_to_images/umap_visualization.png)

### 2. Expression of CD3D Across Clusters
Input a gene (e.g., CD3D) and visualize its expression across clusters on the UMAP plot.

![CD3D Gene Expression](./path_to_images/cd3d_expression.png)

### 3. Heatmap of Top Marker Genes Across Clusters
Explore the top differentially expressed genes with this heatmap, which shows marker genes across clusters based on the selected clustering resolution.

![Heatmap of Marker Genes](./path_to_images/heatmap_marker_genes.png)

## How to Run

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your_username/Advanced-scRNA-seq-Analysis-App.git
   cd Advanced-scRNA-seq-Analysis-App
   ```

2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the App**:
   ```bash
   streamlit run app.py
   ```

4. **Access the App**:
   - Open the local URL provided in the terminal to access the app in your browser.

## How to Use

1. **Load the Data**:
   - Make sure your processed scRNA-seq data file (`processed_pbmc_data.h5ad`) is available in the project directory.

2. **Adjust Clustering Resolution**:
   - Use the slider in the sidebar to modify the clustering resolution. The UMAP plot will update based on your selection.

3. **Gene Expression Visualization**:
   - Input a gene name (e.g., CD3D) to visualize its expression levels across clusters.

4. **Analyze Top Marker Genes**:
   - Select a cluster from the sidebar to perform differential expression analysis and display a heatmap of the top marker genes for that cluster.

## Dependencies

- [Streamlit](https://streamlit.io/)
- [Scanpy](https://scanpy.readthedocs.io/)
- [Matplotlib](https://matplotlib.org/)
- [Seaborn](https://seaborn.pydata.org/)
- [Pandas](https://pandas.pydata.org/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request if you have any improvements or suggestions.

---
