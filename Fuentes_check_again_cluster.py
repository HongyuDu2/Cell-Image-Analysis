import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans
import scipy.stats as stats
import umap
import matplotlib.pyplot as plt
import seaborn as sns

# import Fuentes and Select Markers
dataFuentes = pd.read_csv("/Users/dhyscuduke/Desktop/HYD/Cell_Image_Project/Final_Cell_Tables/Fuentes_table.csv")
dataFuentes_markers = ['CD4', 'CD3', 'CD8', 'MPO',  'CD68', 'CD20', 'CD163',
                          'DC-SIGN', 'CD14', 'CD16', 'CD11c', 'CD11b', 'CD56',  'CD45',  'HLADR',
                          'CD38', 'CD31', 'Podoplanin', 'SMA',  'Vimentin', 'PAN-KERATIN', 'cluster_label', 'final_cell_type']
dataFuentes_clustering = dataFuentes[dataFuentes_markers]

# generate the cluster table by 'cell_type' and 'cluster_label'
cluster_means = dataFuentes_clustering.groupby(['cluster_label', 'final_cell_type']).mean().reset_index()
cluster_sizes = dataFuentes_clustering.groupby(['cluster_label', 'final_cell_type']).size().reset_index(name='cell_count')
cluster_summary_dataFuentes = pd.merge(cluster_means, cluster_sizes, on=['cluster_label', 'final_cell_type'])

# output cluster table
cluster_summary_dataFuentes.to_csv("/Users/dhyscuduke/Desktop/clustering_dataFuentes.csv", index=False)

# input cluster table
cluster_dataFuentes = pd.read_excel("/Users/dhyscuduke/Desktop/HYD/Cell_Image_Project/Check_QC_Again/Final Check_Fuentes/clustering_dataFuentes.xlsx")

# combine new cell types to original cell table
dataFuentes_merged = pd.merge(dataFuentes, cluster_dataFuentes[['cluster_label', 'final_cell_type', 'final_cell_type_2']],
                     on=['cluster_label', 'final_cell_type'], how='left')

# final cell table
dataFuentes_merged.drop(columns=['imm_or_non', 'cluster_label_imm_or_non',
       'cluster_label', 'cell_type', 'cell_type_after_qc', 'batch',
       'cell_type_after_qc2', 'final_cell_type'], inplace=True)
dataFuentes_merged['final_cell_type_2'].unique()

data_distal = dataFuentes_merged[dataFuentes_merged['final_cell_type_2'] == 'Distal tubules']
data_distal.shape

# Histogram
plt.hist(data_distal['PAN-KERATIN'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
np.mean(data_distal['PAN-KERATIN'])

data_uniden = dataFuentes_merged[dataFuentes_merged['final_cell_type_2'] == 'Unidentified']
data_uniden.shape

# Histogram
plt.hist(data_uniden['PAN-KERATIN'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
np.mean(data_distal['PAN-KERATIN'])

# unidentified to distal tubules
dataFuentes_merged.loc[
    (dataFuentes_merged['final_cell_type_2'] == 'Unidentified') & (dataFuentes_merged['PAN-KERATIN'] > 0), 
    'final_cell_type_2'
] = 'Distal tubules'

dataFuentes_merged.to_csv("/Users/dhyscuduke/Desktop/Table_Fuentes.csv", index=False)
dataFuentes_merged.shape

# generate the fov-cell type table
fov_cell = pd.pivot_table(dataFuentes_merged, index='fov', columns='final_cell_type_1', aggfunc='size', fill_value=0)
fov_cell = fov_cell.reset_index()
fov_cell.to_csv("/Users/dhyscuduke/Desktop/fov_cell_dataFuentes.csv", index=False)

# only mainitain specific cell types
dataFuentes_merged['cell_type_test'] = dataFuentes_merged['final_cell_type_1']
target_strings = ['M1 macrophages-', 'M2-like macrophages-', 'Monocytes-', 'T cells-', 'Unidentified-']
dataFuentes_merged['cell_type_test'] = dataFuentes_merged['cell_type_test'].apply(lambda x: x if any(s in x for s in target_strings) else 'Unsure')
dataFuentes_merged.to_csv("/Users/dhyscuduke/Desktop/dataFuentes_merged.csv", index=False)

dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'CD163+-13'), 'final_cell_type_1'] = 'CD163+'
dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'CD163+-43'), 'final_cell_type_1'] = 'CD163+'
dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'CD163+-89'), 'final_cell_type_1'] = 'Unidentified'
dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'CD163+-139'), 'final_cell_type_1'] = 'Unidentified'

dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'DC-27'), 'final_cell_type_1'] = 'Unidentified'
dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'DC-63'), 'final_cell_type_1'] = 'Unidentified'
dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'DC-116'), 'final_cell_type_1'] = 'DC'
dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'DC-136'), 'final_cell_type_1'] = 'Unidentified'
dataFuentes_merged.loc[(dataFuentes_merged['final_cell_type_1'] == 'DC-164'), 'final_cell_type_1'] = 'DC'

dataFuentes_merged = pd.merge(dataFuentes_merged, cluster_dataFuentes[['cluster_label', 'final_cell_type', 'final_cell_type_2']],
                     on=['cluster_label', 'final_cell_type'], how='left')
dataFuentes_merged['final_cell_type_2'].unique()

Fuentes_FSGS = dataFuentes_merged[(dataFuentes_merged['fov'].str.contains('_FSGS-'))]
Fuentes_LN = dataFuentes_merged[(dataFuentes_merged['fov'].str.contains('_LN-'))]
Fuentes_IGA = dataFuentes_merged[(dataFuentes_merged['fov'].str.contains('_IGA-'))]
Fuentes_N = dataFuentes_merged[(dataFuentes_merged['fov'].str.contains('_N-'))]
Fuentes_FSGS = Fuentes_FSGS[Fuentes_FSGS['cell_type'] != 'Unidentified']
Fuentes_LN = Fuentes_LN[Fuentes_LN['cell_type'] != 'Unidentified']
Fuentes_IGA = Fuentes_IGA[Fuentes_IGA['cell_type'] != 'Unidentified']
Fuentes_N = Fuentes_N[Fuentes_N['cell_type'] != 'Unidentified']

# create umap
dataFuentes_merged = dataFuentes_merged[dataFuentes_merged['final_cell_type_2'] != 'Unidentified']
columns_to_use = ['CD4', 'CD3', 'CD8', 'CD20', 'MPO', 'CD68', 'CD163', 'DC-SIGN',
                      'CD14', 'CD16', 'CD56', 'CD11b', 'CD11c', 'CD45', 'CD45RO', 'HLADR', 'CD38',
                      'Vimentin', 'CD31', 'Podoplanin', 'SMA', 'PAN-KERATIN']

# Fit and transform UMAP
print('Computing umaps embeddings ...')
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean')
umap_embedding = umap_model.fit_transform(Fuentes_umap)

# Create a DataFrame for the UMAP result
umap_df = pd.DataFrame(umap_embedding, columns=[f'UMAP_{i + 1}' for i in range(2)])
print('umaps complete')

dataFuentes_merged.reset_index(drop=True)
dataFuentes_merged = dataFuentes_merged.reset_index(drop=True)
umap_df['cell_type'] = dataFuentes_merged['final_cell_type_2']
umap_df['fov'] = dataFuentes_merged['fov']

umap_df_5245 = umap_df[(umap_df['fov'].str.contains('S5245'))]
umap_df_3513 = umap_df[(umap_df['fov'].str.contains('S3513'))]
umap_df_5641 = umap_df[(umap_df['fov'].str.contains('S5641'))]
umap_df_N = umap_df[(umap_df['fov'].str.contains('_N-'))]
len(umap_df_5245['cell_type'].unique())

# Color
custom_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf', '#9e9e9e', '#f7b6d2', '#ff9896', '#c5b0d5', '#c49c94', '#ecc8a7',
    '#c7c7c7', '#d84a35', '#f0e442', '#31a354', '#756bb1', '#a6cee3', '#b2df8a'
]
unique_batches = umap_df['cell_type'].unique() 
num_colors = len(unique_batches)
palette = custom_colors[:num_colors]
color_mapping = dict(zip(unique_batches, palette))

# plot
plt.figure(figsize=(10, 8))
sns.scatterplot(x='UMAP_1', y='UMAP_2', hue='cell_type', data=umap_df_5641, palette=color_mapping, s=2.5)
plt.title('UMAP')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left') 
plt.show()



# Bar Plot
dataFuentes_merged.columns
Samples = ['S5245', 'S3513', 'S5641']
counts = {Sample: dataFuentes_merged['fov'].str.contains(Sample).sum() for Sample in Samples}
counts_df = pd.DataFrame(list(counts.items()), columns=['Sample', 'Count'])
plt.figure(figsize=(8, 6))
sns.barplot(x='Sample', y='Count', data=counts_df)