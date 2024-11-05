import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans
import scipy.stats as stats
import umap
import matplotlib.pyplot as plt
import seaborn as sns

# import SLE40 and Select Markers
dataSLE40 = pd.read_csv("/Users/dhyscuduke/Desktop/HYD/Cell_Image_Project/Final_Cell_Tables/SLE40_table.csv")
dataSLE40_markers = ['CD4', 'CD3', 'CD8', 'MPO',  'CD68', 'CD20', 'CD163',
                          'DC-SIGN', 'CD14', 'CD16', 'CD11c', 'CD11b', 'CD56',  'CD45',  'HLADR',
                          'CD38', 'CD31', 'Podoplanin', 'SMA',  'Vimentin', 'PAN-KERATIN', 'cluster_label', 'cell_type_after_qc2']
dataSLE40_clustering = dataSLE40[dataSLE40_markers]

# generate the cluster table by 'cell_type' and 'cluster_label'
cluster_means = dataSLE40_clustering.groupby(['cluster_label', 'cell_type_after_qc2']).mean().reset_index()
cluster_sizes = dataSLE40_clustering.groupby(['cluster_label', 'cell_type_after_qc2']).size().reset_index(name='cell_count')
cluster_summary_dataSLE40 = pd.merge(cluster_means, cluster_sizes, on=['cluster_label', 'cell_type_after_qc2'])

# output cluster table
cluster_summary_dataSLE40.to_csv("/Users/dhyscuduke/Desktop/clustering_dataSLE40.csv", index=False)

# input cluster table
cluster_dataSLE40 = pd.read_excel("/Users/dhyscuduke/Desktop/clustering_dataSLE40.xlsx")

# Combine new cell types to original cell table
dataSLE40_merged = pd.merge(dataSLE40, cluster_dataSLE40[['cluster_label', 'cell_type_after_qc2', 'final_cell_type_2']],
                     on=['cluster_label', 'cell_type_after_qc2'], how='left')
dataSLE40_merged.columns

# final cell table
dataSLE40_merged.drop(columns=['imm_or_non', 'cluster_label_imm_or_non',
       'cluster_label', 'cell_type', 'cell_type_after_qc', 'batch',
       'cell_type_after_qc2'], inplace=True)
dataSLE40_merged['final_cell_type_2'].unique()

########################### QC - CD31
data_Endothelial = dataSLE40_merged[dataSLE40_merged['final_cell_type_2'] == 'Endothelial cells']
data_Endothelial.shape
# Histogram
plt.hist(data_Endothelial['CD31'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

plt.hist(data_Endothelial['PAN-KERATIN'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

plt.hist(data_Endothelial['Vimentin'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
np.mean(data_Endothelial['PAN-KERATIN'])

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'Endothelial cells') & (dataSLE40_merged['CD31'] < 1.25) & (dataSLE40_merged['Vimentin'] > 0.5), 
    'final_cell_type_2'
] = 'Germinal Centers'

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'Endothelial cells') & (dataSLE40_merged['CD31'] < 1.25) & (dataSLE40_merged['PAN-KERATIN'] > 0), 
    'final_cell_type_2'
] = 'Germinal Centers'


########################### QC - CD163
data_163 = dataSLE40_merged[dataSLE40_merged['final_cell_type_2'] == 'CD163+']
data_163.shape
# Histogram
plt.hist(data_163['CD163'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

plt.hist(data_163['CD45'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'CD163+') & (dataSLE40_merged['CD163'] < 3), 
    'final_cell_type_2'
] = 'Unidentified'

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'CD163+') & (dataSLE40_merged['CD163'] >= 3) & (dataSLE40_merged['CD45'] >= 2.5), 
    'final_cell_type_2'
] = 'Leukocytes'


############################# QC - CD4
data_CD4 = dataSLE40_merged[dataSLE40_merged['final_cell_type_2'] == 'CD4T']
data_CD4 = data_CD4[data_CD4['fov'] == '[FOV3-1] UCD151_top_sample_LN17_16-S955_2023-10-16_2-nrm']
data_CD4.shape

plt.hist(data_CD4['CD4'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'CD4T') & (dataSLE40_merged['fov'] == '[FOV3-1] UCD151_top_sample_LN17_16-S955_2023-10-16_2-nrm') & (dataSLE40_merged['CD4'] < 3), 
    'final_cell_type_2'
] = 'Leukocytes'


############################# QC - M2
dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'M2-like macrophages') & (dataSLE40_merged['fov'] == '[FOV30-1] UCD154_bottom_sample_N-5_2023-10-25_3-nrm'), 
    'final_cell_type_2'
] = 'Unidentified'

data_M2 = dataSLE40_merged[dataSLE40_merged['final_cell_type_2'] == 'M2-like macrophages']
data_M2.shape

plt.hist(data_M2['CD68'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'M2-like macrophages') & (dataSLE40_merged['CD68'] < 3), 
    'final_cell_type_2'
] = 'CD163+'

############################## QC - Germinal, PAN, Uniden
data_Germinal = dataSLE40_merged[dataSLE40_merged['final_cell_type_2'] == 'Germinal Centers']
data_Germinal.shape

data_PAN = dataSLE40_merged[dataSLE40_merged['final_cell_type_2'] == 'Distal tubules']
data_PAN.shape

data_Un = dataSLE40_merged[dataSLE40_merged['final_cell_type_2'] == 'Unidentified']
data_Un.shape

plt.hist(data_Germinal['Vimentin'], bins=100, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

plt.hist(data_PAN['PAN-KERATIN'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'Germinal Centers') & (dataSLE40_merged['Vimentin'] < -0.5), 
    'final_cell_type_2'
] = 'Unidentified'

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'Distal tubules') & (dataSLE40_merged['PAN-KERATIN'] < -0.5), 
    'final_cell_type_2'
] = 'Unidentified'

plt.hist(data_Un['Vimentin'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

plt.hist(data_Un['PAN-KERATIN'], bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram Example")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'Unidentified') & (dataSLE40_merged['Vimentin'] > -0.5), 
    'final_cell_type_2'
] = 'Germinal Centers'

dataSLE40_merged.loc[
    (dataSLE40_merged['final_cell_type_2'] == 'Unidentified') & (dataSLE40_merged['PAN-KERATIN'] > -0.5), 
    'final_cell_type_2'
] = 'Distal tubules'

dataSLE40_merged['final_cell_type_2'].unique()

# save the final table
dataSLE40_merged.to_csv("/Users/dhyscuduke/Desktop/Table_SLE40.csv", index=False)
dataSLE40_merged.shape


# generate the fov-cell type table
fov_cell = pd.pivot_table(dataSLE40_merged, index='fov', columns='final_cell_type_2', aggfunc='size', fill_value=0)
fov_cell = fov_cell.reset_index()
fov_cell.to_csv("/Users/dhyscuduke/Desktop/fov_cell_dataSLE40.csv", index=False)
# Only mainitain specific cell types
dataSLE40_merged['cell_type_test'] = dataSLE40_merged['final_cell_type_2']
target_strings = ['CD11c+', 'CD16+', 'CD14+', 'CD163+', 'T cells', 'M1 macrophages', 'M2', 'Monocytes']
dataSLE40_merged['cell_type_test'] = dataSLE40_merged['cell_type_test'].apply(lambda x: x if any(s in x for s in target_strings) else 'Unsure')
dataSLE40_merged.to_csv("/Users/dhyscuduke/Desktop/dataSLE40_merged.csv", index=False)

# ump
columns_to_use = ['CD4', 'CD3', 'CD8', 'CD20', 'MPO', 'CD68', 'CD163', 'DC-SIGN',
                      'CD14', 'CD16', 'CD56', 'CD11b', 'CD11c', 'CD45', 'CD45RO', 'HLADR', 'CD38',
                      'Vimentin', 'CD31', 'Podoplanin', 'SMA', 'PAN-KERATIN']

SLE40_umap = dataSLE40_merged[columns_to_use]

# Fit and transform UMAP
# print('Computing umaps embeddings ...')
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean')
umap_embedding = umap_model.fit_transform(SLE40_umap)

# Create a DataFrame for the UMAP result
umap_df = pd.DataFrame(umap_embedding, columns=[f'UMAP_{i + 1}' for i in range(2)])
# print('umaps complete')
umap_df

dataSLE40_merged.reset_index(drop=True)
dataSLE40_merged = dataSLE40_merged.reset_index(drop=True)
umap_df['cell_type'] = dataSLE40_merged['final_cell_type_2']
umap_df['fov'] = dataSLE40_merged['fov']
umap_df = umap_df[umap_df['cell_type'] != 'Unidentified']
umap_df

umap_df_FSGS = umap_df[(umap_df['fov'].str.contains('_FSGS-'))]
umap_df_LN = umap_df[(umap_df['fov'].str.contains('_LN-'))]
umap_df_IGA = umap_df[(umap_df['fov'].str.contains('_IGA-'))]
umap_df_N = umap_df[(umap_df['fov'].str.contains('_N-'))]
len(umap_df['cell_type'].unique())

# Color
custom_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf', '#9e9e9e', '#f7b6d2', '#ff9896', '#c5b0d5', '#c49c94', '#ecc8a7',
    '#c7c7c7', '#d84a35', '#f0e442', '#31a354', '#756bb1', '#a6cee3', '#b2df8a', '#6a3d9a',
    '#ffed6f', '#b15928', '#cab2d6'
]
unique_batches = umap_df['cell_type'].unique()
num_colors = len(unique_batches)
palette = custom_colors[:num_colors]
color_mapping = dict(zip(unique_batches, palette))

# plot-UMAP
plt.figure(figsize=(10, 8))
sns.scatterplot(x='UMAP_1', y='UMAP_2', hue='cell_type', data=umap_df, palette=color_mapping, s=0.5)
plt.title('UMAP')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()


