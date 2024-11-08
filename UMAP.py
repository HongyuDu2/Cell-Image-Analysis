import umap
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import torch as th

def generate_umap(df, umap_cols, color, n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean',
                  sample_size=2000000):

    # Extract numerical columns
    df_umap = df[umap_cols]

    # Sample the data if it's too large
    if len(df_umap) > sample_size:
        df_sampled = df_umap.sample(n=sample_size, random_state=42)
        print(f"Data is too large; sampled {sample_size} rows for UMAP.")
    else:
        df_sampled = df_umap
        print(f"Data is small enough; using all {len(df_umap)} rows for UMAP.")

    # Fit and transform UMAP
    print('Computing umaps embeddings ...')
    umap_model = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, metric=metric)
    umap_embedding = umap_model.fit_transform(df_sampled)

    # Create a DataFrame for the UMAP result
    umap_df = pd.DataFrame(umap_embedding, columns=[f'UMAP_{i + 1}' for i in range(n_components)])
    print('umaps complete')

    # Add UMAP columns to the original DataFrame
    df_umap_full = df.copy()
    df_umap_full.loc[df_sampled.index, f'UMAP_1'] = umap_df['UMAP_1']
    df_umap_full.loc[df_sampled.index, f'UMAP_2'] = umap_df['UMAP_2']

    print('umaps dataset complete')

    umap_df['batch'] = df[color]

    # Plotting
    sns.scatterplot(x='UMAP_1', y='UMAP_2', hue='batch', data=umap_df, s=0.5)
    plt.tight_layout()
    plt.title('batch effect')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    #plt.savefig(os.path.join(base_dir, "umap_5_cell_table_with_cell_type_only_identified_sle40_final.png"), bbox_inches='tight')
    plt.show()

    # Return the UMAP DataFrame
    return umap_df




    df_sorted = pd.read_csv("C:/Users/duho/Desktop/Umap_data.csv")

    columns_to_use = ['CD4', 'CD3', 'CD8', 'CD20', 'MPO', 'CD68', 'CD163', 'DC-SIGN',
                      'CD14', 'CD16', 'CD56', 'CD11b', 'CD11c', 'CD45', 'CD45RO', 'HLADR', 'CD38',
                      'Vimentin', 'CD31', 'Podoplanin', 'SMA', 'PAN-KERATIN']

    #colors = ['batch']

    df_sorted_shuffle = df_sorted.sample(frac=1).reset_index(drop=True)

    umap_result_df = generate_umap(df_sorted_shuffle, columns_to_use, color='batch', n_neighbors=15, min_dist=0.1,
                                   n_components=2, metric='euclidean', sample_size=2000000)





    df_umap = df_sorted_shuffle[columns_to_use]
    df_sampled = df_umap
                    

    # Fit and transform UMAP
    print('Computing umaps embeddings ...')
    umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean')
    umap_embedding = umap_model.fit_transform(df_sampled)

    # Create a DataFrame for the UMAP result
    umap_df = pd.DataFrame(umap_embedding, columns=[f'UMAP_{i + 1}' for i in range(2)])
    print('umaps complete')

    umap_df['batch'] = df_sorted_shuffle['batch']

    umap_df['batch'] = df.loc[df_sampled.index, color]

    # Ensure unique colors for each batch
    unique_batches = umap_df['batch'].unique()
    palette = sns.color_palette('hsv', len(unique_batches))
    color_mapping = dict(zip(unique_batches, palette))

    # Plotting
    sns.scatterplot(x='UMAP_1', y='UMAP_2', hue='batch', data=umap_df, palette=color_mapping, s=0.5)
    plt.tight_layout()
    plt.title('batch effect')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    #plt.savefig(os.path.join(base_dir, "umap_5_cell_table_with_cell_type_only_identified_sle40_final.png"), bbox_inches='tight')
    plt.show()





    ######### Randomly selected 200,000 dataset in each batch and plot Umap
    sample_size_per_batch = 199886

    # Use groupby and sample
    sampled_df = umap_df.groupby('batch', group_keys=False).apply(
        lambda x: x.sample(min(len(x), sample_size_per_batch)))

    # Check dataset
    print(sampled_df)
    print(sampled_df['batch'].value_counts())

    sampled_df_shuffle = sampled_df.sample(frac=1).reset_index(drop=True)

    # Plotting
    sns.scatterplot(x='UMAP_1', y='UMAP_2', hue='batch', data=sampled_df_shuffle, palette=color_mapping, s=0.5)
    plt.tight_layout()
    plt.title('batch effect')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.show()









