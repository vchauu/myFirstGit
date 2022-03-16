import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sn


def main():
    # Read the input
    sham_1 = pd.read_csv("data/I_R_txt/Cortex-Sham-Rep1_WT-MCAO.txt", sep='\t',
                         skiprows=2, index_col=0
                         )
    sham_2 = pd.read_csv("data/I_R_txt/Cortex-Sham-Rep2_WT-MCAO.txt", sep='\t',
                         skiprows=2, index_col=0
                         )
    sham_3 = pd.read_csv("data/I_R_txt/Cortex-Sham-Rep3_WT-MCAO.txt", sep='\t',
                         skiprows=2, index_col=0
                         )
    sham_4 = pd.read_csv("data/I_R_txt/Cortex-Sham-Rep4_WT-MCAO.txt", sep='\t',
                         skiprows=2, index_col=0
                         )
    ir2_1 = pd.read_csv("data/I_R_txt/Cortex-I_R-2h-Rep1_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    ir2_2 = pd.read_csv("data/I_R_txt/Cortex-I_R-2h-Rep2_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    ir2_3 = pd.read_csv("data/I_R_txt/Cortex-I_R-2h-Rep3_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    ir2_4 = pd.read_csv("data/I_R_txt/Cortex-I_R-2h-Rep4_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    ir8_1 = pd.read_csv("data/I_R_txt/Cortex-I_R-8h-Rep1_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    ir8_2 = pd.read_csv("data/I_R_txt/Cortex-I_R-8h-Rep2_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    ir8_3 = pd.read_csv("data/I_R_txt/Cortex-I_R-8h-Rep3_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    ir8_4 = pd.read_csv("data/I_R_txt/Cortex-I_R-8h-Rep4_WT-MCAO.txt", sep='\t',
                        skiprows=2, index_col=0
                        )
    sham = pd.concat([sham_1, sham_2, sham_3, sham_4], axis=1)
    ir2 = pd.concat([ir2_1, ir2_2, ir2_3, ir2_4], axis=1)
    ir8 = pd.concat([ir8_1, ir8_2, ir8_3, ir8_4], axis=1)
    hm_data = pd.concat([sham, ir2, ir8], axis=1)

    # Read the platform file to convert gene ID to gene names
    platform = pd.read_csv("data/I_R_txt/Platform_Stroke.txt", sep='\t',
                           skiprows=28, index_col=0
                           )
    genes = platform['Symbol']

    hm_data = pd.concat([hm_data, genes], axis=1)
    hm_data.set_index('Symbol', inplace=True)
    groups = ['Sham', 'Sham', 'Sham', 'Sham', 'IR-2h', 'IR-2h', 'IR-2h',
              'IR-2h', 'IR-8h', 'IR-8h', 'IR-8h', 'IR-8h']
    """
    hm_data.columns = groups
    sn.clustermap(hm_data, z_score=None, standard_scale=1, cmap="coolwarm",
                  col_cluster=False, row_cluster=True
                  )
    plt.show()
"""
    # Do PCA
    transposed_data = hm_data.T
    x = transposed_data.values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(x)
    principal_df = pd.DataFrame(data=principal_components,
                                columns=['principal component 1',
                                         'principal component 2']
                                )
    principal_df['Groups'] = groups

    # Visualize 2D projection
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 Component PCA', fontsize=20)

    targets = ['Sham', 'IR-2h', 'IR-8h']
    colors = ['r', 'g', 'b']
    for target, color in zip(targets, colors):
        indices_to_keep = principal_df['Groups'] == target
        ax.scatter(principal_df.loc[indices_to_keep, 'principal component 1'],
                   principal_df.loc[indices_to_keep, 'principal component 2'],
                   c=color, s=50
                   )
    ax.legend(targets)
    ax.grid()
    plt.show()

    # Calculate explained variance
    # print(pca.explained_variance_ratio_)
    loadings = pca.components_
    num_pc = pca.n_features_
    pc_list = ["PC" + str(i) for i in list(range(1, num_pc + 1))]
    loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
    loadings_df['variable'] = hm_data.index
    loadings_df = loadings_df.set_index('variable')
    loadings_df.sort_values(by='PC1', key=abs, inplace=True)
    print(loadings_df)

if __name__ == "__main__":
    main()
