##
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
import seaborn as sns
from statsmodels.stats.multitest import multipletests


##
# This function reads the input files. This function takes the name of the
# brain area of interest (striatum or cortex), the names of the experimental
# groups, the names of the replicates and return a multi-indexed DataFrame
# with data organized into their corresponding experimental group and replicate.
# Inputs:
# name: string, the data name
# group: array of strings, names of experimental groups
# replicate: array of strings, names of replicates
# Output: a multi-indexed DataFrame
# Assumptions for inputs:
# name: a string, either "Cortex" or "Striatum" (case sensitive)
# group and replicate are non-empty arrays of strings
# group is ordered alphabetically so that the data are assigned the correct
# group.
def read_input(name, group, replicate):
    iterables = [group, replicate]
    df = pd.DataFrame()
    # Read the files in the dataset
    i = 0  # index to keep track of sample name
    for file in sorted(os.listdir("data/I_R_txt")):
        if file.endswith(".txt") and file.startswith(name):
            df[str(i)] = pd.read_csv("data/I_R_txt/" + file, sep='\t',
                                     skiprows=2, index_col=0
                                     )
            i += 1
    ind = pd.MultiIndex.from_product(iterables, names=["group", "replicate"])
    df.columns = ind
    return df


##
# Do PLS-DA
# This function performs PLS-DA analysis on a given DataFrame and number of
# principal components of interest and returns a DataFrame containing the
# scores of each sample when projected onto the new PCs.
# Inputs:
# df: a multi-indexed DataFrame
# n_pc: int, number of principal components axes
# n_groups: number of experimental groups
# n_samples: number of samples in each experimental group
# Output: a DataFrame containing the scores of the PLS-DA analysis
# Assumptions for inputs:
# df is a multi-indexed DataFrame, with each sample assigned to a specific
# experimental group and replicate
# n_pc is an int in the range (0, min(n_samples, n_features))
def do_plsda(df, n_pc, n_groups, n_samples):
    # Assign each sample to a group
    y = []
    for i in range(0, n_groups):
        y += [i] * n_samples
    plsr = PLSRegression(n_components=n_pc, scale=False)
    plsr.fit(df.values.T, y)
    scores = pd.DataFrame(plsr.x_scores_)
    scores.index = df.columns
    return plsr, scores


##
# This function plots the PCA or PLS-DA clustering from a given DataFrame and
# color
# scheme.
# Inputs:
# name: title of the plot
# scores: DataFrame containing the PLS-DA scores
# targets: an array of strings containing all the experimental group names
# colors: an array of strings containing all the corresponding colors to the
# experimental groups.
# analysis: string, the name of analysis performed
# ax1: int, which PC or LV on the x-axis
# ax2: int, which PC or LV on the y-axis
# Output: a scatter plot, returns None
# Assumptions for inputs:
# name is a string
# scores is a DataFrame generated from PLS-DA analysis
# targets is an array of strings containing the names of the experimental
# groups. Note that these names should be the same as the multi-indices in
# the scores DataFrame
# colors is an array of string indicating colors
# analysis is either 'plsda' or 'pca'
def plot_scores(name, scores, targets, colors, analysis, ax1, ax2):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    if analysis == 'pca':
        ax.set_xlabel('Principal Component ' + str(ax1), fontsize=15)
        ax.set_ylabel('Principal Component ' + str(ax2), fontsize=15)
        ax.set_title('PCA ' + name, fontsize=20)
    if analysis == 'plsda':
        ax.set_xlabel('Scores on LV ' + str(ax1), fontsize=15)
        ax.set_ylabel('Scores on LV ' + str(ax2), fontsize=15)
        ax.set_title('PLS-DA ' + name, fontsize=20)

    # Plot the data points on a scatter plot
    for target, color in zip(targets, colors):
        ax.scatter(scores.loc[target, ax1 - 1], scores.loc[target, ax2 - 1],
                   c=color, s=50
                   )
    ax.legend(targets)
    plt.show()
    return


##
# This function performs principal component analysis (PCA) on the given
# DataFrame and number of principal components (PCs) . This function returns
# the PCA result and a DataFrame containing the loadings of the PCs.
# Inputs:
# df: a DataFrame of the input
# n_pc: int, number of PCs
# Output:
# PCA results
# A DataFrame containing the PC loadings
def do_pca(df, n_pc):
    transposed_data = df.T
    x = transposed_data.values
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=n_pc)
    principal_components = pca.fit_transform(x)
    principal_df = pd.DataFrame(data=principal_components)
    principal_df.index = df.columns
    return pca, principal_df


##
# This function plots a scree plot and a cumulative explained variance from a
# given PCA result object.
# Inputs:
# pca: an object of type PCA from a PCA calculation
# name: string, title for the plot
# Output: a scree plot, an array containing explained variance ratio
# Assumptions for the inputs:
# pca is the result returned from calling PCA from the scikit package
def plot_pca_var(pca, name):
    pc_values = np.arange(pca.n_components_) + 1
    y = pca.explained_variance_ratio_

    # Plot the scree plot
    plt.figure(figsize=(8, 8))
    plt.plot(pc_values, y, 'o-', linewidth=2, color='b')
    plt.title('Scree Plot ' + name, fontsize=18)
    plt.xlabel('Principal Component', fontsize=18)
    plt.ylabel('Variance Explained', fontsize=18)
    plt.show()
    return y


##
# This function returns a DataFrame containing the loadings of each feature
# on each principal component.
# Inputs:
# pca: PCA result return from running PCA
# df: DataFrame containing original data
# Output: a DataFrame containing the loadings of genes in the original data
def calc_loadings(pca, df):
    loadings = pca.components_
    num_pc = pca.n_features_
    pc_list = ["PC" + str(i) for i in list(range(1, num_pc + 1))]
    loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
    loadings_df.set_index(df.index, inplace=True)
    return loadings_df


##
# This function returns the top genes with their loadings for the principal
# component of choice, ordered from most up-regulated to most down-regulated.
# Inputs:
# loadings_df: DataFrame containing loading factors
# n_genes: int, number of most differentially expressed genes to extract
# filter_by: string, the principal component to sort the genes by
# Output: DataFrame containing the loadings of the top genes
def calc_top_genes_loadings(loadings_df, n_genes, filter_by):
    loadings_df.sort_values(by=filter_by, key=abs, inplace=True, ascending=False
                            )
    top_genes_df = loadings_df.head(n_genes).copy()
    top_genes_df.sort_values(by=filter_by, inplace=True, ascending=False)
    return top_genes_df


##
# This function converts gene indices following the ILMN ID to gene names.
# This function takes in a DataFrame with ILMN indices and return a DataFrame
# with gene names as indices.
# Inputs:
# df: DataFrame to convert indices
# Output: the original DataFrame with an extra column with the gene names
# Assumption:
# platform has ILMN ID as indices
def convert_to_gene_names(df):
    # Platform to convert gene IDs to gene names
    platform = pd.read_csv("data/I_R_txt/Platform_Stroke.txt", sep='\t',
                           skiprows=28, index_col=0
                           )
    symbols = platform.loc[df.index.to_numpy(), 'Symbol'].to_numpy()
    df['gene_name'] = symbols
    return df


##
# This function plots a heatmap from a DataFrame containing gene expression
# from different samples.
# Inputs:
# df: DataFrame containing original data
# name: string, name for the graph title
# Output: a heat map for gene expression
# Assumptions for input:
# df is a DataFrame containing gene expression levels for all samples
def plot_hm(df, title):
    kws = dict(cbar_kws=dict(ticks=[-1.5, 0, 1.5], orientation='horizontal'),
               figsize=(8, 8)
               )
    sns.set(font_scale=0.8)
    g = sns.clustermap(df, z_score=0, cmap="coolwarm", col_cluster=False,
                       row_cluster=True, yticklabels=True, **kws
                       )

    g.fig.suptitle("Heatmap of " + title, size='large', weight='bold')
    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position(
        [x0, 0.9, g.ax_row_dendrogram.get_position().width, 0.02]
        )
    g.ax_cbar.set_title('z-score')
    g.ax_cbar.tick_params(axis='x', length=10)
    for spine in g.ax_cbar.spines:
        g.ax_cbar.spines[spine].set_color('crimson')
        g.ax_cbar.spines[spine].set_linewidth(2)
    return g


##
# This function performs multiple hypothesis testing to identify genes with
# differential expression at two time points. This function takes a
# multi-indexed DataFrame, whose columns contain information about the
# experimental group of each sample.
# Inputs:
# df: DataFrame containing gene expression of all samples
# g1: str, name of 1st group
# g2: str, name of 2nd group
# g3: str, name of 3rd group
# g4: str, name of 4th group
# Output: a new DataFrame containing the gene IDs and the p-values
def do_mht(df, g1, g2, g3, g4):
    num_genes = df.shape[0]
    results = []
    for i in range(0, num_genes):
        p = stats.f_oneway(df.loc[df.index[i], g1], df.loc[df.index[i], g2],
                           df.loc[df.index[i], g3], df.loc[df.index[i], g4]
                           )[1]
        results.append(p)
    p_corrected = multipletests(results, alpha=0.05, method='fdr_bh')[1]
    p_df = pd.DataFrame(data=p_corrected, index=df.index, columns=['p'])
    return p_df
