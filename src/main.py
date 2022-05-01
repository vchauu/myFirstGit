##
from src.process_data import *

##
# Read inputs
group = ['ir24', 'ir2', 'ir8', 'sham']
replicate = ['one', 'two', 'three', 'four']
cortex = read_input("Cortex", group, replicate)
striatum = read_input("Striatum", group, replicate)

##
# Set the groups and their colors
colors = ['r', 'g', 'b', 'darkgoldenrod']

##
# Do PCA
pca_cor_full = do_pca(cortex, 16)
pca_stri_full = do_pca(striatum, 16)
pca_cor = pca_cor_full[0]
pca_stri = pca_stri_full[0]

##
# Plot PCA plots
pc_df_cor = pca_cor_full[1]
pc_df_stri = pca_stri_full[1]
plot_scores("Striatum", pc_df_stri, group, colors, 'pca', 1, 2)
plot_scores("Cortex", pc_df_cor, group, colors, 'pca', 1, 2)

##
# Plot variance explained by PCs
pca_explained_var_cor = plot_pca_var(pca_cor, "Cortex")
pca_explained_var_stri = plot_pca_var(pca_stri, "Striatum")

##
# pca_cor = do_pca(cortex, 5)[0]
pca_stri = do_pca(striatum, 4)[0]
##
# Rank striatal genes by their loadings on PC2
loadings_stri = calc_loadings(pca_stri, striatum)
loadings_stri.sort_values(by='PC2', key=abs, inplace=True, ascending=False
                          )
convert_to_gene_names(loadings_stri)
loadings_stri['Ranking'] = range(1, loadings_stri.shape[0] + 1)
##
# Perform one-way ANOVA to identify genes with highly differential expressions
num_genes = striatum.shape[0]
results = []
for i in range(0, num_genes):
    p = stats.f_oneway(striatum.loc[striatum.index[i], 'sham'],
                       striatum.loc[striatum.index[i], 'ir2'],
                       striatum.loc[striatum.index[i], 'ir8'],
                       striatum.loc[striatum.index[i], 'ir24']
                       )[1]
    results.append(p)
p_corrected = multipletests(results, alpha=0.05, method='fdr_bh')[1]
##
p_df = pd.DataFrame(data=p_corrected, index=striatum.index, columns=['p'])
convert_to_gene_names(p_df)
sig_genes = p_df[p_df['p'] < 0.05]
##
sig_genes.to_csv("results/mht_striatum.csv")

##
# Make a heatmap of the highly differential genes
striatum = striatum.reindex(
    columns=striatum.columns.reindex(['sham', 'ir2', 'ir8', 'ir24'], level=0)[0]
    )
convert_to_gene_names(striatum)
hm_genes = striatum.loc[sig_genes.index, :]
hm_genes.set_index('gene_name', inplace=True)
hm_genes_plot = plot_hm(hm_genes,
                "genes with highly differential expression in the striatum")
hm_genes_plot.savefig("results/deg_striatum.png")
##
# Read the protein atlas file to extract genes encoding proteins that
# circulate in the blood stream
all_genes = pd.read_csv("data/proteinatlas.tsv", sep='\t')
genes_blood = all_genes.filter(
    ['Gene', 'CCD Protein', 'Blood concentration - Conc. blood MS [pg/L]',
     'Blood concentration - Conc. blood IM [pg/L]']
    )
genes_blood.dropna(thresh=2, inplace=True)

##
# Find the genes encoding blood-circulating proteins
blood_arr = genes_blood['Gene'].to_numpy()
stri_diff = sig_genes['gene_name'].str.upper().to_numpy()
stri_blood = np.intersect1d(blood_arr, stri_diff)

##
# Get ILMN IDs of protein-encoding proteins
ilmn_ids = sig_genes.loc[sig_genes['gene_name'].str.upper().isin(
    stri_blood)].index
# Make a heatmap of genes encoding blood-circulating proteins
hm_blood = striatum.loc[ilmn_ids, :]
hm_blood.set_index('gene_name', inplace=True)
hm_blood_plot = plot_hm(hm_blood,
                        "protein-encoding genes with highly differential " +
                        "expression in the striatum"
                        )
hm_blood_plot.savefig("results/blood_proteins_striatum.png")

##
# Get loadings rankings of the identified genes
rankings = loadings_stri.loc[ilmn_ids, :].sort_values(by='Ranking')
rankings.loc[:, ('gene_name', 'Ranking')].to_csv("results/deg_rankings.csv")
##
# Obtain regression calculations
plsr_cor = do_plsda(cortex, 9, len(group), len(replicate))[0]
plsr_stri = do_plsda(striatum, 9, len(group), len(replicate))[0]

scores_cor = pd.DataFrame(plsr_cor.x_scores_)
scores_cor.index = cortex.columns

scores_stri = pd.DataFrame(plsr_stri.x_scores_)
scores_stri.index = striatum.columns

##
# Plot PLS-DA graphs
plot_scores("Striatum", scores_stri, group, colors, 'plsda', 1, 2)
plot_scores("Cortex", scores_cor, group, colors, 'plsda', 1, 2)
