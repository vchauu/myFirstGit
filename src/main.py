import scprep as sc
import pandas as pd
import matplotlib.pyplot as plt

def main():
    data = sc.io.load_mtx(
        "data/E-GEOD-98816.aggregated_filtered_normalised_counts.mtx",
        cell_axis='column',
        gene_names=
        "data/E-GEOD-98816.aggregated_filtered_normalised_counts.mtx_rows",
        cell_names=
        "data/E-GEOD-98816.aggregated_filtered_normalised_counts.mtx_cols",
        sparse=None)
    sc.plot.plot_gene_set_expression(data, genes=None, figsize=(8,8),
                                     xlabel='Gene expression',
                                     title='Vascular transcriptomes',
                                     dpi=300)


if __name__ == "__main__":
    main()
