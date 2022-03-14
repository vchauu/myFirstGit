import scprep as sc

def main():
    data = sc.io.load_mtx(
        "data/Brain vascular transcriptomes/Normalized values.mtx",
        cell_axis='column',
        gene_names=
        "data/Brain vascular transcriptomes/Normalized rows.mtx_rows",
        cell_names=
        "data/Brain vascular transcriptomes/Normalized columns.mtx_cols",
        sparse=None)

    sc.plot.plot_gene_set_expression(data, genes=None, figsize=(8, 8),
                                     xlabel='Gene expression',
                                     title='Vascular transcriptomes', dpi=300
                                     )

if __name__ == "__main__":
    main()
