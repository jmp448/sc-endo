#!/bin/python
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import limix

COUNTS_LOC = "/work-zfs/abattle4/lab_data/sc_endo_diff/counts.tsv"
METADATA_LOC = "/work-zfs/abattle4/lab_data/sc_endo_diff/cell_metadata_cols.tsv"
GENENAMES_LOC = "/home-2/jpopp4@jhu.edu/work/josh/sc-endoderm-differentiation/data/genenames.txt"

def load_data():
    counts = np.loadtxt(COUNTS_LOC, delimiter="\t", skiprows=1, usecols=range(1, 36045))
    metadata = pd.read_csv(METADATA_LOC, sep="\t", header=0)
    gene_names = pd.read_csv(GENENAMES_LOC, sep=' ', header=0)
    adata = anndata.AnnData(X=counts, obs=metadata)
    adata.var_names = gene_names.keys()
    return adata

def dim_reduce(adata):
    sc.pp.scale(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=500)
    sc.tl.pca(adata, svd_solver='arpack')
    pve1 = adata.uns['pca']['variance_ratio'][0]
    pve2 = adata.uns['pca']['variance_ratio'][1]
    sc.pl.pca(adata, color='day', save=".png")
    return adata

def main():
    results_file="./data/sc.endo.h5ad"
    adata = load_data()
    adata = dim_reduce(adata)
    adata.write(results_file)


if __name__=="__main__":
    main()
