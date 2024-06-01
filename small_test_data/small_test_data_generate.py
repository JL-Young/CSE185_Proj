import os
import scanpy as sc, anndata as ad 
import random
#sc.logging.print_versions()
DATADIR = os.getcwd()
print(os.getcwd())
ds = "GSM5114474_M3_E7"
adata = sc.read_10x_mtx(DATADIR, prefix=ds+"_", cache=True)

sc.pp.filter_cells(adata, min_counts = 1000, inplace = True)
sc.pp.filter_cells(adata, min_genes = 200, inplace = True)
sc.pp.filter_genes(adata, min_counts = 15, inplace = True)
sc.pp.filter_genes(adata, min_cells = 5, inplace = True)

print(adata) 
#n_obs × n_vars = 2087 × 13384
#keep 50 samples and 1000 genes

adata.obs["rand_obs"] = [random.random() for i in range(adata.n_obs)]
adata.var["rand_vars"] = [random.random() for i in range(adata.n_vars)]

adata = adata[adata.obs.rand_obs < 0.028, adata.var.rand_vars < 0.07]
print(adata)
adata.write_h5ad("small_test_data2.h5ad")