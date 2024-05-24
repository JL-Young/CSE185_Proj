# CSE185_Proj: t-sKNEE

This is a package that implements tsne analysis and plotting as a CSE185 project. The main tsne calculation function `t-sKNEE` takes in an anndata object and outputs a n_obs x 2 matrix that contains the x and y coordinates for plotting for each sample. The tsne plotting function `tsKNEE_plot` takes an anndata obects and and generate the tsne plot. 

t-sKNEE implements a subset of functions of `scanpy.tl.tsne` and `sc.pl.tsne`. 

## Installation

The user needs to have the following packages installed: `scanpy`, `matplotlib.pyplot`, `numpy`. This line of code can be run in command line for for installing these packages.

```
pip install scanpy matplolib numpy
```



## Data preprocessing

The input adata object needs to be **Leiden clustered** in both `t-sKNEE` and `tsKNEE_plot`. The anndata object needs to have a column named `'leiden'` in `anndata.obs` dataframe storing the cluster information for the sample. 

The desired quality control should be done previous to using t-sKNEE. 

## Testing 

## Benchmarking