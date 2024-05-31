# CSE185_Proj: t-sKNEE

This package implements tsne analysis and plotting. The main tsne calculation function `tsKNEE` takes in an anndata object and outputs a n_obs x 2 matrix that contains the x and y coordinates for each sample. The tsne plotting function `tsKNEE_plot` takes an anndata obect and generates a tsne plot. 

t-sKNEE implements a subset of functions of `scanpy.tl.tsne` and `sc.pl.tsne`. For more information about these functions visit [scanpy](https://scanpy.readthedocs.io/en/stable/api/tools.html) page.

## Installation
First, clone the repository using the following command `git clone https://github.com/JL-Young/CSE185_Proj.git`
The user needs to have the following libraries installed: `matplotlib.pyplot`, `numpy`, `scanpy` as well as `leidenalg`. The following lines of code can be run in command line for for installing these packages.

```
pip install matplotlib as plt
pip install numpy as np
pip install scanpy as sc
pip install leidenalg
```
Once the required libraries are installed, you can install `t_sKNEE` with the following command.
```python setup.py install```

If you do not have admin access, the packages can be installed using the following commands.
```
pip install --user matplolib numpy
pip install --user numpy as np
pip install --user scanpy as sc
pip install --user leidenalg
python setup.py install --user
```

## Basic usage

There are two functions you can run wihtin `tsKNEE`- `tsKNEE` and `tsKNEE_plot`. The only required parameters for both functions is an anndata object. See Data Processing section for more information on the anndata object.

The basic usage of `tsKNEE` is: 
```
tsKNEE(adata, T=1000, perp = 30)
```
- `T=1000`: the number of iterations tsKNEE goes through to plot samples with according to a similarity matrix
- `perp=30`: the perplexity 

The basic usage of `tsKNEE_plot` is: 
```
tsKNEE_plot(adata, perp = 30, xlabel = "tsne1", ylabel = "tsne 2", title = "", save = None)
```
- `perp = 30`: the perplexity to run `tsKNEE` if necessary
- `xlabel = "tsne1"`: the x-axis label for the graph
- `ylabel = "tsne 2"`: the y-axis label for the graph
- `title = ""`: the graph title
- `save = None`: to save graph in a file

## Data preprocessing

The input adata object needs to be **Leiden clustered** in both `tsKNEE` and `tsKNEE_plot`. The anndata object needs to have a column named `leiden` in `anndata.obs` dataframe storing the cluster information for the sample. For `tsKNEE_plot`, the input anndata object must also have `X_tsne` within `anndata.obsm` which is the output of the `tsKNEE` function. 

The desired quality control should be done previous to using t-sKNEE. 

## Testing 

To run `tsKNEE` and `tsKNEE_plot` on a small test example run all cells up to cell 8 in `test_tsKNEE.ipynb`. Cell 8 gives coordinates of each sample in the test data for tsKNEE plotting. 

## Benchmarking

## Contributors

This repository was generated by Jane Li, Jeyasri Venkatasurbamani, and James Young with inspiration from the [Medium](https://towardsdatascience.com/understanding-t-sne-by-implementing-2baf3a987ab3).
