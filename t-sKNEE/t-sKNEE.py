import matplotlib.pyplot as plt
import scanpy as sc

def tsKNEE(adata, perp = 30):



# can have the coordinates already generated, or just have a adata object (can run tsne in this function)
def mytsne_plot(adata, coor = None, perp = 30, xlabel = "tsne1", ylabel = "tsne 2", title = "", save = None):
    # import the n_obs x 2 (coordinates)
    if "leiden" not in adata.obs: 
        raise Exception("anndata object is not clustered with Leiden or the Leiden cluster values are not stored in a column named leiden in adata.obs.")
    # if coor = NONE: 
    #     coor = tsKNEE(adata, perp)

    plt.scatter(coors[:,0], coors[:,1], c = adata.obs["leiden"])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show

    if save is not None: 
        plt.savefig(save)

