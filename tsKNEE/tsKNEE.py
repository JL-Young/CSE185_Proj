import matplotlib.pyplot as plt, matplotlib.colors as mcolors
import random
import scanpy as sc
import numpy as np 
import scanpy as sc, anndata as ad

#find euclidean distance
def pairwise_distances(X):
    return np.sum((X[None, :] - X[:, None])**2, 2)

#normalize gaussian distribution
def g_conditional(dists, sigmas):
    e = np.exp(-dists / (2 * np.square(sigmas.reshape((-1,1)))))
    np.fill_diagonal(e, 0.)
    e += 1e-8
    return e / e.sum(axis=1).reshape([-1,1])

#adjust matrix to find the best sigma 
def perp(condi_matr):
    ent = -np.sum(condi_matr * np.log2(condi_matr), 1)
    return 2 ** ent

#return sigmas corresponding to the perplexity input
def find_sigmas(dists, perplexity):
    found_sigmas = np.zeros(dists.shape[0])
    for i in range(dists.shape[0]):
        func = lambda sig: perp(p_conditional(dists[i:i+1, :], np.array([sig])))
        found_sigmas[i] = search(func, perplexity)
    return found_sigmas

#binary search to find the best sigma
def search(func, goal, tol=1e-10, max_iters=1000, lowb=1e-20, uppb=10000):
    for _ in range(max_iters):
        guess = (uppb + lowb) / 2.
        val = func(guess)

        if val > goal:
            uppb = guess
        else:
            lowb = guess

        if np.abs(val - goal) <= tol:
            return guess

    warnings.warn(f"\nSearch couldn't find goal, returning {guess} with value {val}")
    return guess

#calculate similarity scores for the data in low dimensions with t-distribution
def tdist(y):
    dists = pairwise_distances(y)
    nom = 1 / (1 + dists)
    np.fill_diagonal(nom, 0.)
    return nom / np.sum(np.sum(nom))

#calculate cost function gradient (verify preservation of high dimensionality in lower dimension) 
def gradient(P, Q, y):
    (n, no_dims) = y.shape
    pq_diff = P - Q
    y_diff = np.expand_dims(y,1) - np.expand_dims(y,0)

    dists = pairwise_distances(y)
    aux = 1 / (1 + dists)
    return 4 * (np.expand_dims(pq_diff, 2) * y_diff * np.expand_dims(aux,2)).sum(1)

#step function to speed up learning 
def m(t):
    return 0.5 if t < 250 else 0.8

#calculate original similarity score using gaussian distribution
def gaussian(X, perp):
    N = X.shape[0]
    dists = pairwise_distances(X)
    sigmas = find_sigmas(dists, perp)
    g_cond = g_conditional(dists, sigmas)
    return (g_cond + g_cond.T) / (2. * N)

#reduce dimensions of input dataset
def tsKNEE(adata, T=1000, perp = 30):
    if "leiden" not in adata.obs: 
        raise Exception("Anndata object is not clustered with Leiden or the Leiden cluster values are not stored in a column named leiden in adata.obs.")
    X = adata.X.toarray()
    N = X.shape[0]
    P = gaussian(X, perp)
    learning_rate=500
    ydim=2
    Y = []
    y = np.random.normal(loc=0.0, scale=1e-4, size=(N,ydim))
    Y.append(y); Y.append(y)

    for t in range(T):
        Q = tdist(Y[-1])
        grad = gradient(P, Q, Y[-1])
        momentum = m(t)
        y = Y[-1] - learning_rate*grad + momentum*(Y[-1] - Y[-2])
        Y.append(y)
        if t % 10 == 0:
            Q = np.maximum(Q, 1e-12)
    adata.obsm['X_tsne'] =  y

#plot dataset
def tsKNEE_plot(adata, xlabel = "tsne1", ylabel = "tsne 2", title = "", save = None):
    # import the n_obs x 2 (coordinates)
    if "leiden" not in adata.obs: 
        raise ValueError("Anndata object is not clustered with Leiden or the Leiden cluster values are not stored in a column named leiden in adata.obs.")
    if "X_tsne" not in adata.obsm:
        raise ValueError("Run tsKNEE on anndata object before continuing")
    if save is not None and not isinstance(save, str):
        raise ValueError("parameter for save needs to be a string or None")
    if not (isinstance(xlabel, str) or isinstance(ylabel, str) or isinstance(title, str)):
        raise ValueError("parameter for xlabel, ylabel, title need to be a string")
    x = [val[0] for val in adata.obsm['X_tsne']]
    y = [val[1] for val in adata.obsm['X_tsne']] 
    colormap = list(mcolors.CSS4_COLORS.values())
    colormap = random.sample(colormap, len(set(adata.obs["leiden"])))
    fig, ax = plt.subplots()
    ax.scatter(x=x,y=y, c = [colormap[int(i)] for i in adata.obs["leiden"]])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if save is not None: 
        fig.savefig(save)
    plt.show()