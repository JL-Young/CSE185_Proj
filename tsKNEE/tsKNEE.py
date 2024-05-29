import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np 
import scanpy as sc, anndata as ad

#finding the euclidean distance 
def pairwise_distances(X):
    return np.sum((X[None, :] - X[:, None])**2, 2)

#normalized gaussian distribution
def p_conditional(dists, sigmas):
    e = np.exp(-dists / (2 * np.square(sigmas.reshape((-1,1)))))
    np.fill_diagonal(e, 0.)
    e += 1e-8
    return e / e.sum(axis=1).reshape([-1,1])

#the standard deviation which is determined by the perplexity 
def perp(condi_matr):
    ent = -np.sum(condi_matr * np.log2(condi_matr), 1)
    return 2 ** ent

#sigmas corresponding to the perplexities 
def find_sigmas(dists, perplexity):
    found_sigmas = np.zeros(dists.shape[0])
    for i in range(dists.shape[0]):
        func = lambda sig: perp(p_conditional(dists[i:i+1, :], np.array([sig])))
        found_sigmas[i] = search(func, perplexity)
    return found_sigmas

#binary search
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

#utilizing a t distribution to reduce the clumping of the data 
def q_joint(y):
    dists = pairwise_distances(y)
    nom = 1 / (1 + dists)
    np.fill_diagonal(nom, 0.)
    return nom / np.sum(np.sum(nom))

#calculating the gradient of the cost function to see if the high dimensionality is preserved in the lower dimension 
def gradient(P, Q, y):
    (n, no_dims) = y.shape
    pq_diff = P - Q
    y_diff = np.expand_dims(y,1) - np.expand_dims(y,0)

    dists = pairwise_distances(y)
    aux = 1 / (1 + dists)
    return 4 * (np.expand_dims(pq_diff, 2) * y_diff * np.expand_dims(aux,2)).sum(1)

#step function to make the learning quicker 
def m(t):
    return 0.5 if t < 250 else 0.8

#dealing with the crowding problem
def p_joint(X, perp):
    N = X.shape[0]
    dists = pairwise_distances(X)
    sigmas = find_sigmas(dists, perp)
    p_cond = p_conditional(dists, sigmas)
    return (p_cond + p_cond.T) / (2. * N)

def tsKNEE(adata, T=1000, perp = 30):
    if "leiden" not in adata.obs: 
        raise Exception("Anndata object is not clustered with Leiden or the Leiden cluster values are not stored in a column named leiden in adata.obs.")
    X = adata.X.toarray()
    N = X.shape[0]
    P = p_joint(X, perp)
    l=500
    ydim=2
    Y = []
    y = np.random.normal(loc=0.0, scale=1e-4, size=(N,ydim))
    Y.append(y); Y.append(y)

    for t in range(T):
        Q = q_joint(Y[-1])
        grad = gradient(P, Q, Y[-1])
        y = Y[-1] - l*grad + m(t)*(Y[-1] - Y[-2])
        Y.append(y)
        if t % 10 == 0:
            Q = np.maximum(Q, 1e-12)
    adata.obsm['X_tsne'] =  y

# can have the coordinates already generated, or just have a adata object (can run tsne in this function)
def tsKNEE_plot(adata, perp = 30, xlabel = "tsne1", ylabel = "tsne 2", title = "", save = None):
    # import the n_obs x 2 (coordinates)
    if "leiden" not in adata.obs: 
        raise Exception("Anndata object is not clustered with Leiden or the Leiden cluster values are not stored in a column named leiden in adata.obs.")
    if "X_tsne" not in adata.obsm:
        raise Exception("Run tsKNEE on anndata object before continuing")
    x = [val[0] for val in adata.obsm['X_tsne']]
    y = [val[1] for val in adata.obsm['X_tsne']]
    plt.scatter(x=x,y=y, c = adata.obs["leiden"])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

    if save is not None: 
        plt.savefig(save)



