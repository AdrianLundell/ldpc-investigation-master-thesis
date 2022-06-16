def count_edges(coeffs):
    x = np.arange(1, len(coeffs) + 1)
    return coeffs @ x

def valid_poly(coeffs, n):
    """
    Calculates coeffs matching code rate by distributing edges evenly over n nodes
        coeffs: A vector of degree counts,
        n: The number of nodes to distribute the edges over
    """

    tot_edges = count_edges(coeffs)
    x = np.zeros(n, int)
    for i in range(tot_edges):
        j = np.mod(i, (n))
        x[j] += 1

    coeffs = np.bincount(x)

    return coeffs[1:]

def edge_perspective(coeffs, n):
    """
    Computes the edge perspective coefficients from degree counts
    """
    x = np.arange(1, len(coeffs) + 1) * coeffs
    x = x / np.sum(x)

    return x
    
#%%
