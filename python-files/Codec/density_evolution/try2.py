#%%
import numpy as np 
import matplotlib.pyplot as plt 

#%% Design grid to fit llr points
#Computed llrs from optimal_thresholds, sigma1=sigma2=1, mu1=-1, mu2=1
llrs = [-2.87585378, -0.81106086, 0.81106086, 2.87585378]

c = np.polyfit([-2,-1,1,2], llrs, 1)
step_size = (llrs[1] - llrs[0])/(len(llrs)*15)
n = 2**9
f_grid = np.arange(-n//2+1, n//2+1)*step_size
f_grid = f_grid + np.min(np.abs(llrs[0] - f_grid))

for llr in llrs:
    print(np.min(np.abs(llrs[0] - f_grid)))
#%%
points = np.log((1+np.exp(-f_grid))/(1-np.exp(-f_grid)))
plt.scatter(points, f_grid, s=1)
plt.xlabel("G points")
plt.ylabel("F points")

# %%
plt.hist(points)
# %%
f_grid = 4
np.log((1+np.exp(-f_grid))/(1-np.exp(-f_grid)))

# %%
# Gamma test
#GAMMA(distribution) is the distribution of gamma(messages)
#Simulated density and computed density should look the same
m = -10 + 20*np.random.rand(10**5)

n = 2**10
m_cdf = np.linspace(0,1, n)
F_grid = np.linspace(-10, 20, n)
G_grid = np.linspace(0, 8, n//2)
G_computed = gamma(m_cdf, F_grid, G_grid)

G0_sampled = -np.log(np.tanh(np.abs(m[m<=0])/2))
G0_hist, bins = np.histogram(G0_sampled, bins=G_grid)
G1_sampled = -np.log(np.tanh(np.abs(m[m>0])/2))
G1_hist, bins = np.histogram(G1_sampled, bins=G_grid)

plt.plot(G_grid, to_pdf(G_computed[0,:])*10**5)
plt.plot(G_grid[:-1], G0_hist)
plt.show()
plt.plot(G_grid, to_pdf(G_computed[1,:])*10**5, "--")
plt.plot(G_grid[:-1], G0_hist, "--")
plt.show()

 #%%Gamma inverse should yield same result back
plt.plot(m_cdf)
plt.plot(gamma_inv(G, F_grid, G_grid))
plt.show()

#%% Test of probability convertion functions
N = 10
unif_cdf = np.linspace(1/N, 1, N, endpoint=True)
unif_pdf = np.full(N, 1/N)

plt.plot(to_cdf(unif_pdf))
plt.plot(unif_cdf)
plt.show()

plt.plot(to_pdf(unif_cdf))
plt.plot(unif_pdf)
plt.show()
#%%
import numpy as np 

def count_edges(coeffs):
    x = np.arange(1, len(coeffs) + 1)
    return coeffs @ x

n_vn = 8
n_cn = 4
L = np.array([0,0,0,0,n_vn])

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
    pass

"""
Idea: 
    Calculate "Base vectors" of rho and lambda with one in one place and zero in all others.
    Create new individuals as linear combinations of base vectors (making sure to get only positive coefficients)
"""
def edge_perspective(coeffs, n):
    """
    Computes the edge perspective coefficients from degree counts
    """
    x = np.arange(1, len(coeffs) + 1) * coeffs
    x = x / np.sum(x)

    return x
    
def mix(x1, x2):

    x1_bound = np.min(x2 / x1)
    x2_bound = np.min(x1 / x2)
    coeff = -x1_bound + np.random.rand()*np.abs(x1_bound - x2_bound)
    


n_vn = 8
n_cn = 4
L = np.array([0,0,0,0, n_vn])