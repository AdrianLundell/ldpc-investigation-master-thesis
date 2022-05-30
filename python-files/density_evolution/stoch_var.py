import numpy as np 
import matplotlib.pyplot as plt 

class stoch_var:
    """Representation of a discrete stochastic distribution with 
    densities in +-inf"""

    pdf : np.array
    cdf : np.array
    grid : np.array

    def __init__(self, pdf, grid):

        self.pdf = pdf
        self.grid = grid
        self.step = abs(grid[2] - grid[1])
        self.cdf = to_cdf(pdf)

    def plot(self):
        """Plots the pdf and cdf of variable"""

        grid = self.grid[1:-1]
        grid = np.hstack([])

        fig, axes = plt.subplots(2,1)
        axes[0].scatter(self.grid, self.pdf)
        axes[0].set_title("PDF")
        axes[1].step(self.grid, self.cdf)
        axes[1].set_title("CDF")
        plt.show()

def to_cdf(pdf):
    """Returns the discrete pdf of a cdf"""
    return np.cumsum(pdf)

def to_pdf(cdf):
    """Returns the discrete cdf of a pdf"""
    cdf = np.hstack((0, cdf, 1))
    pdf = cdf[1:] - cdf[:-1]
    return pdf[:-1]