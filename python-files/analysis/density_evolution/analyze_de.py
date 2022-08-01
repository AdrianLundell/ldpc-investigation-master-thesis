import numpy as np
import matplotlib.pyplot as plt
import sys
import yaml

sys.path.append("../../density_evolution")
import de_utils as de_u
config_file = "../../density_evolution/config.yml"
with open(config_file, "r") as ymlfile:
    cfg = yaml.safe_load(ymlfile)

cfg["density_evolution"]["dmc_file"] = "../../compute_dmc/data/4_04_AWGN.csv"
de_u.set_cfg(cfg)

cfg_de = cfg.get('density_evolution')
cfg_cont = cfg_de.get('ga_continuous')
cfg_disc = cfg_de.get('ga_discrete')

######################################
# Discrete analysis
######################################
def node_to_edge(lam_node, rho_node):
    rho_edge = np.arange(1, len(rho_node)+1) * rho_node / \
        np.sum(np.arange(1, len(rho_node)+1) * rho_node)
    lam_edge = np.arange(1, len(lam_node)+1) * lam_node / \
        np.sum(np.arange(1, len(lam_node)+1) * lam_node)
    return lam_edge, rho_edge        
    
# Changeable params
fname = "../../density_evolution/data/optimal_soft.npz"
min_rber = 0.005  # cfg_de.get("min_rber")
max_rber = 0.0  # cfg_de.get("max_rber")
step_rber = 0.001
plot_rber = 0.02
plot_pdf = False
plot_iteration = 0
plot_convergence = False


# Load data
data = np.load(fname)
population = data["population"]
fitness = data["fitness"]
last_gen = int(data["generation"][0])+1
best_idx = int(data["best_idx"][0])
lam_node = data["lam_node"]
rho_node = data["rho_node"]
best_rber = data["best_rber"][0]
lam_edge, rho_edge = node_to_edge(lam_node,rho_node)
"""  
def get_degree_distrubutions(i):
    cn_degrees = np.sum(i, 0)
    vn_degrees = np.sum(i, 1)
    rho_node = np.bincount(vn_degrees)[1:]
    rho_edge = np.arange(1, len(rho_node)+1) * rho_node / \
        np.sum(np.arange(1, len(rho_node)+1) * rho_node)
    lam_node = np.bincount(cn_degrees)[1:]
    lam_edge = np.arange(1, len(lam_node)+1) * lam_node / \
        np.sum(np.arange(1, len(lam_node)+1) * lam_node)
    rho_node = rho_node/np.sum(rho_node)
    lam_node = lam_node/np.sum(lam_node)
    return lam_edge, rho_edge
"""
def convert_edge_to_node(p_edge):
    p_node = p_edge/(np.arange(1, len(p_edge)+1) *
                    np.sum(1/np.arange(1, len(p_edge)+1)*p_edge))
    return p_node
def symmetric_density_evolution(cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, tol, n_iter=cfg_de.get("de_iter"),  plot=plot_pdf):

    pl = cdf
    p0 = cdf
    pl_old = 0
    i = 0
    diff = np.inf
    error = np.inf
    cbp = np.zeros(n_iter)
    
    while True:
        pdf_l = de_u.to_pdf(pl)
        pdf_l = np.where(pdf_l > 1e-10, pdf_l, 0)
        prob_l = np.sum(pdf_l)

        if i == 0 or i == 1:
            coeff = pl.size/f_grid.size
            start = coeff*f_grid[0]
            step = f_grid[1]-f_grid[0]
            f_grid_l = np.arange(start, -start, step)
            if i == 0:
                rber = np.sum(pdf_l[:int(pdf_l.size/2)])
        
        cbp[i] = de_u.compute_cbp(pdf_l, f_grid_l,False)

        if plot==True and i == plot_iteration and np.abs(rber-plot_rber)<step_rber:
            
            title = f"Iteration: {i}. RBER={rber*100:.3f}%. CBP_l = {cbp[i]}. PDF total probability = {prob_l:.3f}."
            fig, (ax1, ax2) = plt.subplots(2)

            ax1.plot(f_grid_l, pl)
            ax1.set_xlabel("Message LLR")
            ax1.set_ylabel("CDF")
            ax1.set_title(title)

            ax2.plot(f_grid_l, pdf_l)
            ax2.set_xlabel("Message LLR")
            ax2.set_ylabel("PDF")
            ax2.set_title(title)

            fig.tight_layout()
            plt.show()
            fig.clf()
            #plt.close()


        
        is_zero = (cbp < tol)
        max_iter = (i == n_iter-1)
        if max_iter:
            return cbp


        x1 = de_u.gamma(pl, f_grid, g_grid)
        x2 = de_u.rho(x1, rho_coeffs)
        x3 = de_u.gamma_inv(x2, f_grid, g_grid)
        x4 = de_u.lambd(x3, lambda_coeffs)
        pl = de_u.conv(x4, p0)

        diff = sum((pl_old - pl)**2)
        pl_old = pl

        zero_index = pl.size//2
        error = pl[zero_index]
        # Stability criterion
        
        
        i += 1




if plot_convergence:
    best_idx = np.argmax(fitness[last_gen, :])
    #lam_edge, rho_edge = get_degree_distrubutions(population[best_idx, :, :])
    #print(f"Lambda edge degree distribution: \n {lam_edge}")
    #print(f"Rho edge degree distribution: \n {rho_edge}")
    print(f"Lambda node degree distribution: \n {lam_node}")
    print(f"Rho node degree distribution: \n {rho_node}")
    print(f"Best RBER: {best_rber}")

    rbers = np.arange(min_rber,max_rber,step_rber)
    epsilons = np.zeros(rbers.size)
    cbp_0 = np.zeros(rbers.size)
    cbp_l = np.zeros(rbers.size)

    

    # Analyze evolution of thresholds 
    for i,rber in enumerate(rbers):
        n_grid = cfg_de.get("n_grid")
        f_grid, g_grid, pdf = de_u.init_pdf(rber, n_grid)
        cdf = de_u.to_cdf(pdf)

        eps = de_u.compute_threshold(pdf, f_grid, lam_edge, rho_edge)
        cbp_0[i] = de_u.compute_cbp(pdf,f_grid,True)
        epsilons[i] = eps

        cbp = symmetric_density_evolution(
            cdf, f_grid, g_grid, rho_edge, lam_edge, tol=1)
        if np.abs(rber-plot_rber) < step_rber:
            cbp_l = cbp
            rber_l =rber
            eps_l = eps

        print(rber)

    plt.subplot(3, 1, 1)
    plt.plot(rbers,cbp_0)
    plt.legend(["CBP_0"])
    plt.xlabel("RBER")

    plt.subplot(3, 1, 2)
    plt.plot(rbers, epsilons)
    plt.legend(["Epsilon"])
    plt.xlabel("RBER")

    plt.subplot(3,1,3)
    iters = np.arange(cbp_l.size)
    plt.plot(iters,cbp_l)
    plt.plot(iters,np.full(cbp_l.size,eps_l),'-')
    legend1 = f"CBP_l when RBER={rber_l:.3f}."
    legend2 = f"epsilon={eps_l}."
    plt.legend([legend1, legend2])
    plt.xlabel("Iterations")

    plt.tight_layout()
    plt.show()


rber_lim = de_u.bisection_search(min_rber, max_rber, rho_edge, lam_edge)
message = f"The boundary for RBER is {rber_lim}"
print(message)
