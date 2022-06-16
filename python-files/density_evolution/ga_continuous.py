# %% Imports and variables
from time import process_time
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
import density_evolution
import de_utils
import random

from config import cfg
cfg_de = cfg.get('density_evolution')
cfg_cont = cfg_de.get('ga_continous')


def start_point(R, dv, dc):
    # lam and rho are assumed to begin at degree 2
    k = 1/np.arange(2, dc+1)
    i = 1/np.arange(2, dv+1)

    # C is the constraint matrix.
    # 1st row: code rate constraint
    # 2nd row: coefficients of rho should sum to 1
    # 3rd row: coefficients of lam should sum to 1
    # lb: coefficients of rho and lam should be positive (>0)
    c1 = np.concatenate((k, -(1-R)*i))
    c2 = np.concatenate((np.ones(dc-1), np.zeros(dv-1)))
    c3 = np.concatenate((np.zeros(dc-1), np.ones(dv-1)))
    C = np.stack((c1, c2, c3))

    d = np.array([[0], [1], [1]])
    lb = np.zeros((dv+dc-2, 1))
    f = np.ones((dv+dc-2, 1))

    sol = linprog(f, A_ub=None, b_ub=None, A_eq=C,
                  b_eq=d, method='interior-point')
    x_0 = sol.x

    return C, x_0


def add_one_degree(C_c, x_0, dc):

    # %% Add degree 1
    x_0 = np.concatenate(([0], x_0[:dc-1], [0], x_0[dc-1:]))
    D = np.size(C_c, 1)
    z_vec = np.zeros((1, D))
    C_c = np.vstack((z_vec, C_c[:dc-1, :], z_vec, C_c[dc-1:, :]))

    return C_c, x_0

# Differential evolution help functions


def compute_x(x_0, C_c, eta):
    return x_0 + np.matmul(C_c, eta)


def poly_eval(p, x):
    exp = np.arange(0, p.size)
    return p*x**exp


def compute_cost(x, dc):
    cost = 0
    in_domain = True
    for xi in x:
        if xi < 0:
            cost += 100 - 10*xi
            in_domain = False

    if cost == 0:
        rho = x[:dc]
        lam = x[dc:]
        min = 1e-5
        max = 0.4
        result = density_evolution.bisection_search(
            min, max, lambda x: eval(x, rho, lam))
        cost = -result

    return cost, in_domain


def eval(x, rho, lam):
    n_grid = 256
    sigma, p0, bins = de_utils.compute_pdf(
        x, 0.5, n_grid, 30)
    f_grid, g_grid, pdf = de_utils.create_pdf(
        p0, bins, n_grid)
    cdf = de_utils.to_cdf(pdf)

    result = density_evolution.symmetric_density_evolution(
        cdf, f_grid, g_grid, rho, lam, plot=False, prnt=False)

    return result


def save_population(eta_Np):
    with open('eta_Np.npy', 'wb') as f:
        np.save(f, eta_Np)


def save_best_individual(C_c, x_0, eta_Np, cost_Np):
    with open('x_best.npy', 'wb') as f:
        idx = np.argmin(cost_Np)
        eta_best = eta_Np[:, idx]
        x_best = compute_x(x_0, C_c, eta_best)
        np.save(f, x_best)

# Optimization


def differential_evolution(C_c, x_0, dc, prnt=True):
    # Np: number of individuals in population
    # gens: Number of generations. Around 5000
    # F: mutation variable. Usually in interval [0.1,1].
    # Cr: recombination probability, [0,1].

    Np = cfg_de.get('Np')
    gens = cfg_de.get('gens')
    F = cfg_de.get('F')
    Cr = cfg_de.get('Cr')
    D = np.size(C_c, 1)

    D = np.size(C_c, 1)

    if prnt:
        header = f"""
        ===================================================================
        Optimizing code ensamble through differential evolution algorithm.
        Number of individuals: {Np}.
        Number of generations: {gens}.
        Parameter settings: F = {F}, Cr = {Cr}
        -------------------------------------------------------------------
        """
        print(header)
        print("Initialising...")

    # x = x_0 + C_c*eta
    # Initialize population
    eta_Np = np.random.rand(D, Np)
    cost_Np = np.zeros(Np)
    cost_gen = np.zeros(gens)
    domain_Np = np.full(Np, False)

    try:
        # Optimize population over generations
        t1_start = process_time()
        for g in range(gens):

            u_Np = np.copy(eta_Np)
            for i in range(Np):
                # Generate selection parameters according to algorithm
                r0 = random.randint(0, Np-1)
                r1 = random.randint(0, Np-1)
                r2 = random.randint(0, Np-1)
                jrand = random.randint(0, D-1)
                while (r0 == i):
                    r0 = random.randint(0, Np-1)
                while (r1 == i or r1 == r0):
                    r1 = random.randint(0, Np-1)
                while (r2 == i or r2 == r0 or r2 == r1):
                    r2 = random.randint(0, Np-1)

                # Mutation and cross over
                cr_arr = np.random.rand(D)
                for j in range(D):
                    if (cr_arr[j] < Cr or j == jrand):
                        u_Np[j, i] = eta_Np[j, r0] + F * \
                            (eta_Np[j, r1]-eta_Np[j, r2])

                # Selection
                if g == 0:
                    x = compute_x(x_0, C_c, eta_Np[:, i])
                    cost, in_domain = compute_cost(x, dc)
                    cost_Np[i] = cost
                    domain_Np[i] = in_domain

                x = compute_x(x_0, C_c, u_Np[:, i])
                cost, in_domain = compute_cost(x, dc)

                if cost < cost_Np[i]:
                    eta_Np[:, i] = u_Np[:, i]
                    cost_Np[i] = cost
                    domain_Np[i] = in_domain

            cost_gen[g] = np.sum(cost_Np)/Np

            # Write current status
            if prnt:
                n_domain = domain_Np.sum()
                if n_domain != 0:
                    min_cost = np.min(cost_Np[domain_Np])
                    rber_best = -min_cost
                else:
                    rber_best = 0
                t1_stop = process_time()
                status = f"{g} generations completed.  Average cost:{cost_gen[g]:.2f}.  Individuals in domain: {n_domain}. Best RBER: {rber_best}"
                print(status, end='\r', flush=True)

        status = f"""
        -------------------------------------------------------------------
        Finished!
        ===================================================================
            """
        print(status)
    finally:
        save_best_individual(C_c, x_0, eta_Np, cost_Np)

        save_population(eta_Np)


def ga_continous(print_terminal):

    R = cfg_de.get('R')
    dv = cfg_de.get('dv')
    dc = cfg_de.get('dc')

    # x = [rho;lam]
    # 1. Compute an initial start point
    C, x_0 = start_point(R, dv, dc)

    # 2. Find the complement of constraint matrix.
    # This gives us the allowed directions to step around in
    C_c = null_space(C)

    # 3. Add degree 1 to x_0 and C_c matrix
    C_c, x_0 = add_one_degree(C_c, x_0, dc)

    # 4. Perform optimization through differential evolution
    differential_evolution(C_c, x_0, dc, prnt=True)
