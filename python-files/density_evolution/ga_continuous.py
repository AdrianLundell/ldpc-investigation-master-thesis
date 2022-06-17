# %% Imports and variables
from time import process_time
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
import de_utils as de_u
import random

from config import cfg
cfg_de = cfg.get('density_evolution')
cfg_cont = cfg_de.get('ga_continuous')
run_id = cfg.get("run_id")


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


def compute_fitness(x, dc):
    fitness = 0
    in_domain = True
    for xi in x:
        if xi < 0:
            fitness += -100 + 10*xi
            in_domain = False

    if fitness == 0:
        rho = x[:dc]
        lam = x[dc:]
        min = cfg_de.get("min_rber")
        max = cfg_de.get("max_rber")
        result = de_u.bisection_search(
            min, max, rho, lam)
        fitness = 100*result

    return fitness, in_domain


def differential_evolution(C_c, x_0, dc):
    # Np: number of individuals in population
    # gens: Number of generations. Around 5000
    # F: mutation variable. Usually in interval [0.1,1].
    # Cr: recombination probability, [0,1].

    Np = cfg_de.get('Np')
    gens = cfg_de.get('generations')
    F = cfg_cont.get('F')
    Cr = cfg_cont.get('Cr')
    D = np.size(C_c, 1)
    load_population = cfg_de.get("load_population")
    de_u.save_params()

    print_terminal = cfg_de.get("print_terminal")
    if print_terminal:
        header = f"""
        ga_continuous.py
        ===================================================================
        Optimizing code ensamble through differential evolution algorithm.
        Number of individuals: {Np}.
        Number of generations: {gens}.
        Parameter settings: F = {F}, Cr = {Cr}
        -------------------------------------------------------------------
        """
        print(header)

    # x = x_0 + C_c*eta
    # Initialize population
    if load_population:
        fname = "data/" + run_id + ".npz"
        data = np.load(fname)
        eta_Np = data["population"]
        fitness_Np = data["fitness"]
    else:
        eta_Np = np.random.rand(D, Np)
        fitness_Np = np.zeros(Np)

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
                    fitness, in_domain = compute_fitness(x, dc)
                    fitness_Np[i] = fitness
                    domain_Np[i] = in_domain

                x = compute_x(x_0, C_c, u_Np[:, i])
                fitness, in_domain = compute_fitness(x, dc)

                if fitness > fitness_Np[i]:
                    eta_Np[:, i] = u_Np[:, i]
                    fitness_Np[i] = fitness
                    domain_Np[i] = in_domain

            if g % 10 == 0:
                de_u.save_fitness(fitness_Np, g)

            ave_fitness = np.sum(fitness_Np)/Np

            # Write current status
            if print_terminal:
                n_domain = domain_Np.sum()
                if n_domain != 0:
                    max_fitness = np.max(fitness_Np[domain_Np])
                    rber_best = max_fitness
                else:
                    rber_best = 0
                t1_stop = process_time()
                status = f"{g} generations completed.  Average fitness:{ave_fitness:.2f}.  Individuals in domain: {n_domain}. Best RBER: {rber_best}"
                print(status, end='\r', flush=True)

        status = f"""
        -------------------------------------------------------------------
        Finished!
        ===================================================================
            """
        print(status)
    finally:
        de_u.save_population(eta_Np, fitness_Np)


def ga_continous():

    R = cfg_de.get('R')
    dv = cfg_cont.get('dv')
    dc = cfg_cont.get('dc')

    # x = [rho;lam]
    # 1. Compute an initial start point
    C, x_0 = start_point(R, dv, dc)

    # 2. Find the complement of constraint matrix.
    # This gives us the allowed directions to step around in
    C_c = null_space(C)

    # 3. Add degree 1 to x_0 and C_c matrix
    C_c, x_0 = add_one_degree(C_c, x_0, dc)

    # 4. Perform optimization through differential evolution
    differential_evolution(C_c, x_0, dc)
