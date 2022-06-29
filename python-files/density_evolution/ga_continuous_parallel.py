# %% Imports and variables
from time import process_time
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
import de_utils as de_u
import random
from multiprocessing import Pool, shared_memory
from multiprocessing.sharedctypes import Value


from config import cfg
cfg_de = cfg.get('density_evolution')
cfg_cont = cfg_de.get('ga_continuous')
run_id = cfg.get("run_id")
n_processes = cfg.get("n_processes")


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

# Function that is parallelized through Pool in multiprocessing
def de_individual(i):
    global eta_name, u_name, fit_name, dom_name, eta_shape, u_shape, fit_shape, dom_shape, Np, Cr, F, C_c, x_0, g_idx, dc
    D = np.size(C_c,1)
    shm_eta = shared_memory.SharedMemory(name = eta_name)
    shm_u = shared_memory.SharedMemory(name=u_name)
    shm_fit = shared_memory.SharedMemory(name=fit_name)
    shm_dom = shared_memory.SharedMemory(name=dom_name)

    eta = np.ndarray(eta_shape, dtype=np.float64, buffer=shm_eta.buf)
    u = np.ndarray(u_shape, dtype=np.float64, buffer=shm_u.buf)
    fitness = np.ndarray(fit_shape, dtype=np.float64, buffer=shm_fit.buf)
    domain = np.ndarray(dom_shape, dtype=bool, buffer=shm_dom.buf)

    g = g_idx.value

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
            u[j, i] = eta[j, r0] + F * \
                (eta[j, r1]-eta[j, r2])

    # Selection
    if g_idx.value == 0:
        x = compute_x(x_0, C_c, eta[:, i])
        fit, in_domain = compute_fitness(x, dc)
        fitness[g, i] = fit
        domain[i] = in_domain
    else:
        x = compute_x(x_0, C_c, u[:, i])
        fit, in_domain = compute_fitness(x, dc)

        if fit > fitness[g-1, i]:
            eta[:, i] = u[:, i]
            fitness[g, i] = fit
            domain[i] = in_domain
        else:
            fitness[g, i] = fitness[g-1, i]

    shm_eta.close()
    shm_u.close()
    shm_fit.close()
    shm_dom.close()


# Used for initializeing multiprocessing
def pool_initializer():
    global eta_name, u_name, fit_name, dom_name, eta_shape, u_shape, fit_shape, dom_shape, Np, Cr, F, C_c, x_0, g_idx, dc



def differential_evolution(C_c, x_0, dc):
    global eta_name, u_name, fit_name, dom_name, eta_shape, u_shape, fit_shape, dom_shape, Np, Cr, F, g_idx
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

    code_rate = cfg_cont.get("R")
    dv = cfg_cont.get("dv")
    dc = cfg_cont.get("dc")

    print_terminal = cfg_de.get("print_terminal")
    header = f"""
    Running ga_continuous.py
    ===================================================================
    Optimizing code ensamble through differential evolution algorithm.
    Code rate: {code_rate:.2f}.
    Number of individuals: {Np}.
    Number of generations: {gens}.
    Max variable node degree: {dv}.
    Max check node degree: {dc}.
    Parameter settings: F = {F}, Cr = {Cr}
    -------------------------------------------------------------------
    """
    if print_terminal:
        print(header)
    else:
        de_u.log(header,'w')

    # x = x_0 + C_c*eta
    # Initialize population
    if load_population:
        fname = "data/" + run_id + ".npz"
        data = np.load(fname)
        eta_Np = data["population"]
        fitness_Np = data["fitness"]
        g_start = int(data["generation"][0])+1
        dim_0 = np.size(fitness_Np, axis=0)
        if gens != np.size(fitness_Np, axis=0):
            fitness_Np_new = np.zeros((gens, Np))
            fitness_Np_new[:dim_0] = fitness_Np
            fitness_Np = fitness_Np_new
    else:
        eta_Np = np.random.rand(D, Np)
        fitness_Np = np.zeros((gens,Np))
        g_start = 0

    domain_Np = np.full(Np, False)

    # Initialize shared_memory in multiprocessing
    shm_eta = shared_memory.SharedMemory(create=True, size=eta_Np.nbytes)
    shm_u = shared_memory.SharedMemory(create=True, size=eta_Np.nbytes)
    shm_fit = shared_memory.SharedMemory(create=True, size=fitness_Np.nbytes)
    shm_dom = shared_memory.SharedMemory(create=True, size=domain_Np.nbytes)
    
    eta = np.ndarray(eta_Np.shape, dtype=np.float64, buffer=shm_eta.buf)
    u = np.ndarray(eta_Np.shape, dtype=np.float64, buffer=shm_u.buf)
    fitness = np.ndarray(fitness_Np.shape, dtype=np.float64, buffer=shm_fit.buf)
    domain = np.ndarray(domain_Np.shape, dtype=bool, buffer=shm_dom.buf)

    eta[:] = eta_Np[:]
    u[:] = eta_Np[:]
    fitness[:] = fitness_Np[:]
    domain[:] = domain_Np[:]

    eta_name = shm_eta.name
    u_name = shm_u.name
    fit_name =  shm_fit.name
    dom_name = shm_dom.name
    
    eta_shape = eta.shape
    u_shape = u.shape 
    fit_shape = fitness.shape
    dom_shape = domain.shape

    g_idx = Value('i', g_start)

    # Use Pool for parallel computing
    pool = Pool(processes=n_processes, initializer=pool_initializer)
    pop_idx = range(Np)
    
    try:
        # Optimize population over generations
        t1_start = process_time()
        for g in range(g_start,gens):

            u[:] = np.copy(eta)
            pool.map(de_individual, pop_idx)
            # Equivalent to:
            #for i in range(Np):
            #    de_individual(i)
            g_idx.value += 1
            
            ave_fitness = np.sum(fitness[g,:])/Np

            # Write current status
            
            n_domain = domain.sum()
            if n_domain != 0:
                rber_max = np.max(fitness[g,domain])
                rber_min = np.min(fitness[g,domain])
            else:
                rber_max = 0
                rber_min = 0
            t1_stop = process_time()
            status = f"{g} generations completed.  Average fitness:{ave_fitness:.2f}.  In domain: {n_domain}. RBER: max {rber_max:.2f}, min: {rber_min:.2f}."
            if print_terminal:    
                print(status)
            else:
                de_u.log(status,'a')

            if g % int(cfg_de.get("save_interval")) == 0:
                de_u.save_population(eta, fitness, g)

        status = f"""
    -------------------------------------------------------------------
    Finished!
    ===================================================================
            """
        if print_terminal:
            print(status)
        else:
            de_u.log(status,'a')
    finally:
        de_u.save_population(eta_Np, fitness_Np, g)
        shm_eta.close()
        shm_u.close()
        shm_fit.close()
        shm_dom.close()
        shm_eta.unlink()
        shm_u.unlink()
        shm_fit.unlink()
        shm_dom.unlink()


def ga_continous_parallel():
    global C_c, x_0, dc

    R = cfg_cont.get('R')
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

# %%

