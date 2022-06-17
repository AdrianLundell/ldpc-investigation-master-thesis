import numpy as np
import de_utils as de_u
import matplotlib.pyplot as plt

from config import cfg
cfg_de = cfg.get('density_evolution')
cfg_disc = cfg_de.get('ga_disc')
run_id = cfg.get("run_id")


def init_population(n_pop, n_cn, n_vn):
    population = np.zeros((n_pop, n_cn, n_vn), int)
    for i in range(n_pop):
        for j in range(n_vn):
            k = np.random.choice(n_cn, 2, replace=False)
            population[i, k, j] = 1

    return population


def vertical_crossover(i1, i2):
    crossover_point = np.random.randint(i1[0, :].size)
    i1_new = np.hstack((i1[:, :crossover_point], i2[:, crossover_point:]))
    i2_new = np.hstack((i1[:, crossover_point:], i2[:, :crossover_point]))
    return i1_new, i2_new


def horizontal_crossover(i1, i2):
    crossover_point = np.random.randint(i1[:, 0].size)
    i1_new = np.vstack((i1[:crossover_point, :], i2[crossover_point:, :]))
    i2_new = np.vstack((i1[crossover_point:, :], i2[:crossover_point, :]))
    return i1_new, i2_new


def mutation(i):
    mutation_i = np.random.randint(0, i[:, 0].size)
    mutation_j = np.random.randint(0, i[0, :].size)
    i[mutation_i, mutation_j] = not i[mutation_i, mutation_j]
    return i


def tournament(population, fitness, n_competitiors=3):
    n_pop = cfg_de.get("Np")
    competitors_index = np.random.choice(n_pop, n_competitiors, replace=False)
    competitors_fitness = fitness[competitors_index]
    winner_index = competitors_index[np.argmax(competitors_fitness)]
    i_new = population[winner_index, :, :]
    fitness_new = fitness[winner_index]
    return i_new, fitness_new


def evaluate(i):
    cn_degrees = np.sum(i, 0)
    vn_degrees = np.sum(i, 1)

    condition1 = np.all(cn_degrees >= 2) and np.all(vn_degrees >= 2)
    condition2 = np.all(cn_degrees <= 25)

    if condition1 and condition2:
        rho_node = np.bincount(vn_degrees)[1:]
        rho_edge = np.arange(1, len(rho_node)+1) * rho_node / \
            np.sum(np.arange(1, len(rho_node)+1) * rho_node)
        lam_node = np.bincount(cn_degrees)[1:]
        lam_edge = np.arange(1, len(lam_node)+1) * lam_node / \
            np.sum(np.arange(1, len(lam_node)+1) * lam_node)

        min = cfg_de.get("min_rber")
        max = cfg_de.get("max_rber")
        fitness = de_u.bisection_search(min, max, rho_edge, lam_edge)
        return fitness

    else:
        return -1


def ga_discrete():
    load_population = cfg_de.get("load_population")

    n_pop = cfg_de.get("Np")
    n_generations = cfg_de.get("generations")
    n_vn = cfg_disc.get("n_vn")
    n_cn = cfg_disc.get("n_cn")
    p_vertical = cfg_disc.get("p_vertical")
    p_horizontal = cfg_disc.get("p_horizontal")
    p_mutation = cfg_disc.get("p_mutation")
    print_terminal = cfg_de.get("print_terminal")

    de_u.save_params()

    if load_population:
        fname = "data/" + run_id + ".npz"
        data = np.load(fname)
        population = data["population"]
        fitness = data["fitness"]
    else:
        population = init_population(n_pop, n_cn, n_vn)
        fitness = np.full(n_pop, -np.inf)

    if print_terminal:
        print(f"""
        ga_disrete.py
        ===================================================================
        Running optimisation of {n_cn}x{n_vn} protograph. 
        Number of individuals: {n_pop}.
        Number of generations: {n_generations}.
        -------------------------------------------------------------------
        """)

    try:
        # Initial evaluation
        for j in range(n_pop):
            if fitness[j] == -np.inf:
                fitness[j] = evaluate(population[j, :, :])

        for i in range(n_generations):

            # Select with elitism
            population_new = np.zeros(population.shape, int)
            fitness_new = np.full(fitness.shape, -np.inf)
            population_new[0] = population[np.argmax(fitness)]
            fitness_new[0] = fitness[np.argmax(fitness)]
            for j in range(1, n_pop):
                population_new[j], fitness_new[j] = tournament(
                    population, fitness)
            population = population_new
            fitness = fitness_new

            # Crossover
            for j in range(1, n_pop-1, 2):
                if np.random.rand() < p_vertical:
                    population[j], population[j +
                                              1] = vertical_crossover(population[j], population[j+1])
                    fitness[j], fitness[j+1] = -np.inf, -np.inf
                if np.random.rand() < p_horizontal:
                    population[j], population[j +
                                              1] = horizontal_crossover(population[j], population[j+1])
                    fitness[j], fitness[j+1] = -np.inf, -np.inf

            # Mutation
            for j in range(1, n_pop):
                if np.random.rand():
                    population[j] = mutation(population[j])
                    fitness[j] = -np.inf

            # Evaluate new individuals
            for j in range(n_pop):
                if fitness[j] == -np.inf:
                    fitness[j] = evaluate(population[j, :, :])

            if i % 10 == 0:
                de_u.save_fitness(fitness, i)

            if print_terminal:
                status = f"Generation {i}/{n_generations}. Best/median/mean/min/variance fitness: {np.max(fitness)}/{np.median(fitness)}/{np.mean(fitness)}/{np.min(fitness)}/{np.var(fitness)}.                                "
                print(status, end='\r', flush=True)

        if print_terminal:
            status = f"""
            -------------------------------------------------------------------
            Finished!
            ===================================================================
                """
            print(status)

    finally:
        de_u.save_population(population, fitness)
