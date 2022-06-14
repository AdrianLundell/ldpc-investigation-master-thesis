#%%
import numpy as np 
import time 

t0 = time.time()
n_cn = 9//3
n_vn = 73//3
n_pop = 50

def init_population(n_pop, n_cn, n_vn):
    population = np.zeros((n_pop, n_cn, n_vn), int)
    for i in range(n_pop):
        for j in range(n_vn):   
            k = np.random.choice(n_cn, 2, replace = False)
            population[i, k, j] = 1
    
    return population

def vertical_crossover(i1, i2):
    crossover_point = np.random.randint(i1[0,:].size)
    i1_new = np.hstack((i1[:,:crossover_point], i2[:,crossover_point:]))
    i2_new = np.hstack((i1[:,crossover_point:], i2[:,:crossover_point]))
    return i1_new, i2_new

def horizontal_crossover(i1, i2):
    crossover_point = np.random.randint(i1[:,0].size)
    i1_new = np.vstack((i1[:crossover_point,:], i2[crossover_point:,:]))
    i2_new = np.vstack((i1[crossover_point:,:], i2[:crossover_point,:]))
    return i1_new, i2_new

def mutation(i):
    mutation_i = np.random.randint(0, i[:,0].size)
    mutation_j = np.random.randint(0, i[0,:].size)
    i[mutation_i,mutation_j] = not i[mutation_i,mutation_j]
    return i

def tournament(population, fitness, n_competitiors = 3):
    competitors_index = np.random.choice(n_pop, n_competitiors, replace = False)
    competitors_fitness = fitness[competitors_index] 
    winner_index = competitors_index[np.argmax(competitors_fitness)]
    i_new = population[winner_index, :,:]
    fitness_new = fitness[winner_index]
    return i_new, fitness_new

def evaluate(i):
    cn_degrees = np.sum(i, 0)
    vn_degrees = np.sum(i, 1)

    condition1 = np.all(cn_degrees >= 2) and np.all(vn_degrees >= 2)
    condition2 = np.all(cn_degrees <= 25)

    if condition1 and condition2:
        rho_node = np.bincount(vn_degrees)
        rho_edge = np.arange(1,len(rho_node)+1) * rho_node / np.sum(np.arange(1,len(rho_node)+1) * rho_node)
        lam_node = np.bincount(cn_degrees)
        rho_edge = np.arange(1,len(lam_node)+1) * lam_node / np.sum(np.arange(1,len(lam_node)+1) * lam_node)

        fitness = np.max(vn_degrees)
        return fitness
    else:
        return -1

#%%
fname = None

if fname is None:
    population = init_population(n_pop, n_cn, n_vn)
    fitness = np.full(n_pop, -np.inf)
    n_pop = 50
    n_generations = 100
    n_vn = 73
    n_cn = 9
    p_vertical = 0.5
    p_horizontal = 0.5 
    p_mutation = 0.01

    fname = "test_ge.npz"
else:
    data = np.load(fname)
    
    population = init_population(n_pop, n_cn, n_vn)
    fitness = np.full(n_pop, -np.inf)
    
    n_pop = data["n_pop"]
    n_generations = data["n_generations"]
    n_vn = data["n_vn"]
    n_cn = data["n_cn"]
    
    p_vertical = data["p_vertical"]
    p_horizontal = data["p_horizontal"]
    p_mutation = data["p_mutation"]

    fname = fname

print(f"""
===================================================================
Optimization
===================================================================
""")

try:
    for i in range(n_generations):
        
        #Evaluate new individuals
        for j in range(n_pop):
            if fitness[j] == -np.inf:
                fitness[j] = evaluate(population[j,:,:])

        #Select with elitism
        population_new = np.zeros(population.shape, int)
        fitness_new = np.full(fitness.shape, -np.inf)
        population_new[0] = population[np.argmax(fitness)]
        fitness_new[0] = fitness[np.argmax(fitness)]
        for j in range(1,n_pop):
            population_new[j], fitness_new[j] = tournament(population, fitness) 
        population = population_new
        fitness = fitness_new

        #Crossover
        for j in range(1, n_pop, 2):
            if np.random.rand() < p_vertical:
                population[j], population[j+1] = vertical_crossover(population[j], population[j+1])
                fitness[j], fitness[j+1] = -np.inf, -np.inf
            if np.random.rand() < p_horizontal:
                population[j], population[j+1] = horizontal_crossover(population[j], population[j+1])
                fitness[j], fitness[j+1] = -np.inf, -np.inf

        #Mutation
        for j in range(1, n_pop):
            if np.random.rand():
                population[j] = mutation(population[j])
                fitness[j] = -np.inf

        status = f"Generation {i}/{n_generations} Best/median/mean fitness: {np.max(fitness)}"
        print(status, end='\r', flush=True)

finally:
    with open(fname, 'wb') as f:
        np.savez(f, population = population, 
                    fitness = fitness, 
                    
                    n_pop = np.array(n_pop),
                    n_generations = np.array(n_generations),
                    n_vn = np.array(n_vn),
                    n_cn = np.array(n_cn),
                    
                    p_vertical = np.array(p_vertical),
                    p_horizontal = np.array(p_horizontal),
                    p_mutation = np.array(p_mutation))
# %%
