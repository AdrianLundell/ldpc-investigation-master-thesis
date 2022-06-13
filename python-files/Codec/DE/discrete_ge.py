#%%
import numpy as np 
import time 

t0 = time.time()
n_cn = 9//3
n_vn = 73//3
n_pop = 50


def init_population(n_pop):
    population = np.zeros((n_pop, n_cn, n_vn))
    for i in range(n_pop):
        for j in range(n_vn):   
            k = np.random.choice(n_cn, 2, replace = False)
            population[i, k, j] = 1
    
    return population

def rho_coeffs(i):
    pass

def lam_coeffs(i):
    pass

def vertical_crossover(i1, i2):
    crossover_point = np.random.randint(i1[0,:].size)
    i1_new = np.vstack((i1[:,:crossover_point], i2[:,crossover_point:]))
    i2_new = np.vstack((i1[:,crossover_point:], i2[:,:crossover_point]))
    return i1_new, i2_new

def horizontal_crossover(i1, i2):
    crossover_point = np.random.randint(i1[:,0].size)
    i1_new = np.hstack((i1[:crossover_point,:], i2[crossover_point:,:]))
    i2_new = np.hstack((i1[crossover_point:,:], i2[:crossover_point,:]))
    return i1_new, i2_new

def mutation(i):
    mutation_index = np.random.randint(i.size)
    i[mutation_index] = not i[mutation_index]
    return i

def tournament(population, fitness, n_competitiors = 3):
    competitors_index = np.random.choice(n_pop, n_competitiors, replace = False)
    competitors_fitness = fitness[competitors_index] 
    winner_index = competitors_index[np.argmax(competitors_fitness)]
    i_new = population[winner_index, :,:]
    return i_new

def evaluate(i):
    cn_degrees = np.sum(i, 0)
    vn_degrees = np.sum(i, 1)

    condition1 = np.all(cn_degrees >= 2) and np.all(vn_degrees >= 2)
    condition2 = np.all(cn_degrees <= 25)

    if condition1 and condition2:
        rho_node = np.bincount(vn_degrees)
        rho_edge = np.arange(1,len(rho_node)) * rho_node / np.sum(np.arange(1,len(rho_node)) * rho_node)
        lam_node = np.bincount(cn_degrees)
        rho_edge = np.arange(1,len(lam_node)) * lam_node / np.sum(np.arange(1,len(lam_node)) * lam_node)

        fitness = 
        return fitness
    else:
        return -1

#%%
n_pop = 1
n_generations = 1
p_vertical = 0.5
p_horizontal = 0.5 
p_mutation = 0.01

population = init_population(n_pop)
fitness = np.zeros(n_pop)

try:
    for i in range(n_generations):
        #Evaluate
        for j in range(n_pop):
            fitness[j] = evaluate(population[j,:,:])

        #Select with elitism
        population_new = np.zeros(population.shape)
        population_new[0] = population[np.argmax(fitness)]
        for j in range(1,n_pop):
            population_new[j] = evaluate(population[j,:,:])   
        population = population_new

        #Crossover
        for j in range(1, n_pop, 2):
            if np.random.rand() < p_vertical:
                population[j], population[j+1] = vertical_crossover(population[j], population[j+1])
            if np.random.rand() < p_horizontal:
                population[j], population[j+1] = horizontal_crossover(population[j], population[j+1])

        #Mutation
        for j in range(1, n_pop):
            if np.random.rand():
                population[j] = mutation(population[j])


finally:
    with open('x_best.npy', 'wb') as f:
        np.save(f, x_best)