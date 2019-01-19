import igraph as ig
import numpy as np
import scipy
import scipy.sparse
from time import sleep
np.set_printoptions(linewidth=np.inf)


def tore( n ):
    g = ig.Graph()
    g.add_vertices(n ** 2)
    edges = []
    for i in range(n):
        for j in range(n):
            vertex = i*n + j
            edges.append((vertex, ((i + 1)%n ) *n  + j))
            edges.append((vertex, (i * n + (1 + j)%n )))
    g.add_edges(edges)
    return g


def unit_square( n ):
    g = ig.Graph()
    g.add_vertices(n ** 2)
    edges = []
    for i in range(n):
        for j in range(n):
            vertex = i*n + j
            if i != n - 1:
                edges.append((vertex, ((i + 1)%n ) *n  + j))
            if j != n - 1:
                edges.append((vertex, (i * n + (1 + j)%n )))
    g.add_edges(edges)
    return g


class Model( object ):


    def __init__( self, alpha, beta, mu, nu, T, dt, graph):
        self.alpha = alpha
        self.beta = beta
        self.nu = nu
        self.mu = mu
        self.graph = graph
        self.T = T
        self.iterations = 0
        self.dt = dt
        self.number =  len(graph.vs)
        self.laplacien = None
        rand = np.random.randn(self.number)
        rand *= 0.2/max(rand)
        self.suspected = np.ones(self.number) / 2 + rand
        self.infected = np.ones(self.number) / 2 - rand
        self.recovered = np.zeros(self.number)


    def update_laplacien( self ):
        neighbor_list = self.graph.neighborhood(order=1)
        data = []
        indices = []
        indptr = [0]
        for vertex, neighbors in enumerate(neighbor_list):
            neighbors.remove(vertex)
            n = len(neighbors)
            data = data + [-1] + [1. / n] * n
            indices = indices + [vertex] + neighbors
            indptr.append(indptr[-1] + n + 1)
        self.laplacien = scipy.sparse.csr_matrix((data, indices, indptr), shape=(self.number, self.number))*self.number


    def iterate(self):
        self.iterations += 1
        old_suspected, old_infected, old_recovered = self.suspected, self.infected, self.recovered
        n = old_suspected + old_infected + old_recovered
        self.suspected = (self.mu * ( n - old_suspected) - self.beta *  old_suspected * old_infected
                          + self.alpha * self.laplacien @ old_suspected) * self.dt + old_suspected
        self.infected = ( - (self.nu + self.mu) * old_infected + self.beta * old_suspected * old_infected
                          + self.alpha * self.laplacien @ old_infected) * self.dt + old_infected
        self.recovered = (self.nu * old_infected - self.mu * old_recovered
                          + self.alpha * self.laplacien @ old_recovered) * self.dt + old_recovered
        self.apply_BndCondition()


    def apply_BndCondition(self):
        pass

class UnitSquareModel( Model ):


    def __init__(self, alpha, beta, mu, nu, T, dt, size):
        self.size = size
        Model.__init__(self, alpha, beta, mu, nu, T, dt, graph = unit_square(size))


    def bndCondition(self, array):
        for j in range(self.size):
            array[j] = array[self.size + j]
            array[self.size * (self.size - 1) + j] = array[self.size * (self.size - 2) + j]
        for i in range(self.size):
            array[i * self.size] = array[i * self.size + 1]
            array[(i + 1) * self.size - 1] = array[(i + 1) * self.size - 2]


    def apply_BndCondition(self):
        self.bndCondition(self.suspected)
        self.bndCondition(self.infected)
        self.bndCondition(self.recovered)




def discretize_function(f):
    return np.array([[f(i/N, j/N) for i in range(N)] for j in range(N)]).flatten()

def f(x, y):
    return np.exp(-.5/sigma**2 * ((x-.9)**2 + (y-.9)**2))

sigma = 0.00001
N = 100
dt = 1/10000
model = UnitSquareModel(0.05, 5000, 1., .2, 10, dt, N)

model.suspected = discretize_function(f)
model.infected = discretize_function(f)

model.update_laplacien()




