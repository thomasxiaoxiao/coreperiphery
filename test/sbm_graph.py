from __future__ import division
import numpy as np
import random as rd
import networkx as nx
import matplotlib.pylab as plt
from bp_function2 import BP_infer as bp2
from sbm_function import sbm

n = 50
Cmat = [(10,3),(3,10)]
G=sbm(Cmat,n)

q =2
A = nx.adjacency_matrix(G)
gamma = [0.5,0.5]
omega1 =[(1,0.1),(0.1,1)]
tmax = 1

ass1,phi1,gamma1,omega1 =bp2(n,q,gamma,omega1,A,tmax)


nx.draw(G,node_color=ass1)
plt.show()

