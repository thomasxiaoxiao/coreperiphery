import networkx as nx
import numpy as np
from beliefprop import bp
import matplotlib.pylab as plt

G = nx.read_gml('karate.gml')
A = np.array(nx.adjacency_matrix(G))

q = 2
n = len(G.nodes())
gamma =[0.5,0.5]
omega =[(0.25,0.04),(0.04,.25)]
tmax = 60
ass,phi,gamma,omega,messages = bp(gamma,omega,A,tmax)

for u in range(5):
  s = ""
  for v in range(5):
    s += "%d -> %d: " % (u,v)
    for r in range(q):
      s += "%f, " % messages[u][v][r]
    s += '\n'
  print s,

nx.draw(G,node_color =ass)
plt.show()
