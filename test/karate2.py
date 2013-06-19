import networkx as nx
import numpy as np
from bp_function import BP_infer as bp
import matplotlib.pylab as plt

G = nx.read_gml('karate.gml')
A = np.array(nx.adjacency_matrix(G))

q = 2
n = len(G.nodes())
gamma =[0.5,0.5]
omega =[(0.25,0.04),(0.04,.25)]
tmax = 40
ass,phi,gamma,omega,messages = bp(n, q, gamma,omega,A,tmax)

print ass
#print phi
print gamma
print omega
for u in range(5):
  s = ""
  for v in range(5):
    s += "%d -> %d: " % (u,v)
    for r in range(q):
      s += "%f, " % messages[u][v][r]
    s += '\n'
  print s,

'''
norm  =ass[0]
for i in range(n):
    if ass[i]==norm:
        ass[i]=0
    else:
        ass[i]=1

'''

nx.draw(G,node_color =ass)
plt.show()
