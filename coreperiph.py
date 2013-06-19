from numpy import zeros
from numpy import linalg
import random

def make_graph(n_c, n_p, prob_cc, prob_cp, prob_pp):
  size = n_c + n_p
  mat = zeros((size, size))
  for i in range(size):
    for j in range(size):
      if i == j:
        continue
      if i < n_c and j < n_c:
        if random.random() < prob_cc:
          mat[i,j] = mat[j,i] = 1
      elif i < n_c or j < n_c:
        if random.random() < prob_cp:
          mat[i,j] = mat[j,i] = 1
      else:
        if random.random() < prob_pp:
          mat[i,j] = mat[j,i] = 1

  return mat

def make_mod(G):
  size = len(G)
  edges = 0
  degrees = [0]*size
  B = G
  for i in range(size):
    for j in range(size):
      degrees[i] += G[i,j]
    edges += degrees[i]
  print edges
  for i in range(size):
    for j in range(size):
      B[i,j] = G[i,j] - degrees[i]*degrees[j]/float(edges)
  return B

G = make_graph(300, 300, .5, .005, .005)
#print Matrix(G).str()
#print G
B = make_mod(G)
#print B
#val = linalg.eig(B)
#print val
m = Matrix(B)
#print m.str()
eig = m.eigenvectors_right()
print eig[0]
