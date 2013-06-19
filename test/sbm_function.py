#generates sbm graph of size n according canonical 2*2 matrix C_mat 
from __future__ import division
import numpy as np
import random as rd
import networkx as nx

def sbm(c_mat,n):
    A = np.zeros((n,n),float)
    for i in range(int(n/2)):
        for j in range(int(n/2)):
            if j!=i:
                buff = rd.random()
                if buff < c_mat[0][0]/n:
                    A[i][j] = 1
        for j in range(int(n/2),n):
            if j!=i:
                buff = rd.random()
                if buff < c_mat[0][1]/n:
                    A[i][j] = 1
                    
    for i in range(int(n/2),n):
        for j in range(int(n/2)):
            if j!=i:
                buff = rd.random()
                if buff < c_mat[1][0]/n:
                    A[i][j] = 1
        for j in range(int(n/2),n):
            if j!=i:
                buff = rd.random()
                if buff < c_mat[1][1]/n:
                    A[i][j] = 1
    G=nx.from_numpy_matrix(A)       
    return G    

