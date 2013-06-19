#By: Thomas
#the O(N^2) version of the bp_algorithm
#implemented with degree guessing
from __future__ import division
import numpy as np
import random as rd

def prob(adj, omega):
    return omega if adj == 0 else 1 - omega
    #return (omega**adj)*np.exp(-omega)

def BP_infer(N,q,gamma,omega,A,tmax):
    error = 10
    run = 0
    cav = np.zeros((N,N),dtype=object)
    phi = np.zeros(N,dtype=object)
    #initialize the matrix
    for i in range(N):
        for j in range(N):
            a = .49
            cav[i][j] = [a,1-a]
            
    d =np.zeros(N,int)
    for u in range(N):
        d[u] = A[u].sum()
        
    #iterations over cavity matrix
    while run < tmax:
        cav_temp = np.copy(cav)
        for u in range(N):
            for v in range(N):
                if v!=u:
                    buff=[]
                    for r in range(q):
                        PI = 0
                        for w in range(N):
                            if w!=u and w!=v:
                                SIGMA = 0
                                for s in range(q):
                                    SIGMA += cav[w][u][s]*prob(A.item(w,u),omega[r][s])
                                PI += np.log(SIGMA)
                        buff.append(np.exp(np.log(gamma[r]) + PI))
                    cav_temp[u][v] = [e/sum(buff) for e in buff]
        cav = np.copy(cav_temp)

        #find the group assigment vector            
        for i in range(N):
            buff =[]
            for r in range(q):
                PI =0
                for k in range(N):
                    if k!=i:
                        SIGMA = 0
                        for s in range(q):
                            SIGMA += cav[k][i][s]*prob(A.item(i,k),omega[r][s])
                        PI +=np.log(SIGMA)
                buff.append(np.exp(np.log(gamma[r])+PI))
            phi[i]=[e/sum(buff) for e in buff]
        
        #print phi
        #calculate the joint marginals
        b =np.zeros((N,N,2),dtype=object)
        for u in range(N):
            for v in range(N):
                if v!=u:
                    Norm = 0
                    buff=np.zeros((q,q),float)
                    for r in range(q):
                        P = 0
                        for s in range(q):
                            P = prob(A.item(u,v),omega[r][s])*cav[u][v][r]*cav[v][u][s]
                            Norm += P
                            buff[r][s] = P
                    for r in range(q):
                        b[u][v][r]= [e/Norm for e in buff[r]]
        #print b[0][1]
        #print b[32][33]
        #update the parameters
        gamma=[]
        for r in range(q):
            SIGMA = 0
            for i in range(N):
                SIGMA += phi[i][r]
            gamma.append(SIGMA)
        gamma = [e/N for e in gamma]
        

        omega = np.zeros((q,q),float)    
        for r in range(q):
            buff=[]
            for s in range(q):
                SIGMA = 0
                for u in range(N):
                    for v in range(N):
                        if v !=u and A.item(u,v)!=0:
                            SIGMA += b[u][v][r][s]
                w = SIGMA/(gamma[r]*gamma[s]*N**2)
                if w > 1:
                  w = 1
                assert 0 <= w <= 1
                buff.append(w)
            omega[r] =buff
        #print omega
        #print run
        run += 1
        
    #ass =[e.index(max(e)) for e in phi]     
    ass =[e[0] for e in phi]     
    return ass,phi,gamma,omega,cav
