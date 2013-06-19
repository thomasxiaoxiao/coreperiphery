#By: Thomas
#the O(N+M) version of the bp_algorithm for non-degree corrected block model.
#implemented with parameter updates
from __future__ import division
import numpy as np
import random as rd
from time import clock, time

def BP_infer(N,q,gamma,omega,A,tmax):
    run = 0
    cav = np.zeros((N,N),dtype=object)
    phi = np.zeros(N,dtype=object)
    #initialize the matrix
    for i in range(N):
        for j in range(N):
            a = rd.random()
            cav[i][j] = [a,1-a]
    #initialize group assignment
    for i in range(N):
            buff =[]
            for r in range(q):
                PI =0
                for k in range(N):
                    if k!=i:
                        SIGMA = 0
                        for s in range(q):
                            SIGMA += cav[k][i][s]*(omega[r][s]**A.item(i,k))*np.exp(-omega[r][s])
                        PI +=np.log(SIGMA)
                buff.append(np.exp(np.log(gamma[r])+PI))
            phi[i]=[e/sum(buff) for e in buff]
    #iterations over cavity matrix
    while run < tmax:
        #calculate the log of the second product term(external field)
        PI2 = []
        for r in range(q):
            buff = 0
            for w in range(N):
                SIGMA = 0
                for s in range(q):
                    SIGMA += phi[w][s]*np.exp(-omega[r][s])
                buff += np.log(SIGMA)
            PI2.append(buff)
        #calculate the message
        for u in range(N):
            for v in range(N):
                if v!=u:
                    buff=[]
                    for r in range(q):
                        PI1 = 0
                        for w in range(N):
                            if w!=v and A.item(u,w)!=0:
                                SIGMA1 = 0
                                SIGMA2 = 0
                                for s in range(q):
                                    SIGMA1 += cav[w][u][s]*(omega[r][s]**A.item(w,u))*np.exp(-omega[r][s])
                                PI1 += np.log(SIGMA1)
                        buff.append(np.exp(np.log(gamma[r]) + PI1 + PI2[r]))
                    cav[u][v] = [e/sum(buff) for e in buff]
        #find the group assigment vector            
        for i in range(N):
            buff =[]
            for r in range(q):
                PI =0
                for k in range(N):
                    if k!=i:
                        SIGMA = 0
                        for s in range(q):
                            SIGMA += cav[k][i][s]*(omega[r][s]**A.item(i,k))*np.exp(-omega[r][s])
                        PI +=np.log(SIGMA)
                buff.append(np.exp(np.log(gamma[r])+PI))
            phi[i]=[e/sum(buff) for e in buff]
       
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
                            P = (omega[r][s]**A.item(u,v))*np.exp(-omega[r][s])*cav[u][v][r]*cav[v][u][s]
                            Norm += P
                            buff[r][s] = P
                    for r in range(q):
                        b[u][v][r]= [e/Norm for e in buff[r]]
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
                buff.append(w)
            omega[r] =buff
        run +=1
    ass =[e.index(max(e)) for e in phi]   
    return ass,phi,gamma,omega
