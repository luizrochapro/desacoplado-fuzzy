#coding: utf-8

import numpy as np
import skfuzzy as fuzz

class UniversoDiscurso:

    def __init__(self, unidis,DF,G,B,v,a,nb,nr,bini,bfim):
        self.unidis = unidis
        self.DF = DF
        self.G = G
        self.B = B
        self.v = v
        self.a = a
        self.nb = nb
        self.nr = nr
        self.bini = bini
        self.bfim = bfim
    
    def calc_dfmax_dxmax(self):
        '''Função que calcula os novos universos de discurso'''
        DX = np.zeros[self.nb,3]
        DT = np.zeros[self.nb,3]
        
        for n in range(0,self.nb):
            DT[n,:] = fuzz.trimf(self.unidis, [self.DF[n,0],self.DF[n,1],self.DF[n,2]])

        for t in range(0,3):
            vb = self.v[:,t]
            ab = self.a[:,t]

            for k in range(0,self.nb):
                H[k,t] = -self.B*np.squared(vb[k])
            
            for r in range(0,self.nr):
                k = self.bini[r]
                m = self.bfim[r]
                dt = ab[k] - ab[m]
                H[k,t] = H[k] - vb[k]*vb[m]*(G[k,m]*np.sin(dt)-B[k,m]*np.cos(dt))
            
            for k in range(0,self.nb):
                L[k,t] = -B[k,k]*vb[k]

            for r in range(0,self.nr):
                k = self.bini[r]
                m = self.bfim[r]
                dt = ab[k] - ab[m]
                L[k,t] = L[k] + vb[m]*(G[k,m]*np.sin(dt)-B[k,m]*mp.cos(dt))              
        
            aux = np.absolute(DT[:])
            dFmax = np.amax(aux, axis = 0)
            ind = np.argmax(aux, axis = 0)
            aux = np.absolute(H[ind])
            dXmax = np.inverse(aux)*dFmax