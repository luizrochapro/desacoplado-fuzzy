#coding: utf-8

import numpy as np
import skfuzzy as fuzz
#from FuzzyMath import *

class UniversoDiscurso:

    def __init__(self,DF,G,B,e,f,nb,nr,bini,bfim,tb):
        self.DF = DF
        self.G = G
        self.B = B
        self.e = e
        self.f = f
        self.nb = nb
        self.nr = nr
        self.bini = bini
        self.bfim = bfim
        self.tb = tb
        self.DP = np.zeros((self.nb,3))
    
    def calc_dfmax_dxmax(self):
        '''Função que calcula os novos universos de discurso'''
        for k in range(0,self.nb): # transferindo objetos FuzzyMath para array bidimensional
            self.DP[k] = self.DF[k].f

        H = np.zeros((self.nb,self.nb))
        L = np.zeros((self.nb,self.nb))
        e,f = [],[]
        for t in range(0,3):
            for k in range(0,self.nb):
                e.append(self.e[k][t])
                f.append(self.f[k][t])
            # o cálculo do H
            for k in range(0,self.nb):
                if self.tb[k] != 2:
                    H[k,t] = 2*e[k] * self.G[k,k]
            
            for r in range(0,self.nr):
                if self.tb[k] != 2:
                    k = self.bini[r]
                    m = self.bfim[r]
                    k = k-1 #correção da dimensao para python
                    m = m-1 #correção da dimensao para python
                    H[k,m] = e[k] * self.G[k,m] + f[k] * self.B[k,m]    
            
            # o cálculo do L depende do tipo de barra
            for k in range(0, self.nb):
                if self.tb[k] == 1: # Barra PV
                    L[k,t] = 2*f[k]
                elif self.tb[k] == 0: # Barra PQ
                    L[k,t] = -2*f[k]*self.B[k,k]
                    
                    for k in range(0,self.nr):
                        k = self.bini[r]
                        m = self.bfim[r]
                        k = k-1 #correção da dimensao para python
                        m = m-1 #correção da dimensao para python
                        L[k,t] = L[k,t] + e[m]*self.G[k,m]-f[m]*self.B[k,m]
                else:
                    None
            for r in range(0, self.nr):
                if self.tb[k] == 1:
                    k = self.bini[r]
                    m = self.bfim[r]
                    k = k-1 #correção da dimensao para python
                    m = m-1 #correção da dimensao para python
                    L[k,t] = 0
                elif self.tb[k] == 0:
                    k = self.bini[r]
                    m = self.bfim[r]
                    k = k - 1
                    m = m - 1
                    L[k,t] = -e[k]*self.G[k,m]-f[k]*self.B[k,m] 
                else:
                    None         
        
        aux = np.absolute(self.DP[:])
        ind  = np.unravel_index(np.argmax(aux, axis=None), aux.shape)
        dFmax = aux[ind]
        aux = np.absolute(H[ind])
        dXmax = (1/aux)*dFmax
        return dFmax, dXmax