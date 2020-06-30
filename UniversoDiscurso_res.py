#coding: utf-8

import numpy as np
import skfuzzy as fuzz
#from FuzzyMath import *

class UniversoDiscurso:

    def __init__(self,DF,G,B,vb,ab,nb,nr,bini,bfim,tb):
        self.DF = DF
        self.G = G
        self.B = B
        self.vb = vb
        self.ab = ab
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

        H = np.zeros((self.nb,1))
        L = np.zeros((self.nb,1))
        # o cálculo do H
        for k in range(0,self.nb):
            #if self.tb[k] != 2:
            H[k,0] = -1* self.B[k,k] * self.vb[k,1]**2
            
        for r in range(0,self.nr):
            #if self.tb[k] != 2:
            k = self.bini[r]
            m = self.bfim[r]
            k = k-1 #correção da dimensao para python
            m = m-1 #correção da dimensao para python
            dt = self.ab[k] - self.ab[m]
            H[k,0] = H[k,0] - self.vb[k,1]*self.vb[m,1]*(self.G[k,m]*np.sin(dt[1])-self.B[k,m]*np.cos(dt[1]))    
        
        # o cálculo do L depende do tipo de barra
        for k in range(0,self.nb):
            #if self.tb[k] != 2:
            L[k,0] = -1* self.B[k,k] * self.vb[k,1]
        
        for r in range(0,self.nr):
            #if self.tb[k] != 2:
            k = self.bini[r]
            m = self.bfim[r]
            k = k-1 #correção da dimensao para python
            m = m-1 #correção da dimensao para python
            dt = self.ab[k] - self.ab[m]
            L[k,0] = L[k,0] + self.vb[m,1]*(self.G[k,m]*np.sin(dt[1])-self.B[k,m]*np.cos(dt[1]))   
            
        aux = np.absolute(self.DP[:,1])
        ind  = np.unravel_index(np.argmax(aux, axis=None), aux.shape)
        dFmax = aux[ind]
        aux = np.absolute(H[ind])
        dXmax = (1/aux)*dFmax
        return dFmax, dXmax