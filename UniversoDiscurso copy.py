#coding: utf-8

import numpy as np
import skfuzzy as fuzz
from FuzzyMath import *

class UniversoDiscurso:

    def __init__(self,DF,G,B,vb,ab,nb,nr,bini,bfim,tb,ar,iter):
        '''
        if ar == 'ativo':
            self.DP = DF
        elif ar=='reativo':
            self.DQ = DF
        '''    
        
        if ar == 'ativo' or ar == 'reativo':
            self.DP = DF
            self.DQ = DF

        self.G = G
        self.B = B
        self.vb = vb
        self.ab = ab
        self.nb = nb
        self.nr = nr
        self.bini = bini
        self.bfim = bfim
        self.tb = tb
        self.ar = ar
        self.iter = iter

    def calc_dfmax_dxmax(self):
        '''Função que calcula os novos universos de discurso'''
        #for k in range(0,self.nb): # transferindo objetos FuzzyMath para array bidimensional
        #    self.DP[k] = self.DF[k].f
        H = []
        L = []
        for k in range(self.nb):
            H.append((FuzzyMath([0,0,0])))
            L.append((FuzzyMath([0,0,0])))

        # o cálculo do H
        for k in range(self.nb):
            H[k] = (self.vb[k]**2) * (-1) * float(self.B[k,k])
        
        for r in range(self.nr):
            k = self.bini[r]
            m = self.bfim[r]
            k = k-1 #correção da dimensao para python
            m = m-1 #correção da dimensao para python
            dt = self.ab[k] - self.ab[m]
            dtm = self.ab[m] - self.ab[k]
            H[k] = H[k] - self.vb[k]*self.vb[m]*(dt.sen()*float(self.G[k,m])-dt.cos()*float(self.B[k,m]))    
            H[m] = H[m] - self.vb[k]*self.vb[m]*(dtm.sen()*float(self.G[k,m])-dtm.cos()*float(self.B[k,m]))
        # o cálculo do L depende do tipo de barra
        for k in range(self.nb):
            #if self.tb[k] != 2:
            L[k] = self.vb[k] * (-1) * float(self.B[k,k])
        
        for r in range(0,self.nr):
            #if self.tb[k] != 2:
            k = self.bini[r]
            m = self.bfim[r]
            k = k-1 #correção da dimensao para python
            m = m-1 #correção da dimensao para python
            dt = self.ab[k] - self.ab[m]
            dtm = self.ab[m] - self.ab[k]
            L[k] = L[k] + self.vb[m]*(dt.sen()*float(self.G[k,m])-dt.cos()*float(self.B[k,m]))
            L[m] = L[m] + self.vb[m]*(dtm.sen()*float(self.G[k,m])-dtm.cos()*float(self.B[k,m]))   
        
        
        DF =[]
        Mat = []
        '''
        if self.ar == 'ativo':            
            for k in range(self.nb):
                DF.append([self.DP[k].f[0], self.DP[k].f[1], self.DP[k].f[2]])
                Mat.append(H[k])
            for k in range(self.nb):
                if self.tb[k] == 2:                
                    Mat[k] = FuzzyMath([0,0,0])

        elif self.ar == 'reativo':
            for k in range(self.nb):
                DF.append([self.DQ[k].f[0], self.DQ[k].f[1], self.DQ[k].f[2]])
                Mat.append(L[k])
            for k in range(self.nb):
                if self.tb[k] != 0:                
                    Mat[k] = FuzzyMath([0,0,0])
        '''
        if self.ar == 'ativo' or self.ar=='reativo':            
            for k in range(self.nb):
                DF.append([self.DP[k].f[0], self.DP[k].f[1], self.DP[k].f[2]])
                Mat.append(H[k])
            for k in range(self.nb):
                if self.tb[k] == 2:                
                    Mat[k] = FuzzyMath([0,0,0])
        '''
        aux = np.absolute(DF[:])
        ind  = np.unravel_index(np.argmax(aux, axis=None), aux.shape)
        dFmax = aux[ind]
        #aux = np.absolute(Mat[ind[0]].f[ind[1]])
        #aux = np.absolute(Mat[ind[0]].f[1])
        #dXmax = (1/(aux))*dFmax
        '''
        aux = []
        for k in range(self.nb):
            aux.append(np.absolute(DF[k][1]))

        #fator = 6 - (self.iter/1000)
        fator = 1
        print('Fator = {0}'.format(fator))
        ind  = np.argmax(aux)
        dFmax = fator*aux[ind]
        aux = np.absolute(Mat[ind].f[1])
        dXmax = (1/(aux))*dFmax
        
        return dFmax, dXmax