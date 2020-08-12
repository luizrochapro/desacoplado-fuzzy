#coding: utf-8

import numpy as np
import skfuzzy as fuzz

class FuzzyInfSystem:
    '''Nesta classe as funções de pertinencia são triangulares e as regras já estão previstas como lineares 
       de acordo com o artigo , IF LN THEN LN'''
    
    def __init__(self, dFmax, dXmax, precisaoF, precisaoX):
        self.dFmax = dFmax
        self.dXmax = dXmax
        self.precisaoF = precisaoF
        self.precisaoX = precisaoX
        #self.stepF = self.precisaoF
        #self.stepX = self.precisaoX
        '''
        num = 2
        while (self.stepF >= self.precisaoF):
            (self.uni_dis_F, self.stepF) = np.linspace(-dFmax, dFmax, num, endpoint=True, retstep=True, dtype=float)
            num+=1
        num = 2
        while (self.stepX >= self.precisaoX):
            (self.uni_dis_X, self.stepX)= np.linspace(-dXmax, dXmax, num, endpoint=True, retstep=True, dtype=float)
            num+=1
        '''
        self.uni_dis_F = np.arange(-self.dFmax, self.dFmax, self.precisaoF)
        self.uni_dis_X = np.arange(-self.dXmax, self.dXmax, self.precisaoX)

    def pert_funcs_df(self):
        '''Função cria funções de pertinência triangulares'''
        ln = fuzz.trimf(self.uni_dis_F, [-self.dFmax, -self.dFmax, -self.dFmax/3])
        mn = fuzz.trimf(self.uni_dis_F, [-self.dFmax, -self.dFmax/2, 0])
        sn = fuzz.trimf(self.uni_dis_F, [-self.dFmax/3, -self.dFmax/6, 0])
        #ssn = fuzz.trimf(self.uni_dis_F, [-self.dFmax/6, -self.dFmax/12, 0]) #acrescentei
        zr = fuzz.trimf(self.uni_dis_F, [-self.dFmax/12, 0, self.dFmax/12])
        #ssp = fuzz.trimf(self.uni_dis_F, [0, self.dFmax/12, self.dFmax/6]) #acrescentei
        sp = fuzz.trimf(self.uni_dis_F, [0, self.dFmax/6, self.dFmax/3])
        mp = fuzz.trimf(self.uni_dis_F, [0, self.dFmax/2, self.dFmax])
        lp = fuzz.trimf(self.uni_dis_F, [self.dFmax/3, self.dFmax, self.dFmax])
        pert_matrix_df = [ln,mn,sn,zr,sp,mp,lp]
        return(pert_matrix_df)

    def pert_funcs_dx(self):
        '''Função cria funções de pertinência triangulares'''
        ln = fuzz.trimf(self.uni_dis_X, [-self.dXmax, -self.dXmax, -self.dXmax/3])
        mn = fuzz.trimf(self.uni_dis_X, [-self.dXmax, -self.dXmax/2, 0])
        sn = fuzz.trimf(self.uni_dis_X, [-self.dXmax/3, -self.dXmax/6, 0])
        #ssn = fuzz.trimf(self.uni_dis_X, [-self.dXmax/6, -self.dXmax/12, 0]) #acrescentei
        zr = fuzz.trimf(self.uni_dis_X, [-self.dXmax/12, 0, self.dXmax/12])
        #ssp = fuzz.trimf(self.uni_dis_X, [0, self.dXmax/12, self.dXmax/6]) #acrescentei
        sp = fuzz.trimf(self.uni_dis_X, [0, self.dXmax/6, self.dXmax/3])
        mp = fuzz.trimf(self.uni_dis_X, [0, self.dXmax/2, self.dXmax])
        lp = fuzz.trimf(self.uni_dis_X, [self.dXmax/3, self.dXmax, self.dXmax])
        pert_matrix_dx = [ln,mn,sn,zr,sp,mp,lp]
        return(pert_matrix_dx)

    def activate_mfs(self, mf, matrix_mfs):
        '''Função que ativa as funções de pertinencia com a funcção de entrada'''
        act_mfs =[]
        for i in range(0, np.size(matrix_mfs, axis = 0)):
            act_mfs.append(fuzz.fuzzy_and(self.uni_dis_F, mf, self.uni_dis_F, matrix_mfs[i])[1])
        return act_mfs

    def calc_mfs_saida(self, act_mfs, matrix_mfs):
        '''Função que tranfere as mfs da entrada para saida'''
        mfs_saida = []
        for i in range(0, np.size(act_mfs, axis=0)):
            mfs_saida.append(np.fmin(np.max(act_mfs[i]), matrix_mfs[i]))     
        return mfs_saida          
    
    def agregar_mfs_saida(self, mfs_saida):
        '''Função que agrega as mfs de saida'''
        mfs_agregadas = mfs_saida[0]
        for i in range(1, np.size(mfs_saida, axis=0)-1):
            mfs_agregadas = np.fmax(mfs_agregadas, mfs_saida[i])
        return mfs_agregadas

    def calc_centroide(self, mfs_agregadas):
        '''função que defuzifica as funções agregadas'''
        return fuzz.defuzz(self.uni_dis_X, mfs_agregadas, 'centroid')