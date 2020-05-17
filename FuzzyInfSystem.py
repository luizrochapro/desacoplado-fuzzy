#coding: utf-8

import numpy as np
import skfuzzy as fuzz

class FuzzyInfSystem:
    
    def __init__(self, dFmax, dXmax, precisaoF, precisaoX):
        self.dFmax = dFmax
        self.dXmax = dXmax
        self.precisaoF = precisaoF
        self.precisaoX = precisaoX
        self.stepF = self.precisaoF
        self.stepX = self.precisaoX
        num = 2
        while (self.stepF >= self.precisaoF):
            (self.uni_dis_F, self.stepF) = np.linspace(-dFmax, dFmax, num, endpoint=True, retstep=True, dtype=float)
            num+=1
        num = 2
        while (self.stepX >= self.precisaoX):
            (self.uni_dis_X, self.stepX)= np.linspace(-dXmax, dXmax, num, endpoint=True, retstep=True, dtype=float)
            num+=1
    
    def pert_funcs_df(self):
        '''Função cria funções de pertinência triangulares'''
        ln = fuzz.trimf(self.uni_dis_F, [-self.dFmax, -self.dFmax, -self.dFmax/3])
        mn = fuzz.trimf(self.uni_dis_F, [-self.dFmax, -self.dFmax/2, 0])
        sn = fuzz.trimf(self.uni_dis_F, [-self.dFmax/3, -self.dFmax/6, 0])
        zr = fuzz.trimf(self.uni_dis_F, [-self.dFmax/12, 0, self.dFmax/12])
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
        zr = fuzz.trimf(self.uni_dis_X, [-self.dXmax/12, 0, self.dXmax/12])
        sp = fuzz.trimf(self.uni_dis_X, [0, self.dXmax/6, self.dXmax/3])
        mp = fuzz.trimf(self.uni_dis_X, [0, self.dXmax/2, self.dXmax])
        lp = fuzz.trimf(self.uni_dis_X, [self.dXmax/3, self.dXmax, self.dXmax])
        pert_matrix_dx = [ln,mn,sn,zr,sp,mp,lp]
        return(pert_matrix_dx)