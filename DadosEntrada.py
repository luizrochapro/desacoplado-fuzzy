#coding: utf-8

import numpy as np
import skfuzzy as fuzz
import matplotlib.pyplot as plt
from FuzzyMath import *

class DadosEntrada:

    def __init__(self,path):
        self.path = path
        self.sbase = 100.0
        self.barras = []
        self.ramos = []
        self.barras_fuzzy = []
        self.nb = 0  # número de barras
        self.tipo_barras = [] #tipos das barras
        self.vb = [] #módulo das tensões de barra
        self.ab = [] #ângulo das tensões de barra
        self.nr = 0 # Número de ramos
        self.bini = [] #barra início do ramo
        self.bfim = [] #barra fim do ramo
        self.unidis = [] #universo de discurso
        self.precisao = [] #precisão do universo de discurso
        self.pliq = [] # PG - PL  armazena o espaço de discurso da Pliq[0] e a MF[1]
        self.qliq = [] # QG - QL  armazena o espaço de discurso da Qliq[0] e a MF[1]
        self.e = [] #parte real da tensão de barra retangular
        self.f =[] #parte imaginaria da tensão de barra retangular

    def setPath(self,path):
        self.path = path


    def getPath(self,path):
        self.path = path


    def definir_tag(self,tag,s):
        if s != "99999" and tag != "none":
            s = tag

        switcher={
            "Sbase": "Sbase",
            "Barras": "Barras",
            "Ramos": "Ramos",
            "Barras_Fuzzy": "Barras_Fuzzy",
            "Universo_Discurso": "Universo_Discurso",
            "Precisao_Universo_Discurso": "Precisao_Universo_Discurso",
            "99999": "none"
        }
        return switcher.get(s,"none")

    def p2r(self, A, phi):
        ret_real =[]
        ret_imag =[]
        if len(A)==len(phi):        
            for k in range(0,len(A)):
                ret_real.append((A[k] * (np.cos(phi[k]) + np.sin(phi[k]) * 1j)).real)
                ret_imag.append((A[k] * (np.cos(phi[k]) + np.sin(phi[k]) * 1j)).imag)
            return ret_real,ret_imag
        else:
            print('Entrada inválida, vetores dimensões diferentes')
            return None

    def carregar_dados(self):
        #função que carrega os dados de entrada nas propriedades dos objetos
        filein = open(self.path,'r')
        tag = "none" # variavel tag armazena em que cartão estamos
        for line in filein:
            columns = line.split()    
            tag = self.definir_tag(tag,columns[0])
            if (tag == "Sbase") and (columns[0] != "Sbase") and (columns[0] != "%"):
                self.sbase = float(columns[0])
            if (tag == "Barras") and (columns[0] != "Barras" and (columns[0] != "%")):
                self.barras.append([float(i) for i in columns])
            if (tag == "Ramos") and (columns[0] != "Ramos" and (columns[0] != "%")):
                self.ramos.append([float(i) for i in columns])    
            if (tag == "Barras_Fuzzy") and (columns[0] != "Barras_Fuzzy" and (columns[0] != "%")):
                self.barras_fuzzy.append([float(i) for i in columns])
            if (tag == "Universo_Discurso") and (columns[0] != "Universo_Discurso" and (columns[0] != "%")):
                self.unidis = float(columns[0])
            if (tag == "Precisao_Universo_Discurso") and (columns[0] != "Precisao_Universo_Discurso" and (columns[0] != "%")):
                self.precisao = float(columns[0])
        self.barras = np.array(self.barras)
        self.ramos = np.array(self.ramos)
        self.barras_fuzzy = np.array(self.barras_fuzzy)
        self.nb = int(np.shape(self.barras)[0])
        self.tipo_barras = self.barras[:,1]
        self.vb = np.reshape(self.barras[:,2:5],(self.nb,3))
        self.ab = np.reshape(self.barras[:,5:8],(self.nb,3))
        self.nr = int(np.shape(self.ramos)[0])
        self.bini = [int(i) for i in self.ramos[:,0]]
        self.bfim = [int(i) for i in self.ramos[:,1]]
        self.unidis = np.linspace((-1*self.unidis), self.unidis, 10001)
        # converter tensão de polar para retangular
        for k in range(0,self.nb):
            x,y = self.p2r(self.vb[k],self.ab[k])   
            self.e.append(x)
            self.f.append(y)
        self.e = np.array(self.e)
        self.f = np.array(self.f)


    def calc_pliq(self):
        '''função que calcula vetor de potências ativas liquidas'''
        for k in range(0,self.nb):
            tri1 = FuzzyMath(np.array([self.barras_fuzzy[k,1],self.barras_fuzzy[k,2],self.barras_fuzzy[k,3]])) #PG
            tri2 = FuzzyMath(np.array([self.barras_fuzzy[k,7],self.barras_fuzzy[k,8],self.barras_fuzzy[k,9]])) #PL    
            self.pliq.append((tri1-tri2)*(1/self.sbase))
        return None

    def calc_qliq(self):
        '''função que calcula vetor de potências reativas liquidas'''
        for k in range(0,self.nb):
            tri1 = FuzzyMath(np.array([self.barras_fuzzy[k,4],self.barras_fuzzy[k,5],self.barras_fuzzy[k,6]])) #PG
            tri2 = FuzzyMath(np.array([self.barras_fuzzy[k,10],self.barras_fuzzy[k,11],self.barras_fuzzy[k,12]])) #PL
            self.qliq.append((tri1-tri2)*(1/self.sbase))
        return None
    
    def calc_pliq_newton(self):
        '''função calcula vetor de potências liquidas de potencia ativa'''
        pliq = []
        for k in range(0,self.nb):
            pliq.append(self.barras_fuzzy[k,2]-self.barras_fuzzy[k,8])
        return np.array(pliq)
    
    def calc_qliq_newton(self):
        qliq = []
        '''função calcula vetor de potências liquidas de potencia reativa'''
        for k in range(0,self.nb):
            qliq.append(self.barras_fuzzy[k,5]-self.barras_fuzzy[k,11])
        return np.array(qliq)
