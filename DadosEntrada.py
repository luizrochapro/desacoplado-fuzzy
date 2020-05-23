#coding: utf-8

import numpy as np
import skfuzzy as fuzz
import matplotlib.pyplot as plt

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


    def carregar_dados(self):
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
        self.vb = np.reshape(self.barras[:,2],(self.nb,1))
        self.ab = np.reshape(self.barras[:,3],(self.nb,1))
        self.nr = int(np.shape(self.ramos)[0])
        self.bini = [int(i) for i in self.ramos[:,0]]
        self.bfim = [int(i) for i in self.ramos[:,1]]
        #self.unidis = np.linspace((-1*self.unidis), self.unidis, int(2*self.unidis/self.precisao), endpoint=False, dtype=float)
        self.unidis = np.linspace((-1*self.unidis),self.unidis,10001)

    def calc_pliq(self):
        '''função que calcula vetor de potências ativas liquidas'''
        for k in range(0,self.nb):
            tri1 = fuzz.trimf(self.unidis,np.array([self.barras_fuzzy[k,1],self.barras_fuzzy[k,2],self.barras_fuzzy[k,3]])/self.sbase) #PG
            tri2 = fuzz.trimf(self.unidis,np.array([self.barras_fuzzy[k,7],self.barras_fuzzy[k,8],self.barras_fuzzy[k,9]])/self.sbase) #PL
            self.pliq.append(fuzz.dsw_sub(self.unidis,tri1,self.unidis,tri2,1000))
        return None

    def calc_qliq(self):
        '''função que calcula vetor de potências reativas liquidas'''
        for k in range(0,self.nb):
            tri1 = fuzz.trimf(self.unidis,np.array([self.barras_fuzzy[k,4],self.barras_fuzzy[k,5],self.barras_fuzzy[k,6]])/self.sbase) #PG
            tri2 = fuzz.trimf(self.unidis,np.array([self.barras_fuzzy[k,10],self.barras_fuzzy[k,11],self.barras_fuzzy[k,12]])/self.sbase) #PL
            self.qliq.append(fuzz.dsw_sub(self.unidis,tri1,self.unidis,tri2,1000))
        return None