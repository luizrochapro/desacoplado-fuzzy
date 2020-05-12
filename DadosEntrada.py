#coding: utf-8

import numpy as np

class DadosEntrada:

    def __init__(self,path):
        self.path = path
        self.sbase = 100.0
        self.barras = []
        self.ramos = []
        self.cargas_fuzzy = []
        self.nb = 0  # número de barras
        self.tipo_barras = [] #tipos das barras
        self.vb = [] #módulo das tensões de barra
        self.ab = [] #ângulo das tensões de barra
        self.nr = 0 # Número de ramos
        self.bini = [] #barra início do ramo
        self.bfim = [] #barra fim do ramo

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
            "Cargas_Fuzzy": "Cargas_Fuzzy",
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
            if (tag == "Cargas_Fuzzy") and (columns[0] != "Cargas_Fuzzy" and (columns[0] != "%")):
                self.cargas_fuzzy.append([float(i) for i in columns])
        self.barras = np.array(self.barras)
        self.ramos = np.array(self.ramos)
        self.cargas_fuzzy = np.array(self.cargas_fuzzy)
        self.nb = int(np.shape(self.barras)[0])
        self.tipo_barras = self.barras[:,1]
        self.vb = np.reshape(self.barras[:,2],(self.nb,1))
        self.ab = np.reshape(self.barras[:,3],(self.nb,1))
        self.nr = int(np.shape(self.ramos)[0])
        self.bini = [int(i) for i in self.ramos[:,0]]
        self.bfim = [int(i) for i in self.ramos[:,1]]