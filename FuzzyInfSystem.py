#coding: utf-8

import numpy as np
import skfuzzy as fuzz

class FuzzyInfSystem:


    def __init__(self, undis_in, undis_out):
        self.undis_in = undis_in
        self.undis_out = undis_out
  
  # montar  as funções de pertinÊncias das regras

  # Generate fuzzy membership functions
    '''
        qual_tt = fuzz.trimf(x_qual,[0, 3, 5])
        qual_lo = fuzz.trimf(x_qual, [0, 0, 5])
        qual_md = fuzz.trimf(x_qual, [0, 5, 10])
        qual_hi = fuzz.trimf(x_qual, [5, 10, 10])
        serv_lo = fuzz.trimf(x_serv, [0, 0, 5])
        serv_md = fuzz.trimf(x_serv, [0, 5, 10])
        serv_hi = fuzz.trimf(x_serv, [5, 10, 10])
        tip_lo = fuzz.trimf(x_tip, [0, 0, 13])
        tip_md = fuzz.trimf(x_tip, [0, 13, 25])
        tip_hi = fuzz.trimf(x_tip, [13, 25, 25])
    '''

  # usar função "and" como no exemplo abaixo para obter a saida
  # uni, aux = fuzz.fuzzy_and(x_qual,qual_md,x_qual,qual_tt)
  # obter o max das funções do passo anterior
  # deffuzificar com o centroide