##############################################################################
#  Algoritmo Newton Raphson                                                  #
#  Eng. Luiz Gonzaga Rocha Junior                                            #
#  contato@luizrocha.eng.br                                                  #
##############################################################################
#coding: utf-8

#import skfuzzy as fuzz
from DadosEntrada import *
#from UniversoDiscurso import *
#from FuzzyMath import *
#from FuzzyInfSystem import *
from Log import *
import numpy as np

#Arquivo de Log
log = Log('saida.log')
log.open_file()

#Dados de entrada
filein = 'sis6.dat'

# Instancia objeto dados
d = DadosEntrada('entradas/{0}'.format(filein))

#carrega arquivo de entrada para o objeto dados
d.carregar_dados()

#calcula as potências liquidas das barras em pu
p_esp = d.calc_pliq_newton()/d.sbase
q_esp = d.calc_qliq_newton()/d.sbase

# normalização em pu dos dados de shunt
bsh_k = d.barras[:,8]/d.sbase

#tipos de barras
npq = len(d.tipo_barras[d.tipo_barras==0])
npvpq = len(d.tipo_barras[d.tipo_barras==0]) + len(d.tipo_barras[d.tipo_barras==1])

# normalização em pu dos dados de ramos
### Cálculo das condutâncias e susceptâncias primitivas de ramos
gkm = 100 * d.ramos[:,2] / (np.square(d.ramos[:,2]) + np.square(d.ramos[:,3]))
bkm = -100 * d.ramos[:,3] / (np.square(d.ramos[:,2]) + np.square(d.ramos[:,3]))
xkm = d.ramos[:,3] / 100  # usado na formação da B' do desacoplado XB
bsh_km = d.ramos[:,4] / (2 * d.sbase) #Susceptância shunt de ramo
akm = d.ramos[:,5] #Tap dos transformadores defasadores (se LT, tap deve ser 1,0)

# ---------- Montagem das matrizes condutância G e susceptância B ---------
G = np.zeros((d.nb,d.nb))  
B = np.zeros((d.nb,d.nb))

for k in range(0,d.nb):
    B[k,k] = bsh_k[k]

for r in range(0,d.nr):
    k = d.bini[r] 
    m = d.bfim[r] 
    k = k-1 #correção da dimensao para python
    m = m-1 #correção da dimensao para python
    G[k,k] = G[k,k] + np.square(akm[r])*gkm[r]     
    G[k,m] = - akm[r]*gkm[r] 
    G[m,k] = - akm[r]*gkm[r] 
    G[m,m] = G[m,m] + gkm[r] 
    B[k,k] = B[k,k] + bsh_km[r] + np.square(akm[r])*bkm[r] 
    B[k,m] = - akm[r]*bkm[r] 
    B[m,k] = - akm[r]*bkm[r] 
    B[m,m] = B[m,m] + bsh_km[r] + bkm[r]


# ------- PROCESSO ITERATIVO DO MÉTODO NEWTON RAPSHON --------------------
## inicializações
iter = 0
tol = 1 # tolerância
'''
p = 0
q = 0
Kp = 1
Kq = 1
'''
#inicio do processo iterativo
while (iter < 100 and tol > 1e-5):
    
    # Calcular P and Q
    P = np.zeros((d.nb,1))
    Q = np.zeros((d.nb,1))
    for i in range(0, d.nb):
        for k in range(0, d.nb):
            P[i] = P[i] + d.vb[i][1] * d.vb[k][1] * (G[i,k]*np.cos(d.ab[i][1]-d.ab[k][1]) + B[i,k]*np.sin(d.ab[i][1]-d.ab[k][1]))
            Q[i] = Q[i] + d.vb[i][1] * d.vb[k][1] * (G[i,k]*np.sin(d.ab[i][1]-d.ab[k][1]) - B[i,k]*np.cos(d.ab[i][1]-d.ab[k][1]))
   
    #calculo dos mismatch
    dpa = p_esp.reshape(-1,1) - P
    dqa = q_esp.reshape(-1,1) - Q
    
    k = 0
    dQ = np.zeros((npq,1))
    for i in range(0,d.nb):
        if d.tipo_barras[i] == 0:
            dQ[k,0] = dqa[i]
            k += 1
    
    dP = dpa[1:d.nb]
    mismatch = np.concatenate((dP,dQ), axis=0) # Mismatch Vector 

    # cálculo da Jacobiana
    J = np.zeros((npq+npvpq,npq+npvpq))
    # H 
    H = np.zeros((npvpq,npvpq))
    for i in range(0, npvpq):
        m = i + 1
        for k in range(0,npvpq):
            n = k + 1 
            if n == m:
                for n in range(0,d.nb):
                    H[i,k] = H[i,k] + d.vb[m][1] * d.vb[n][1]*(-G[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]) + B[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]))
                H[i,k] = H[i,k] - np.square(d.vb[m][1])*B[m,m]
            else:
                H[i,k] = d.vb[m][1]* d.vb[n][1] * (G[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]) - B[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]))

    # N 
    N = np.zeros((npvpq,npq))
    for i in range(0, npvpq):
        m = i + 1
        for k in range(0, npq):
            n = np.where(d.tipo_barras==0)[0][k] #pega só as barras que são PQ
            if n == m:
                for n in range(0,d.nb):
                    N[i,k] = N[i,k] + d.vb[n][1]*(G[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]) + B[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]))
                N[i,k] = N[i,k] + 2*d.vb[m][1]*G[m,m]
            else:
                N[i,k] = d.vb[m][1] * (G[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]) + B[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]))

    # M 
    M = np.zeros((npq,npvpq))
    for i in range(0, npq):
        m = np.where(d.tipo_barras==0)[0][i] #pega só as barras que são PQ
        for k in range(0,npvpq):
            n = k + 1 
            if n == m:
                for n in range(0,d.nb):
                    M[i,k] = M[i,k] + d.vb[m][1] * d.vb[n][1]*(G[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]) + B[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]))
                M[i,k] = M[i,k] - np.square(d.vb[m][1])*G[m,m]
            else:
                M[i,k] = d.vb[m][1]* d.vb[n][1] * (-G[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]) - B[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]))

    # L
    L = np.zeros((npq,npq))
    for i in range(0, npq):
        m = np.where(d.tipo_barras==0)[0][i] #pega só as barras que são PQ
        for k in range(0, npq):
            n = np.where(d.tipo_barras==0)[0][k] #pega só as barras que são PQ
            if n == m:
                for n in range(0,d.nb):
                    L[i,k] = L[i,k] + d.vb[n][1]*(G[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]) - B[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]))
                L[i,k] = L[i,k] - 2*d.vb[m][1]*B[m,m]
            else:
                L[i,k] = d.vb[m][1] * (G[m,n]*np.sin(d.ab[m][1]-d.ab[n][1]) - B[m,n]*np.cos(d.ab[m][1]-d.ab[n][1]))


    J = np.concatenate((np.concatenate((H,M), axis=0),np.concatenate((N,L), axis=0)), axis=1) # Mismatch Vector 

    dx = np.matmul(np.linalg.inv(J), mismatch) #vetor de correção
    dtetha = dx[0:npvpq]
    dv = dx[npvpq:]

    #atualiza vetor de estado
    i = 0
    for k in range(0,d.nb):
        if k in np.where(d.tipo_barras!=2)[0]:
            d.ab[k][1] = d.ab[k][1] + dtetha[i]
            i += 1
    
    i = 0
    for k in range(0,d.nb):
        if k in np.where(d.tipo_barras==0)[0]:
            d.vb[k][1] = d.vb[k][1] + dv[i]   
            i+=1

    iter+=1
    tol = np.max(np.absolute(mismatch))

print('Nº iterações final = {0}'.format(iter))
print(d.vb)
print(d.ab)