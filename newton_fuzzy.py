from newton import newton
import numpy as np
from Imprimir import *

#filein = 'sis3_2.dat'
#filein = 'sis13.dat'
filein = 'sis6.dat'
#filein = 'sis14.dat'

v, ang, J, Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Y, d, npvpq, npq = newton(filein)

G = np.real(Y)
B = np.imag(Y)
'''
# cálculo do delta z
dZPg = np.zeros((d.nb,3))
for k in range(0,d.nb):
    for i in range(0,3):
        dZPg[k,i] = (d.pg[k,i] - d.pg[k,1])/d.sbase

dZPl = np.zeros((d.nb,3))
for k in range(0,d.nb):
    for i in range(0,3):
        dZPl[k,i] = (d.pl[k,i] - d.pl[k,1])/d.sbase

dZQg = np.zeros((d.nb,3))
for k in range(0,d.nb):
    for i in range(0,3):
        dZQg[k,i] = (d.qg[k,i] - d.qg[k,1])/d.sbase

dZQl = np.zeros((d.nb,3))
for k in range(0,d.nb):
    for i in range(0,3):
        dZQl[k,i] = (d.ql[k,i] - d.ql[k,1])/d.sbase
'''
dZP = np.zeros((npq,3))
for k in range(0,npq):
    for i in range(0,3):
        if d.tipo_barras[k]==0:
            dZP[k,i] = ((d.pg[k,1] - d.pl[k,1]) - (d.pg[k,i] - d.pl[k,i]))/d.sbase

dZQ = np.zeros((npvpq,3))
for k in range(0,npvpq):
    for i in range(0,3):
        if d.tipo_barras[k]!=2:
            dZQ[k,i] = ((d.qg[k,1] - d.ql[k,1]) - (d.qg[k,i] - d.ql[k,i]))/d.sbase

#concatenar os vetores dZP e dZQ
dZPdZQ = np.concatenate((dZP, dZQ), axis=0)

#inverso da jacobiana
Jinv = np.linalg.inv(J)

#calculo dos incrementos
dvec = np.matmul(Jinv,dZPdZQ)
dv = dvec[:npq]
dang = dvec[npq:]


#atualizar módulo de tensão 
m = 0
for k in range(0, d.nb):
    if d.tipo_barras[k] == 0:
        for i in range(0, 3):
            d.vb[k,i] = d.vb[k,1] + dv[m,i]
        m+=1

#atualizar angulos
m = 0
for k in range(0, d.nb):
    if d.tipo_barras[k] != 2:
        for i in range(0,3):
            d.ab[k,i] = d.ab[k,1] + dang[m,i]
        m+=1


# derivadas parciais para fluxo ativo
dPdVda = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPdVda[0,m] = -2 * d.vb[i,1] * G[i,k] + d.vb[k,1] * (G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dPdVda[1,m] = d.vb[i,1] * G[i,k] + d.vb[k,1] * (G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]))      # derivada / Vk
    dPdVda[2,m] = d.vb[i,1] * d.vb[k,1] * (-G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))              # derivada / theta i
    dPdVda[3,m] = d.vb[i,1] * d.vb[k,1] * (G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) - B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))               # derivada / theta k

# derivadas parciais para fluxo reativo
dQdVda = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dQdVda[0,m] = -2 * d.vb[i,1] * B[i,k] + d.vb[k,1] * (B[i,k] * np.cos(d.ab[i,1] - d.ab[k,1]) + G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dQdVda[1,m] = d.vb[i,1] * G[i,k] + d.vb[k,1] * (G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]))        # derivada / Vk
    dQdVda[2,m] = d.vb[i,1] * d.vb[k,1] * (-B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) + G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))                # derivada / theta i
    dQdVda[3,m] = d.vb[i,1] * d.vb[k,1] * (B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) - G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))                 # derivada / theta k


# criar os vetores de imprecisao de tensão e angulo  deltaTheta e deltaAngulo fuzzy
dvf = np.zeros((d.nb,3))
for k in range(d.nb):
    dvf[k,0] = d.vb[k,0]
    dvf[k,1] = 0
    dvf[k,2] = d.vb[k,2]

dangf = np.zeros((d.nb,3))
for k in range(d.nb):
    dangf[k,0] = d.ab[k,0]
    dangf[k,1] = 0
    dangf[k,2] = d.ab[k,2]    

#calculo das distribuições do fluxo de potencia ativa
dPik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPik[m,:] = (dPdVda[0,m] * dvf[i]) + (dPdVda[1,m] * dvf[k]) + (dPdVda[2,m] * dangf[i]) + (dPdVda[3,m] * dangf[k])

#calculo das distribuições do fluxo de potencia reativa
dQik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dQik[m,:] = (dQdVda[0,m] * dvf[i]) + (dQdVda[1,m] * dvf[k]) + (dQdVda[2,m] * dangf[i]) + (dQdVda[3,m] * dangf[k])

#calculo do fluxo fuzzy nos ramos
Pik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Pik[m,:] = Pij[i,k] + dPik[m]   

#calculo do fluxo fuzzy nos ramos
Qik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Qik[m,:] = Qij[i,k] + dQik[m]      


#calculo das derivadas parciais da perdas 
dPerdasik = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPerdasik[0,m] = 2*d.vb[i][1]*d.vb[k][1]*G[i,k]*np.sin(d.ab[i][1]-d.ab[k][1])
    dPerdasik[1,m] = -2*d.vb[i][1]*d.vb[k][1]*G[i,k]*np.sin(d.ab[i][1]-d.ab[k][1])
    dPerdasik[2,m] = 2*G[i,k]*d.vb[i][1]-2*d.vb[k][1]*G[i,k]*np.cos(d.ab[i][1]-d.ab[k][1])
    dPerdasik[3,m] = 2*G[i,k]*d.vb[k][1]-2*d.vb[i][1]*G[i,k]*np.cos(d.ab[i][1]-d.ab[k][1])

#calculo das distribuições das perdas ativas
dPerPik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPerPik[m,:] = (dPerdasik[0,m] * dvf[i]) + (dPerdasik[1,m] * dvf[k]) + (dPerdasik[2,m] * dangf[i]) + (dPerdasik[3,m] * dangf[k])

#calculo das perdas ativas
PerdasPik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    #i -= 1
    #k -= 1
    PerdasPik[m,:] = Lpij[m] + dPerPik[m]   


#S = np.matmul(dPerdasik , Jinv) 

#imprimir resultados no arquivo de saída
printer = Imprimir('saídas/resultados.lst', d)
printer.write_results(Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Pik, Qik, PerdasPik)

print('Fim processamento')