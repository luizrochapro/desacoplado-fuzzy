from newton import newton
import numpy as np
from Imprimir import *
from FuzzyMath import *

#filein = 'sis3_2.dat'
filein = 'sis13_2.dat'
#filein = 'sis6.dat'
#filein = 'sis14.dat'

v, ang, J, Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Y, d, npvpq, npq = newton(filein)

G = np.real(Y)
B = np.imag(Y)


dZP = np.zeros((npvpq,3))
for k in range(npvpq):
    if d.tipo_barras[k]!=2:
        a = FuzzyMath(d.pg[k])
        b = FuzzyMath(d.pl[k])
        c = a - b
        dZP[k] = c.asArray()
        for i in range(3):
            dZP[k,i] = (dZP[k,i] - Pi[k])/d.sbase

dZQ = np.zeros((npq,3))
for k in range(npq):
    if d.tipo_barras[k]==0:
        a = FuzzyMath(d.qg[k])
        b = FuzzyMath(d.ql[k])
        c = a - b
        dZQ[k] = c.asArray()
        for i in range(3):
            dZQ[k,i] = (dZQ[k,i] - Qi[k])/d.sbase

#concatenar os vetores dZP e dZQ
dZPdZQ = np.concatenate((dZP, dZQ), axis=0)

#inverso da jacobiana
Jinv = np.linalg.inv(J)

#calculo dos incrementos
aux = np.zeros((1,3))
dvec = np.zeros((npvpq+npq,3))
for i in range(npvpq+npq): #dimensao da Jinv
    for j in range(npvpq+npq):
        for k in range(3):
            aux[0,k] = Jinv[i,j] * dZPdZQ[j,k]
            if Jinv[i,j] < 0:
                a = aux[0,0]
                c = aux[0,2]
                aux[0,0] = c
                aux[0,2] = a                
        dvec[i,:] = dvec[i,:] + aux

dang = dvec[:npvpq]
dv = dvec[npvpq:]


#atualizar angulos
m = 0
for k in range(0, d.nb):
    if d.tipo_barras[k] != 2:
        for i in range(0,3):
            d.ab[k,i] = d.ab[k,1] + dang[m,i]
        m+=1


#atualizar módulo de tensão 
m = 0
for k in range(0, d.nb):
    if d.tipo_barras[k] == 0:
        for i in range(0, 3):
            d.vb[k,i] = d.vb[k,1] + dv[m,i]
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

# criar matrizes Jinv expandidas para todas barras
JinvExpT = np.zeros((d.nb,npvpq+npq))
JinvExpV = np.zeros((d.nb,npvpq+npq))
m = 0
n = npvpq
for k in range(d.nb):
    if d.tipo_barras[k]==1 or d.tipo_barras[k]==0:
        JinvExpT[k,:] = Jinv[m,:]
        if d.tipo_barras[k]==0:
            JinvExpV[k,:] = Jinv[n,:]
            n += 1
        m += 1

#calculo matriz sensibilidade para potência ativa
Ep = np.zeros((d.nr,npvpq+npq))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Ep[m,:] = (dPdVda[0,m] * JinvExpV[i,:]) + (dPdVda[1,m] * JinvExpV[k,:]) + (dPdVda[2,m] * JinvExpT[i,:]) + (dPdVda[3,m] * JinvExpT[k,:])

#calculo matriz sensibilidade para potência reativa
Eq = np.zeros((d.nr,npvpq+npq))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Eq[m,:] = (dQdVda[0,m] * JinvExpV[i,:]) + (dQdVda[1,m] * JinvExpV[k,:]) + (dQdVda[2,m] * JinvExpT[i,:]) + (dQdVda[3,m] * JinvExpT[k,:])

#calculo das distribuições do fluxo de potencia ativa
dPik = np.zeros((d.nr,3))
for m in range(d.nr):
    for j in range(npvpq+npq):
        aux = Ep[m,j] * dZPdZQ[j,:]
        if Ep[m,j] < 0:
            a = aux[0]
            c = aux[2]
            aux[0] = c
            aux[2] = a 
        dPik[m,:] = dPik[m,:] + aux

#calculo das distribuições do fluxo de potencia reativa
dQik = np.zeros((d.nr,3))
for m in range(d.nr):
    for j in range(npvpq+npq):
        aux = Eq[m,j] * dZPdZQ[j,:]
        if Eq[m,j] < 0:
            a = aux[0]
            c = aux[2]
            aux[0] = c
            aux[2] = a 
        dQik[m,:] = dQik[m,:] + aux

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

#calculo matriz sensibilidade para perdas
S = np.zeros((d.nr,npvpq+npq))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    S[m,:] = (dPerdasik[0,m] * JinvExpV[i,:]) + (dPerdasik[1,m] * JinvExpV[k,:]) + (dPerdasik[2,m] * JinvExpT[i,:]) + (dPerdasik[3,m] * JinvExpT[k,:])


#calculo das distribuições das perdas ativas
dPerPik = np.zeros((d.nr,3))
for m in range(d.nr):
    for j in range(npvpq+npq):
        aux = S[m,j] * dZPdZQ[j,:]
        if S[m,j] < 0:
            a = aux[0]
            c = aux[2]
            aux[0] = c
            aux[2] = a 
        dPerPik[m,:] = dPerPik[m,:] + aux
'''
#calculo das distribuições das perdas ativas
dPerPik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPerPik[m,:] = (dPerdasik[0,m] * dvf[i]) + (dPerdasik[1,m] * dvf[k]) + (dPerdasik[2,m] * dangf[i]) + (dPerdasik[3,m] * dangf[k])
'''
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
printer = Imprimir('saídas/resultados_{0}'.format(filein), d)
printer.write_results(Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Pik, Qik, PerdasPik)

print('Fim processamento')