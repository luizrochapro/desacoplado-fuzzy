from newton import newton
import numpy as np
from Imprimir import *
from FuzzyMath import *

#filein = 'sis3_2.dat'
filein = 'sis3_3.dat'
#filein = 'sis13_2.dat'
#filein = 'sis13.dat'
#filein = 'sis6.dat'
#filein = 'sis14.dat'
#filein = 'sis14radial.dat'

def multDisp(d,v):
    aux = np.zeros((1,3))
    if d < 0:
        aux[0,0] = v[2] 
        aux[0,1] = v[1]
        aux[0,2] = v[0]
        return d*aux
    return d*v

sbase, v, ang, J, Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Y, d, npvpq, npq, gkm, bkm, bsh_k, bsh_km, akm = newton(filein)

G = np.real(Y)
B = np.imag(Y)


dZP = np.zeros((npvpq,3))
m = 0
for k in range(d.nb):
    if d.tipo_barras[k]!=2:
        a = FuzzyMath(d.pg[k])
        b = FuzzyMath(d.pl[k])
        c = a - b
        dZP[m] = c.asArray()
        for i in range(3):
            dZP[m,i] = (dZP[m,i] - Pi[k])/d.sbase
        m+=1

dZQ = np.zeros((npq,3))
m = 0
for k in range(d.nb):
    if d.tipo_barras[k]==0:
        a = FuzzyMath(d.qg[k])
        b = FuzzyMath(d.ql[k])
        c = a - b
        dZQ[m] = c.asArray()
        for i in range(3):
            dZQ[m,i] = (dZQ[m,i] - Qi[k])/d.sbase
        m+=1

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
    elif d.tipo_barras[k] == 2:
        d.ab[k,0] = d.ab[k,1]
        d.ab[k,2] = d.ab[k,1]

#atualizar módulo de tensão 
m = 0
for k in range(0, d.nb):
    if d.tipo_barras[k] == 0:
        for i in range(0, 3):
            d.vb[k,i] = d.vb[k,1] + dv[m,i]
        m+=1
    elif d.tipo_barras[k]!=0:
        d.vb[k,0] = d.vb[k,1]
        d.vb[k,2] = d.vb[k,1]

# derivadas parciais para fluxo ativo
'''
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPdVda[0,m] = -2 * d.vb[i,1] * G[i,k] + d.vb[k,1] * (G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dPdVda[1,m] = d.vb[i,1] * (G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]))      # derivada / Vk
    dPdVda[2,m] = d.vb[i,1] * d.vb[k,1] * (-G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))              # derivada / theta i
    dPdVda[3,m] = d.vb[i,1] * d.vb[k,1] * (G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) - B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))               # derivada / theta k
'''
dPkmdVda = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPkmdVda[0,m] = 2 * d.vb[i,1] * gkm[m] - akm[m] * d.vb[k,1] * (bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) + gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dPkmdVda[1,m] = -akm[m]*d.vb[i,1] * (gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]) + bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]))      # derivada / Vk
    dPkmdVda[2,m] = -akm[m] * d.vb[i,1] * d.vb[k,1] * (-1*gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) + bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]))              # derivada / theta i
    dPkmdVda[3,m] = -akm[m] * d.vb[i,1] * d.vb[k,1] * (gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) - bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]))               # derivada / theta k

dPmkdVda = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPmkdVda[0,m] = - akm[m] * d.vb[k,1] * (-bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) + gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dPmkdVda[1,m] = 2*gkm[m]*np.square(akm[m])*d.vb[k,1] - akm[m]*d.vb[i,1] * (gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]) - bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]))      # derivada / Vk
    dPmkdVda[2,m] = -akm[m] * d.vb[i,1] * d.vb[k,1] * (-1*gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) - bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]))              # derivada / theta i
    dPmkdVda[3,m] = -akm[m] * d.vb[i,1] * d.vb[k,1] * (gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) + bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]))               # derivada / theta k


# derivadas parciais para fluxo reativo
'''
#para dar o valor do trabalho de mestrado
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dQdVda[0,m] = d.vb[k,1] * (G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dQdVda[1,m] = - d.vb[i,1] * (G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vk        # derivada / Vk
    dQdVda[2,m] = d.vb[i,1] * d.vb[k,1] * (-B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) + G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))                # derivada / theta i
    dQdVda[3,m] = d.vb[i,1] * d.vb[k,1] * (B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) - G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))                 # derivada / theta k
'''
dQkmdVda = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dQkmdVda[0,m] = -2 * d.vb[i,1] * (bkm[m] + bsh_km[m]) - akm[m] * d.vb[k,1] * (gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) - bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dQkmdVda[1,m] = -akm[m] * d.vb[i,1] * (gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) - bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vk        # derivada / Vk
    dQkmdVda[2,m] = -akm[m] * d.vb[i,1] * d.vb[k,1] * (gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]) + bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]))                # derivada / theta i
    dQkmdVda[3,m] = -akm[m] * d.vb[i,1] * d.vb[k,1] * (-gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]) - bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]))                 # derivada / theta k

dQmkdVda = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dQmkdVda[0,m] = akm[m] * d.vb[k,1] * (gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) + bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vi
    dQmkdVda[1,m] = -2*np.square(akm[m]) * d.vb[k,1] * (bkm[m]+bsh_km[m]) + akm[m]* d.vb[i,1]*(gkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]) + bkm[m] * np.cos(d.ab[i,1]-d.ab[k,1])) # derivada / Vk        # derivada / Vk
    dQmkdVda[2,m] = akm[m] * d.vb[i,1] * d.vb[k,1] * (gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]) - bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]))                # derivada / theta i
    dQmkdVda[3,m] = akm[m] * d.vb[i,1] * d.vb[k,1] * (-gkm[m] * np.cos(d.ab[i,1]-d.ab[k,1]) + bkm[m] * np.sin(d.ab[i,1]-d.ab[k,1]))                 # derivada / theta k

# criar os vetores de imprecisao de tensão e angulo  deltaTheta e deltaAngulo fuzzy
dvf = np.zeros((d.nb,3))
for k in range(d.nb):
    dvf[k,0] = d.vb[k,0] - d.vb[k,1]
    dvf[k,1] = d.vb[k,1] - d.vb[k,1]
    dvf[k,2] = d.vb[k,2] - d.vb[k,1]

dangf = np.zeros((d.nb,3))
for k in range(d.nb):
    dangf[k,0] = d.ab[k,0] - d.ab[k,1]
    dangf[k,1] = d.ab[k,1] - d.ab[k,1]
    dangf[k,2] = d.ab[k,2] - d.ab[k,1]

#calculo das distribuições do fluxo de potencia ativa
dPik = np.zeros((d.nr,3))
aux = np.zeros((1,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    #dPik[m,:] = dPdVda[0,m]*dvf[i,:] + dPdVda[1,m]*dvf[k,:] + dPdVda[2,m]*dangf[i,:] + dPdVda[3,m]*dangf[k,:]
    dPik[m,:] = multDisp(dPkmdVda[0,m],dvf[i,:]) + multDisp(dPkmdVda[1,m],dvf[k,:]) + multDisp(dPkmdVda[2,m],dangf[i,:]) + multDisp(dPkmdVda[3,m],dangf[k,:])
    #testa se está com valores invertidos e corrigi
    #if dPik[m,2] < dPik[m,0]:
    #    aux[0,0] = dPik[m,2] 
    #    aux[0,1] = dPik[m,1] 
    #    aux[0,2] = dPik[m,0]
    #    dPik[m,:] = aux

#calculo das distribuições do fluxo de potencia ativa
dPki = np.zeros((d.nr,3))
aux = np.zeros((1,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    #dPik[m,:] = dPdVda[0,m]*dvf[i,:] + dPdVda[1,m]*dvf[k,:] + dPdVda[2,m]*dangf[i,:] + dPdVda[3,m]*dangf[k,:]
    dPki[m,:] = multDisp(dPmkdVda[0,m],dvf[i,:]) + multDisp(dPmkdVda[1,m],dvf[k,:]) + multDisp(dPmkdVda[2,m],dangf[i,:]) + multDisp(dPmkdVda[3,m],dangf[k,:])


#calculo das distribuições do fluxo de potencia reativa
dQik = np.zeros((d.nr,3))
aux = np.zeros((1,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    #dQik[m,:] = dQdVda[0,m]*dvf[i,:] + dQdVda[1,m]*dvf[k,:] + dQdVda[2,m]*dangf[i,:] + dQdVda[3,m]*dangf[k,:]
    dQik[m,:] = multDisp(dQkmdVda[0,m],dvf[i,:]) + multDisp(dQkmdVda[1,m],dvf[k,:]) + multDisp(dQkmdVda[2,m],dangf[i,:]) + multDisp(dQkmdVda[3,m],dangf[k,:])
    #testa se está com valores invertidos e corrigi
    #if dQik[m,2] < dQik[m,0]:
    #    aux[0,0] = dQik[m,2] 
    #    aux[0,1] = dQik[m,1] 
    #    aux[0,2] = dQik[m,0]
    #    dQik[m,:] = aux  

dQki = np.zeros((d.nr,3))
aux = np.zeros((1,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    #dQik[m,:] = dQdVda[0,m]*dvf[i,:] + dQdVda[1,m]*dvf[k,:] + dQdVda[2,m]*dangf[i,:] + dQdVda[3,m]*dangf[k,:]
    dQki[m,:] = multDisp(dQmkdVda[0,m],dvf[i,:]) + multDisp(dQmkdVda[1,m],dvf[k,:]) + multDisp(dQmkdVda[2,m],dangf[i,:]) + multDisp(dQmkdVda[3,m],dangf[k,:])

#calculo do fluxo fuzzy nos ramos
Pik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Pik[m,:] = Pij[i,k] + dPik[m]*sbase   

Pki = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Pki[m,:] = Pij[k,i] + dPki[m]*sbase   

#calculo do fluxo fuzzy nos ramos
Qik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Qik[m,:] = Qij[i,k] + dQik[m]*sbase 

Qki = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Qki[m,:] = Qij[k,i] + dQki[m]*sbase 

PerdasPik = np.zeros((d.nr,3))
for m in range(d.nr):
    PerdasPik[m,:] = Pik[m] + Pki[m]   


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
    Ep[m,:] = (dPkmdVda[0,m] * JinvExpV[i,:]) + (dPkmdVda[1,m] * JinvExpV[k,:]) + (dPkmdVda[2,m] * JinvExpT[i,:]) + (dPkmdVda[3,m] * JinvExpT[k,:])

#calculo matriz sensibilidade para potência reativa
Eq = np.zeros((d.nr,npvpq+npq))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Eq[m,:] = (dQkmdVda[0,m] * JinvExpV[i,:]) + (dQkmdVda[1,m] * JinvExpV[k,:]) + (dQkmdVda[2,m] * JinvExpT[i,:]) + (dQkmdVda[3,m] * JinvExpT[k,:])

'''
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
    Pik[m,:] = Pij[i,k] + dPik[m]*sbase   

#calculo do fluxo fuzzy nos ramos
Qik = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    Qik[m,:] = Qij[i,k] + dQik[m]*sbase      

'''
#calculo das derivadas parciais da perdas 
dPerdasik = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPerdasik[0,m] = 2*d.vb[i][1]*gkm[m] -2*d.vb[k][1]*gkm[m]*np.square(akm[m])*np.cos(d.ab[i][1]-d.ab[k][1])
    dPerdasik[1,m] = 2*d.vb[k][1]*gkm[m]*np.square(akm[k]) -2*d.vb[i][1]*gkm[m]*np.square(akm[m])*np.cos(d.ab[i][1]-d.ab[k][1])
    dPerdasik[2,m] = 2*gkm[m]*np.square(akm[m])*d.vb[i][1]*d.vb[k][1]*np.sin(d.ab[i][1]-d.ab[k][1])
    dPerdasik[3,m] = -2*gkm[m]*np.square(akm[m])*d.vb[i][1]*d.vb[k][1]*np.sin(d.ab[i][1]-d.ab[k][1])


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
    for j in range(npvpq):
        aux = S[m,j] * dZP[j,:]
        #aux = S[m,j] * dPik[m]
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
PerdasPikN = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    #i -= 1
    #k -= 1
    PerdasPikN[m,:] = Lpij[m] + dPerPik[m]*sbase   

########################
#Pki = np.zeros((d.nr,3))

dPe = np.zeros((d.nr,3))
aux = np.zeros((1,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPe[m,:] = multDisp(dPerdasik[0,m],dvf[i,:]) + multDisp(dPerdasik[1,m],dvf[k,:]) + multDisp(dPerdasik[2,m],dangf[i,:]) + multDisp(dPerdasik[3,m],dangf[k,:])

#calculo das perdas ativas
PerdasPikM = np.zeros((d.nr,3))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    #i -= 1
    #k -= 1
    PerdasPikM[m,:] = Lpij[m] + dPe[m]*sbase   

#%%
#######################
#S = np.matmul(dPerdasik , Jinv) 

#imprimir resultados no arquivo de saída
printer = Imprimir('saídas/resultados_{0}'.format(filein), d)
printer.write_results(Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Pik, Pki, Qik, Qki, PerdasPik)
print('Perdas com Matriz sensibilidade S')
print(PerdasPikN)
print('Perdas com delta V e delta Teta')
print(PerdasPikM)
print('Fim processamento')