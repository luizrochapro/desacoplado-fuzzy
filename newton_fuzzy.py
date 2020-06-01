from newton import newton
import numpy as np

v, ang, J, Pij, Qij, Pg, Qg, Lij, Y, d = newton('sis3_2.dat')

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
dZP = np.zeros((d.nb,3))
for k in range(0,d.nb):
    for i in range(0,3):
        dZP[k,i] = ((d.pg[k,1] - d.pl[k,1]) - (d.pg[k,i] - d.pl[k,i]))/d.sbase

dZQ = np.zeros((d.nb,3))
for k in range(0,d.nb):
    for i in range(0,3):
        dZQ[k,i] = ((d.qg[k,1] - d.ql[k,1]) - (d.qg[k,i] - d.ql[k,i]))/d.sbase


#calculo dos incrementos
dv = np.matmul(np.linalg.inv(J), dZP)
dang = np.matmul(np.linalg.inv(J), dZQ)

#atualizar módulo de tensão 
for k in range(0, d.nb):
    for i in range(0,3):
        d.vb[k,i] = d.vb[k,1] + dv[k,i]

#atualizar angulos
for k in range(0, d.nb):
    for i in range(0,3):
        d.ab[k,i] = d.ab[k,1] + dang[k,i]

# derivadas parciais para fluxo ativo
dPdVda = np.zeros((4,d.nr))
for m in range(d.nr):
    i = d.bini[m]
    k = d.bfim[m]
    i -= 1
    k -= 1
    dPdVda[0,m] = -2 * d.vb[i,1] * G[i,k] + d.vb[k,1] * (G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]))
    dPdVda[1,m] = d.vb[i,1] * G[i,k] + d.vb[k,1] * (G[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]))
    dPdVda[2,m] = d.vb[i,1] * d.vb[k,1] * (-G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) + B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))
    dPdVda[3,m] = d.vb[i,1] * d.vb[k,1] * (G[i,k] * np.sin(d.ab[i,1]-d.ab[k,1]) - B[i,k] * np.cos(d.ab[i,1]-d.ab[k,1]))



print(d.vb)
print(d.ab)
print(dPdVda)