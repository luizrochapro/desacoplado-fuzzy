import numpy as np
from FuzzyMath import *

def calc_dpcalc(nb,nr,bini,bfim,G,B,v,ab,base):
    dpcalc = []
    dPd = np.zeros((2,nb))
    #calculo do Pk
    Pcalc = np.zeros((nb,1))
    for k in range(nb):
        Pcalc[k,0] = np.square(v[k][1])*G[k,k]
    for r in range(nr):
        k = bini[r]
        m = bfim[r]
        k -= 1
        m -= 1
        Pcalc[k,0] = Pcalc[k,0] + v[k][1]*v[m][1]*(G[k,m]*np.cos(ab[k,1]-ab[m,1])+B[k,m]*np.sin(ab[k,1]-ab[m,1]))

    for k in range(nb):
        dPd[0,k] = 2*v[k][1]*G[k,k]
        dPd[1,k] = 0

    for r in range(nr):
        k = bini[r]
        m = bfim[r]
        k -= 1
        m -= 1
        dPd[0,k] = dPd[0,k] + v[m][1] * (G[k,m] * np.cos(ab[k,1]-ab[m,1]) + B[k,m] * np.sin(ab[k,1]-ab[m,1]))
        dPd[1,k] = dPd[1,k] + v[k][1]*v[m][1]*(-1*G[k,m] * np.sin(ab[k,1]-ab[m,1]) + B[k,m] * np.cos(ab[k,1]-ab[m,1]))
    
    for k in range(nb):
        dpcalc.append(FuzzyMath((dPd[0,k]/base) * v[k] + (dPd[1,k]/base) * ab[k]))

    #atualizar dispersao da potencia calculada com a dispersao das derivadas
    for k in range(nb):
        a = dpcalc[k].f[0]-dpcalc[k].f[1]    
        c = dpcalc[k].f[2]-dpcalc[k].f[1]
        dpcalc[k].f[0] = (Pcalc[k]/base) + a
        dpcalc[k].f[1] = Pcalc[k]/base
        dpcalc[k].f[2] = (Pcalc[k]/base) + c 

    return dpcalc


def calc_dqcalc(nb,nr,bini,bfim,G,B,v,ab,base):
    dqcalc = []
    dQd = np.zeros((2,nb))
    #calculo do Qk
    Qcalc = np.zeros((nb,1))
    for k in range(nb):
        Qcalc[k,0] = -1*np.square(v[k][1])*B[k,k]
    for r in range(nr):
        k = bini[r]
        m = bfim[r]
        k -= 1
        m -= 1
        Qcalc[k,0] = Qcalc[k,0] + v[k][1]*v[m][1]*(G[k,m]*np.sin(ab[k,1]-ab[m,1])-B[k,m]*np.cos(ab[k,1]-ab[m,1]))

    for k in range(nb):
        dQd[0,k] = -2*v[k][1]*B[k,k]
        dQd[1,k] = 0

    for r in range(nr):
        k = bini[r]
        m = bfim[r]
        k -= 1
        m -= 1
        dQd[0,k] = dQd[0,k] + v[m][1] * (G[k,m] * np.sin(ab[k,1]-ab[m,1]) - B[k,m] * np.cos(ab[k,1]-ab[m,1]))
        dQd[1,k] = dQd[1,k] + v[k][1]*v[m][1]*(G[k,m] * np.cos(ab[k,1]-ab[m,1]) + B[k,m] * np.sin(ab[k,1]-ab[m,1]))
    
    for k in range(nb):
        dqcalc.append(FuzzyMath((dQd[0,k]/base) * v[k] + (dQd[1,k]/base) * ab[k]))

    #atualizar dispersao da potencia calculada com a dispersao das derivadas
    for k in range(nb):
        a = dqcalc[k].f[0]-dqcalc[k].f[1]    
        c = dqcalc[k].f[2]-dqcalc[k].f[1]
        dqcalc[k].f[0] = (Qcalc[k]/base) + a
        dqcalc[k].f[1] = Qcalc[k]/base
        dqcalc[k].f[2] = (Qcalc[k]/base) + c 

    return dqcalc