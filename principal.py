##############################################################################
#  Algoritmo desacoplado rápido                                              #
#  Eng. Luiz Gonzaga Rocha Junior                                            #
#  contato@luizrocha.eng.br                                                  #
##############################################################################
#coding: utf-8

import skfuzzy as fuzz
from DadosEntrada import *
from UniversoDiscurso import *
from FuzzyMath import *
from FuzzyInfSystem import *
from Log import *
import matplotlib.pyplot as plt
#Arquivo de Log
log = Log('saida.log')
log.open_file()

#Dados de entrada
#filein = 'sis6.dat'
#filein = 'sis6_varang.dat'
#filein = 'sis3_1.dat'
#filein = 'sis6_2.dat'
#
#filein = 'sis2.dat'
#filein = 'sis3_2.dat' #fator = 1
#filein = 'sis6_vbcte.dat' #fator = 4.92 
#filein = 'sis13.dat'
#filein = 'sis14.dat'  #fator = 3.22 ou 3.24
filein ='sis14_5.dat'

# Instancia objeto dados
d = DadosEntrada('entradas/{0}'.format(filein))

#carrega arquivo de entrada para o objeto dados
d.carregar_dados()

#calcula as potências liquidas das barras em pu
d.calc_pliq()
d.calc_qliq()

#log.write_log_div()
#log.write_log(d.pliq)
#log.write_log(d.qliq)

# normalização em pu dos dados de shunt
bsh_k = d.barras[:,14]/d.sbase

# normalização em pu dos dados de ramos
### Cálculo das condutâncias e susceptâncias primitivas de ramos
gkm = 100 * d.ramos[:,2] / (np.square(d.ramos[:,2]) + np.square(d.ramos[:,3]))
bkm = -100 * d.ramos[:,3] / (np.square(d.ramos[:,2]) + np.square(d.ramos[:,3]))
xkm = d.ramos[:,3]/100 # usado na formação da B' do desacoplado XB
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


# ===================== CÁLCULO DO SUBSISTEMA 1 ===========================

## ------- MONTAGEM E INVERSA DAS MATRIZES DE SENSIBILIDADE B' e B" -------
B1L = np.zeros((d.nb,d.nb))
B2L = np.zeros((d.nb,d.nb))

for r in range(0, d.nr):
    k = d.bini[r] 
    m = d.bfim[r] 
    k = k-1 #correção da dimensao para python
    m = m-1 #correção da dimensao para python
    B1L[k,k] = B1L[k,k] + B[k,m]
    B1L[m,m] = B1L[m,m] + B[m,k]
    B1L[k,m] = - B[k,m]
    B1L[m,k] = - B[m,k]
    B2L[k,k] = B2L[k,k] + 1/xkm[r]
    B2L[m,m] = B2L[m,m] + 1/xkm[r]
    B2L[k,m] = -1/xkm[r]
    B2L[m,k] = -1/xkm[r]

# Adequação das matrizes B' e B'' para desconsiderar barras V-Teta e P-V
for k in range (0,d.nb):
    if d.tipo_barras[k] == 1:
        B2L[k,k] = 1.e+10

    if d.tipo_barras[k] == 2:
        B1L[k,k] = 1.e+10
        B2L[k,k] = 1.e+10

B1L = np.linalg.inv(B1L)
B2L = np.linalg.inv(B2L)

#calcular dispersão do módulo da tensão e ângulo
d.calc_disp(B1L,B2L)

#Converter em valores fuzzy
vb, ab = [], []
for k in range(d.nb):
    vb.append(FuzzyMath(d.vb[k]))
    ab.append(FuzzyMath(d.ab[k]))


# ------- PROCESSO ITERATIVO DO MÉTODO DESACOPLADO RÁPIDO VERSÃO XB -------
## inicializações
iter = 0
p = 0
q = 0
Kp = 1
Kq = 1
DPant = []
DQant = []
diffP = []
DPAcum =[]
DXPAcum = []
DQAcum = []
DXQAcum = []

errP = 1e-4
errQ = 1e-4
#errP = 0.003
#errQ = 0.003
presF = 0.01
presX = 0.0000001
maxiter = 1500
dpPrint = []
dqPrint = []
mod_dp_diff = []
mod_dq_diff = []
dpPrintDiff = []
dqPrintDiff = []

for k in range(d.nb):
    DPant.append(FuzzyMath([0,0,0]))
    DQant.append(FuzzyMath([0,0,0]))
    diffP.append(FuzzyMath([0,0,0]))

while (iter < maxiter):
    
    # -------------- CALCULO DO SUBPROBLEMA ATIVO OU P-TETA ---------------
    # Cálculo do vetor de equações básicas (resíduos) DELTA_P
    # OBS: lembrar que injeção de potência ativa especificada = Pliq
    
    #Marcação de início de iteração no log
    log.write_log_iter_mark(iter) 
    log.write_log_space()
    log.write_log(">>> P liq")
    log.write_log_list_fuzzy(d.pliq)
    log.write_log_space()
    log.write_log(">>> Q liq")
    log.write_log_list_fuzzy(d.qliq)
    log.write_log_space()
    
    DP = []
    DP2 = []
    PCALC = []
    
    for k in range(0,d.nb):
        DP.append(d.pliq[k] - ((vb[k]**2) * float(G[k,k])))
        PCALC.append((vb[k]**2) * (-1) * float(G[k,k]))

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python       
        dt = ab[k] - ab[m]
        dtm = ab[m] - ab[k]
        DP[k] = DP[k] - (dt.cos()*float(G[k,m]) + dt.sen()*float(B[k,m])) * (vb[k]) * (vb[m])
        DP[m] = DP[m] - (dtm.cos()*float(G[m,k])+ dtm.sen()*float(B[m,k])) * (vb[m]) * (vb[k])

        PCALC[k] = PCALC[k] - (dt.cos()*float(G[k,m]) + dt.sen()*float(B[k,m])) * (vb[k]) * (vb[m])
        PCALC[m] = PCALC[m] - (dtm.cos()*float(G[m,k])+ dtm.sen()*float(B[m,k])) * (vb[m]) * (vb[k])

    for k in range(d.nb):
        DP2.append(d.pliq[k]-PCALC[k])

    #for i in range(d.nb):
    #    DP.append(d.pliq[i]-PCALC[i])

    '''
    for i in range(d.nb):
        for r in range(d.nr):
            k = d.bini[r] 
            m = d.bfim[r] 
            k-=1
            m-=1
            dt = ab[k] - ab[m]
            if i == k:
    '''
    
    log.write_log(">>> DP ")
    for i in range(d.nb):
        log.write_log(str(DP[i].f))

    #calcula DP
    #for k in range(0,d.nb):
    #    DP.append(d.pliq[k] - PCALC[k])

    # Artifício para a barra V-Teta não interfir no teste de convergência
    for k in range(d.nb):
        if d.tipo_barras[k] == 2: 
            DP[k] = FuzzyMath(np.array([0,0,0]))

    ###############################################
    for k in range(d.nb):
        #ant = DPant[k].f[2]- DPant[k].f[0]
        #dep = DP[k].f[2] - DP[k].f[0]
        #diffP[k] = dep - ant
        diffP[k] = DPant[k] - DP[k]

    log.write_log(">>> diff DP ")
    for i in range(d.nb):
        log.write_log(str(diffP[i].f))
    log.write_log(">>> teste conv DP ")
    for i in range(d.nb):
        log.write_log(str(((DP[i].f[2]-DP[i].f[0])-(DPant[i].f[2]-DPant[i].f[0]))/2))
    ################################################

    # Teste de convergência e obtenção de nova estimativa para os ângulos
    # de tensões de barra
    mod_dp=[]
    for i in range(d.nb):
        #if d.tipo_barras[k]!= 2:
        mod_dp_diff.append(((DP[i].f[2]-DP[i].f[0])-(DPant[i].f[2]-DPant[i].f[0]))/2)
        #mod_dp.append(((DP[i].f[2]-DP[i].f[0])-(DPant[i].f[2]-DPant[i].f[0]))/2)
        mod_dp.append((DP[i].f[1]-DPant[i].f[1]))
        #mod_dp.append(DP[i].f[1])
    dpPrint.append(np.amax(np.absolute(mod_dp))) 
    dpPrintDiff.append(np.amax(np.absolute(mod_dp_diff)))   
    print('dp max = {0}'.format(np.amax(np.absolute(mod_dp))))
    if np.amax(np.absolute(mod_dp)) <= errP: # Teste de convergência subproblema P-Teta 
        Kp = 0
        if Kq == 0:
            break  # Sai do processo iterativo
    else:  #Correção dos ângulos de tensões de barra
        #d.ab = d.ab + np.matmul(B1L,DP)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        #--

    #    log.write_log(">>> DP")
    #    for i in range(0,6):
    #        log.write_log(str(DP[i].f))

        
        dx = []
        som, mom, lom = np.zeros((d.nb,1)), np.zeros((d.nb,1)), np.zeros((d.nb,1))
        aux = np.zeros((d.nb,1))
        #vbA = FuzzyMath.convToArray(vb)
        #abA = FuzzyMath.convToArray(ab)

        ud = UniversoDiscurso(DP, G, B, vb, ab, d.nb, d.nr, d.bini, d.bfim, d.tipo_barras,'ativo', iter)
        dfmax, dxmax = ud.calc_dfmax_dxmax()       
        print(">>  dfmax={0} >>> dxmax={1}".format(dfmax,dxmax)) 
        f = FuzzyInfSystem(dfmax, dxmax, presF, presX)
        x = f.pert_funcs_df()
        y = f.pert_funcs_dx()
        for k in range(0,d.nb):
            if np.sum(DP[k].f) != 0:
                dp = fuzz.trimf(f.uni_dis_F, np.sort(np.round(np.array(DP[k].f),1)))
                act_mfs = f.activate_mfs(dp, x)
                mfs_saida = f.calc_mfs_saida(act_mfs, y)
                mfs_agregadas = f.agregar_mfs_saida(mfs_saida)
                #dx.append(f.calc_centroide(mfs_agregadas))
                som[k]=fuzz.defuzzify.defuzz(f.uni_dis_X,mfs_agregadas,'som')
                mom[k]=fuzz.defuzzify.defuzz(f.uni_dis_X,mfs_agregadas,'mom')
                lom[k]=fuzz.defuzzify.defuzz(f.uni_dis_X,mfs_agregadas,'lom')
            else:
                dx.append(0)
                som[k]=0
                mom[k]=0
                lom[k]=0
            aux[k] = DP[k].f[1]
        DPAcum.append(np.amax(np.absolute(aux)))
        DXPAcum.append(np.amax(np.absolute(mom)))
        log.write_log_space()                
        log.write_log(">>> DX  >>> P >>> delta theta")
        log.write_log(str(dx))
        log.write_log_space() 
        if iter == 15:
            print('iteração {0} '.format(iter))

        #atualizar os angulos
        for k in range(0,d.nb):
            ab[k].f[0] = ab[k].f[0] + som[k] #dx[k]
            ab[k].f[1] = ab[k].f[1] + mom[k]  #mom[k]
            ab[k].f[2] = ab[k].f[2] + lom[k] #dx[k]
            #d.e[k] =d.e[k]+ dx[k]
        #--

        log.write_log_space()
        log.write_log(">>> Ângulos")
        log.write_log(str(FuzzyMath.convToArray(ab)*180/np.pi))
        log.write_log_space()

        p = p + 0.5
    

    # -------------- CALCULO DO SUBPROBLEMA REATIVO OU Q-V ----------------
    # Cálculo do vetor de equaçõees básicas (resíduos) DELTA_Q
    # OBS: lembrar que injeição de potência reativa especificada = Qliq

    #DQ = np.zeros((d.nb,1))
    DQ = []
    #QCALC = []

    for k in range(0,d.nb):
            DQ.append(d.qliq[k] + (vb[k]**2)*float(B[k,k]))

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python
        dt = ab[k] - ab[m]
        dtm = ab[m] - ab[k]
        DQ[k] = DQ[k] - (vb[k]*vb[m])*(dt.sen()*float(G[k,m])-dt.cos()*float(B[k,m]))
        DQ[m] = DQ[m] - (vb[m]*vb[k])*(dtm.sen()*float(G[m,k])-dtm.cos()*float(B[m,k]))


    log.write_log(">>> DQ ")
    for i in range(d.nb):
        log.write_log(str(DQ[i].f))

    #calcula DQ
    #for k in range(0,d.nb):
    #    DQ.append(d.qliq[k] - QCALC[k])

    # Artifício para as barras V-Teta e P-V não interfir no teste de convergência
    for k in range(0,d.nb):
        if d.tipo_barras[k] != 0:
            DQ[k] = FuzzyMath(np.array([0,0,0]))
    
    # Teste de convergÊncia e obtenção de nova estimativa para os módulos de tensões de barra
    mod_dq=[]
    for i in range(0, d.nb):
        #if d.tipo_barras[k]==0:
        mod_dq_diff.append(((DQ[i].f[2]-DQ[i].f[0])-(DQant[i].f[2]-DQant[i].f[0]))/2)
        #mod_dq.append(((DQ[i].f[2]-DQ[i].f[0])-(DQant[i].f[2]-DQant[i].f[0]))/2)
        mod_dq.append((DQ[i].f[1]-DQant[i].f[1]))
        #mod_dq.append(DQ[i].f[1])
    dqPrint.append(np.amax(np.absolute(mod_dq)))
    dqPrintDiff.append(np.amax(np.absolute(mod_dq_diff)))
    print('dq max = {0}'.format(np.amax(np.absolute(mod_dq))))
    if np.amax(np.absolute(mod_dq)) <= errQ: # Teste de convergência subproblema Q-V 
        Kq = 0
        if Kp == 0:
            break #Sai do processo iterativo
    else: #Correção dos módulos de tensões de barra
        #d.vb = d.vb + np.matmul(B2L,DQ) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
        #--
               
        log.write_log(">>> DX  >>> Q >>> delta v ")
        log.write_log(str(dx))


        dx = []
        som, mom, lom = np.zeros((d.nb,1)), np.zeros((d.nb,1)), np.zeros((d.nb,1))
        aux = np.zeros((d.nb))
        #vbA = FuzzyMath.convToArray(vb)
        #abA = FuzzyMath.convToArray(ab)

        ud = UniversoDiscurso(DQ, G, B, vb, ab, d.nb, d.nr, d.bini, d.bfim, d.tipo_barras,'reativo',iter)
        dfmax, dxmax = ud.calc_dfmax_dxmax()        
        print(">>  dfmax={0} >>> dxmax={1}".format(dfmax,dxmax)) 
        f = FuzzyInfSystem(dfmax,dxmax,presF,presX)
        x = f.pert_funcs_df()
        y = f.pert_funcs_dx()
        for k in range(0,d.nb):
            if np.sum(DQ[k].f) != 0:
                dq = fuzz.trimf(f.uni_dis_F, np.sort(np.round(np.array(DQ[k].f),1)))
                act_mfs = f.activate_mfs(dq, x)
                mfs_saida = f.calc_mfs_saida(act_mfs, y)
                mfs_agregadas = f.agregar_mfs_saida(mfs_saida)
                #dx.append(f.calc_centroide(mfs_agregadas))
                som[k]=fuzz.defuzzify.defuzz(f.uni_dis_X,mfs_agregadas,'som')
                mom[k]=fuzz.defuzzify.defuzz(f.uni_dis_X,mfs_agregadas,'mom')
                lom[k]=fuzz.defuzzify.defuzz(f.uni_dis_X,mfs_agregadas,'lom')
            else:
                dx.append(0)
                som[k]=0
                mom[k]=0
                lom[k]=0
            aux[k] = DQ[k].f[1]
        DQAcum.append(np.amax(np.absolute(aux)))
        DXQAcum.append(np.amax(np.absolute(mom)))

        #atualizar os módulos de tensão
        for k in range(0,d.nb):
            vb[k].f[0] = vb[k].f[0] + som[k] #dx[k]
            vb[k].f[1] = vb[k].f[1] + mom[k] #dx[k]
            vb[k].f[2] = vb[k].f[2] + lom[k] #lom[k] #dx[k]

            #d.f[k] =d.f[k]+ dx[k]
        #--             

        log.write_log(">>> Módulo de Tensão")
        log.write_log(str(FuzzyMath.convToArray(vb)))


        q = q + 0.5
    
    iter = p + q
    print('iter={}'.format(iter))
    DPant = DP
    DQant = DQ

log.write_log_space()
log.write_log(">>> Ângulos em graus")
#log.write_log(str(FuzzyMath.convToArray(ab)*180/np.pi))
log.write_log(str(np.round(FuzzyMath.convToArray(ab)*180/np.pi,5)))
log.write_log_space()

log.write_log(">>> Módulo de Tensão")
log.write_log(str(np.round(FuzzyMath.convToArray(vb),5)))

#fechando arquivo de log
log.close_file()

x = np.arange(len(dpPrint))
plt.plot(x, dpPrint, 'b', linewidth=1.5)
plt.title('DP - DPAnt')
plt.show()

y = np.arange(len(dqPrint))
plt.plot(y, dqPrint, 'b', linewidth=1.5)
plt.title('DQ - DQAnt')
plt.show()

x = np.arange(len(dpPrintDiff))
plt.plot(x, dpPrintDiff, 'b', linewidth=1.5)
plt.title('base - base anterior /2 de DP')
plt.show()

y = np.arange(len(dqPrintDiff))
plt.plot(y, dqPrintDiff, 'b', linewidth=1.5)
plt.title('base - base anterior /2 de DQ')
plt.show()

x = np.arange(len(DPAcum))
plt.plot(x, DPAcum, 'b', linewidth=1.5)
plt.title('DP')
plt.show()

y = np.arange(len(DXPAcum))
plt.plot(y, DXPAcum, 'b', linewidth=1.5)
plt.title('mom para DP')
plt.show()

x = np.arange(len(DQAcum))
plt.plot(x, DQAcum, 'b', linewidth=1.5)
plt.title('DQ')
plt.show()

y = np.arange(len(DXQAcum))
plt.plot(y, DXQAcum, 'b', linewidth=1.5)
plt.title('mom para DQ')
plt.show()

# ===================== CÁLCULO DO SUBSISTEMA 2 ===========================
# OBS: Pliq de barras PV e PQ, e Qliq de barras PQ, são dados do problema
# mas sãoo aqui recalculados para se certificar que os resultados obtidos na
# solução do subsistema 1 estáo corretos (se igual aos dados -> ok)
'''
Qsh_k = np.zeros((d.nb,1))
Pkm = np.zeros((d.nr,1))
Qkm = np.zeros((d.nr,1))
Pmk = np.zeros((d.nr,1))
Qmk = np.zeros((d.nr,1))
Pperdas = np.zeros((d.nr,1))
Qperdas = np.zeros((d.nr,1))
for k in range(0,d.nb):
    Pliq[k] = G[k,k]*np.square(d.vb[k])*d.sbase
    Qliq[k] = -B[k,k]*np.square(d.vb[k])*d.sbase
    Qsh_k[k] = np.square(d.vb[k])*bsh_k[k]*d.sbase
for r in range(0,d.nr):
    k = d.bini[r] 
    m = d.bfim[r] 
    k = k-1 #correção da dimensao para python
    m = m-1 #correção da dimensao para python
    tp = akm[r]
    dt = d.ab[k] - d.ab[m]
    #Calculo das injeções líquidas de potência ativa e reativa de barra
    Pliq[k] = Pliq[k] + (d.vb[k]*d.vb[m]*(G[k,m]*np.cos(dt)+B[k,m]*np.sin(dt)))*d.sbase
    Qliq[k] = Qliq[k] + (d.vb[k]*d.vb[m]*(G[k,m]*np.sin(dt)-B[k,m]*np.cos(dt)))*d.sbase
    Pliq[m] = Pliq[m] + (d.vb[m]*d.vb[k]*(G[m,k]*np.cos(-dt)+B[m,k]*np.sin(-dt)))*d.sbase
    Qliq[m] = Qliq[m] + (d.vb[m]*d.vb[k]*(G[m,k]*np.sin(-dt)-B[m,k]*np.cos(-dt)))*d.sbase
    # Calculo dos fluxos de potência ativa e reativa de ramos
    Pkm[r] = ((tp*np.square(d.vb[k]))*gkm[r] - tp*d.vb[k]*d.vb[m]*gkm[r]*np.cos(dt)- tp*d.vb[k]*d.vb[m]*bkm[r]*np.sin(dt))*d.sbase
    Qkm[r] = (-(tp*np.square(d.vb[k]))*(bsh_km[r]+bkm[r]) + tp*d.vb[k]*d.vb[m]*bkm[r]*np.cos(dt) - tp*d.vb[k]*d.vb[m]*gkm[r]*np.sin(dt))*d.sbase
    Pmk[r] = (np.square(d.vb[m])*gkm[r] - tp*d.vb[m]*d.vb[k]*gkm[r]*np.cos(dt) + tp*d.vb[m]*d.vb[k]*bkm[r]*np.sin(dt))*d.sbase
    Qmk[r] = (-np.square(d.vb[m])*(bsh_km[r]+bkm[r]) + tp*d.vb[m]*d.vb[k]*bkm[r]*np.cos(dt)+tp*d.vb[m]*d.vb[k]*gkm[r]*np.sin(dt))*d.sbase
    #Calculo das perdas de ramos
    Pperdas[r] = Pkm[r]+Pmk[r]
    Qperdas[r] = Qkm[r]+Qmk[r]
#Imprimir vetor de tensão de barra e ângulo de tensao
print(np.around(d.vb,4))
print(np.around((180/np.pi)*d.ab,3))
'''