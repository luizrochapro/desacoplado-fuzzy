##############################################################################
#  Algoritmo Desacoplado                                                     #
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
d.calc_pliq()
d.calc_qliq()

#log.write_log_div()
#log.write_log(d.pliq)
#log.write_log(d.qliq)

# normalização em pu dos dados de shunt
bsh_k = d.barras[:,8]/d.sbase

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

# ------- PROCESSO ITERATIVO DO MÉTODO DESACOPLADO RÁPIDO VERSÃO XB -------
## inicializações
iter = 0
p = 0
q = 0
Kp = 1
Kq = 1

while (iter < 30):
    
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
    DPCALC = []
    '''
    for k in range(0,d.nb):
            e = FuzzyMath(d.e[k]) # criar número fuzzy a partir dos pontos do triangulo
            f = FuzzyMath(d.f[k]) # criar número fuzzy a partir dos pontos do triangulo
            DP.append(d.pliq[k] - ((e * e) + (f * f)) * float(G[k,k]))

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python       
        ek = FuzzyMath(d.e[k])
        em = FuzzyMath(d.e[m])
        fk = FuzzyMath(d.f[k])
        fm = FuzzyMath(d.f[m])
        DP[k] = DP[k] - ek * (em * float(G[k,m]) - fm * float(B[k,m])) + fk * (fm * float(G[k,m]) + em * float(B[k,m]))
    '''
    for k in range(0,d.nb):
            e = FuzzyMath(d.e[k]) # criar número fuzzy a partir dos pontos do triangulo
            f = FuzzyMath(d.f[k]) # criar número fuzzy a partir dos pontos do triangulo
            DPCALC.append(((e * e) + (f * f)) * float(G[k,k]))

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python       
        ek = FuzzyMath(d.e[k])
        em = FuzzyMath(d.e[m])
        fk = FuzzyMath(d.f[k])
        fm = FuzzyMath(d.f[m])
        DPCALC[k] = DPCALC[k] + (ek * ((em * float(G[k,m])) - (fm * float(B[k,m]))) + fk * ((fm * float(G[k,m])) + (em * float(B[k,m]))))
        DPCALC[m] = DPCALC[m] + (em * ((ek * float(G[m,k])) - (fk * float(B[m,k]))) + fm * ((fk * float(G[m,k])) + (ek * float(B[m,k]))))

    log.write_log(">>> DP CALC")
    for i in range(0,6):
        log.write_log(str(DPCALC[i].f))

    #calcula DP
    for k in range(0,d.nb):
        DP.append(d.pliq[k] - DPCALC[k])

    # Artifício para a barra V-Teta não interfir no teste de convergência
    for k in range(0,d.nb):
        if d.tipo_barras[k] == 2: 
            DP[k] = FuzzyMath(np.array([0,0,0]))

    # Teste de convergência e obtenção de nova estimativa para os ângulos
    # de tensões de barra
    mod_dp=[]
    for k in range(0, d.nb):
        mod_dp.append(np.absolute(DP[k].f[1]))
        
    if np.amax(mod_dp) <= 1e-4: # Teste de convergência subproblema P-Teta 
        Kp = 0
        if Kq == 0:
            break  # Sai do processo iterativo
    else:  #Correção dos ângulos de tensões de barra
        #d.ab = d.ab + np.matmul(B1L,DP)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        #--

        log.write_log(">>> DP")
        for i in range(0,6):
            log.write_log(str(DP[i].f))

        
        dx = []
        ud = UniversoDiscurso(DP, G, B, d.e, d.f, d.nb, d.nr, d.bini, d.bfim, d.tipo_barras)
        dfmax,dxmax = ud.calc_dfmax_dxmax()        
        f = FuzzyInfSystem(dfmax,dxmax,0.1,0.01)
        x = f.pert_funcs_df()
        y = f.pert_funcs_dx()
        for k in range(0,d.nb):
            if np.sum(DP[k].f) != 0:
                dp = fuzz.trimf(f.uni_dis_F, np.sort(np.round(np.array(DP[k].f),1)))
                act_mfs = f.activate_mfs(dp, x)
                mfs_saida = f.calc_mfs_saida(act_mfs, y)
                mfs_agregadas = f.agregar_mfs_saida(mfs_saida)
                dx.append(f.calc_centroide(mfs_agregadas))
            else:
                dx.append(0)
        log.write_log_space()                
        log.write_log(">>> DX  >>> P >>> delta theta")
        log.write_log(str(dx))
        log.write_log_space() 
        
        #atualizar os angulos
        for k in range(0,d.nb):
            #d.ab[k] = d.ab[k] + dx[k]
            d.e[k] =d.e[k]+ dx[k]
        #--

        log.write_log_space()
        log.write_log(">>> Ângulos")
        log.write_log(str(d.ab))
        log.write_log_space()

        p = p + 0.5
    

    # -------------- CALCULO DO SUBPROBLEMA REATIVO OU Q-V ----------------
    # Cálculo do vetor de equaçõees básicas (resíduos) DELTA_Q
    # OBS: lembrar que injeição de potência reativa especificada = Qliq

    #DQ = np.zeros((d.nb,1))
    DQ = []
    DQCALC = []
    '''
    for k in range(0,d.nb):
            e = FuzzyMath(d.e[k]) # criar número fuzzy a partir dos pontos do triangulo
            f = FuzzyMath(d.f[k]) # criar número fuzzy a partir dos pontos do triangulo
            DQ.append(d.qliq[k] + ((e * e) + (f * f)) * float(B[k,k]))

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python
        ek = FuzzyMath(d.e[k])
        em = FuzzyMath(d.e[m])
        fk = FuzzyMath(d.f[k])
        fm = FuzzyMath(d.f[m])
        DQ[k] = DQ[k] - fk * (em * float(G[k,m]) - fm * float(B[k,m])) - ek * (fm * float(G[k,m]) + em * float(B[k,m]))
    '''   

    for k in range(0,d.nb):
            e = FuzzyMath(d.e[k]) # criar número fuzzy a partir dos pontos do triangulo
            f = FuzzyMath(d.f[k]) # criar número fuzzy a partir dos pontos do triangulo
            DQCALC.append(((e * e) + (f * f)) * float(-1*B[k,k]))

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python
        ek = FuzzyMath(d.e[k])
        em = FuzzyMath(d.e[m])
        fk = FuzzyMath(d.f[k])
        fm = FuzzyMath(d.f[m])
        DQCALC[k] = DQCALC[k] + fk * ((em * float(G[k,m])) - (fm * float(B[k,m]))) - ek * ((fm * float(G[k,m])) + (em * float(B[k,m])))
        DQCALC[m] = DQCALC[m] + fm * ((ek * float(G[m,k])) - (fk * float(B[m,k]))) - em * ((fk * float(G[m,k])) + (ek * float(B[m,k])))

    log.write_log(">>> DQ CALC")
    for i in range(0,6):
        log.write_log(str(DQCALC[i].f))

    #calcula DQ
    for k in range(0,d.nb):
        DQ.append(d.qliq[k] - DQCALC[k])

    # Artifício para as barras V-Teta e P-V não interfir no teste de convergência
    for k in range(0,d.nb):
        if d.tipo_barras[k] == 2:
            DQ[k] = FuzzyMath(np.array([0,0,0]))
    
    # Teste de convergÊncia e obtenção de nova estimativa para os módulos de tensões de barra
    mod_dq=[]
    for k in range(0, d.nb):
        mod_dq.append(np.absolute(DQ[k].f[1]))

    if np.amax(mod_dp) <= 1e-4: # Teste de convergência subproblema Q-V 
        Kq = 0
        if Kp == 0:
            break #Sai do processo iterativo
    else: #Correção dos módulos de tensões de barra
        #d.vb = d.vb + np.matmul(B2L,DQ) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
        #--
               
        log.write_log(">>> DX  >>> Q >>> delta v ")
        log.write_log(str(dx))

        log.write_log(">>> DQ")
        for i in range(0,6):
            log.write_log(str(DQ[i].f))

        dx = []
        ud = UniversoDiscurso(DQ, G, B, d.e, d.f, d.nb, d.nr, d.bini, d.bfim, d.tipo_barras)
        dfmax,dxmax = ud.calc_dfmax_dxmax()        
        f = FuzzyInfSystem(dfmax,dxmax,0.1,0.01)
        x = f.pert_funcs_df()
        y = f.pert_funcs_dx()
        for k in range(0,d.nb):
            if np.sum(DQ[k].f) != 0:
                dq = fuzz.trimf(f.uni_dis_F, np.sort(np.round(np.array(DQ[k].f),2)))
                act_mfs = f.activate_mfs(dq, x)
                mfs_saida = f.calc_mfs_saida(act_mfs, y)
                mfs_agregadas = f.agregar_mfs_saida(mfs_saida)
                dx.append(f.calc_centroide(mfs_agregadas))
            else:
                dx.append(0)

        #atualizar os módulos de tensão
        for k in range(0,d.nb):
            #d.vb[k] = d.vb[k] + dx[k]
            d.f[k] =d.f[k]+ dx[k]
        #--             

        log.write_log(">>> Módulo de Tensão")
        log.write_log(str(d.vb))


        q = q + 0.5
    
    iter = p + q

#fechando arquivo de log
log.close_file()
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