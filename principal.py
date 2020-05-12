##############################################################################
#  Algoritmo desacoplado rápido                                              #
#  Eng. Luiz Gonzaga Rocha Junior                                            #
#  contato@luizrocha.eng.br                                                  #
##############################################################################

from DadosEntrada import *

#Dados de entrada
filein = 'sis3.dat'

# Instancia objeto dados
d = DadosEntrada('entradas/{0}'.format(filein))

#carrega arquivo de entrada para o objeto dados
d.carregar_dados()

# normalização em pu dos dados de barra
Pliq = (d.barras[:,4]-d.barras[:,6])/d.sbase # P especificado
Qliq = (d.barras[:,5]-d.barras[:,7])/d.sbase # Q especificado
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

while (iter < 100):
    
    # -------------- CALCULO DO SUBPROBLEMA ATIVO OU P-TETA ---------------
    # Cálculo do vetor de equações básicas (resíduos) DELTA_P
    # OBS: lembrar que injeção de potência ativa especificada = Pliq

    DP = np.zeros((d.nb,1))

    for k in range(0,d.nb):
            DP[k] = Pliq[k] - G[k,k]*(np.square(d.vb[k]))

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python
        dt = d.ab[k] - d.ab[m]
        DP[k] = DP[k] - d.vb[k]*d.vb[m]*(G[k,m]*np.cos(dt)+B[k,m]*np.sin(dt))
        DP[m] = DP[m] - d.vb[m]*d.vb[k]*(G[m,k]*np.cos(-dt)+B[m,k]*np.sin(-dt))
    
    
    # Artifício para a barra V-Teta não interfir no teste de convergência
    for k in range(0,d.nb):
        if d.tipo_barras[k] == 2: 
            DP[k] = 0

    # Teste de convergência e obtenção de nova estimativa para os ângulos
    # de tensões de barra
    if np.amax(np.absolute(DP)) <= 1e-4: # Teste de convergência subproblema P-Teta 
        Kp = 0
        if Kq == 0:
            break  # Sai do processo iterativo
    else:  #Correção dos ângulos de tensões de barra
        d.ab = d.ab + np.matmul(B1L,DP)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        p = p + 0.5
    

    # -------------- CALCULO DO SUBPROBLEMA REATIVO OU Q-V ----------------
    # Cálculo do vetor de equaçõees básicas (resíduos) DELTA_Q
    # OBS: lembrar que injeição de potência reativa especificada = Qliq

    DQ = np.zeros((d.nb,1))

    for k in range(0,d.nb):
        DQ[k] = Qliq[k] + B[k,k]*np.square(d.vb[k])

    for r in range(0,d.nr):
        k = d.bini[r] 
        m = d.bfim[r] 
        k = k-1 #correção da dimensao para python
        m = m-1 #correção da dimensao para python
        dt = d.ab[k] - d.ab[m]
        DQ[k] = DQ[k] - d.vb[k]*d.vb[m]*(G[k,m]*np.sin(dt)-B[k,m]*np.cos(dt))
        DQ[m] = DQ[m] - d.vb[m]*d.vb[k]*(G[m,k]*np.sin(-dt)-B[m,k]*np.cos(-dt))
 
    # Artifício para as barras V-Teta e P-V não interfir no teste de convergência
    for k in range(0,d.nb):
        if d.tipo_barras[k] != 0:
            DQ[k] = 0.0
    
    # Teste de convergÊncia e obtenção de nova estimativa para os módulos de tensões de barra
    if np.amax(np.absolute(DQ)) <= 1e-4: # Teste de convergência subproblema Q-V 
        Kq = 0
        if Kp == 0:
            break #Sai do processo iterativo
    else: #Correção dos módulos de tensões de barra
        d.vb = d.vb + np.matmul(B2L,DQ) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
              
        q = q + 0.5
    
    iter = p + q


# ===================== CÁLCULO DO SUBSISTEMA 2 ===========================
# OBS: Pliq de barras PV e PQ, e Qliq de barras PQ, são dados do problema
# mas sãoo aqui recalculados para se certificar que os resultados obtidos na
# solução do subsistema 1 estáo corretos (se igual aos dados -> ok)

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