##############################################################################
#  Algoritmo Newton Raphson                                                  #
#  Eng. Luiz Gonzaga Rocha Junior                                            #
#  contato@luizrocha.eng.br                                                  #
##############################################################################
#coding: utf-8

from DadosEntrada import *
#from Log import *

def newton(filein):
    #Arquivo de Log
    #log = Log('saida.log')
    #log.open_file()

    #Dados de entrada
    #filein = 'sis14.dat'
    #filein = 'sis6.dat'

    # Instancia objeto dados
    d = DadosEntrada('entradas/{0}'.format(filein))

    #carrega arquivo de entrada para o objeto dados
    d.carregar_dados()

    #calcula as potências liquidas das barras em pu
    p_esp = d.calc_pliq_newton()/d.sbase
    q_esp = d.calc_qliq_newton()/d.sbase

    # normalização em pu dos dados de shunt
    #bsh_k = d.barras[:,8]/d.sbase

    #tipos de barras
    npq = len(d.tipo_barras[d.tipo_barras==0])
    npvpq = len(d.tipo_barras[d.tipo_barras==0]) + len(d.tipo_barras[d.tipo_barras==1])


    # Formaçaao da matriz Z
    z = np.zeros((d.nr),dtype=complex)
    b = np.zeros((d.nr),dtype=complex)
    y = np.zeros((d.nr),dtype=complex)
    Y = np.zeros((d.nb,d.nb),dtype=complex)
    G = np.zeros((d.nb,d.nb))  
    B = np.zeros((d.nb,d.nb)) 
    a = d.ramos[:,5]

    for k  in range(0,d.nr):
        z[k] = complex(d.ramos[k][2],d.ramos[k][3])
        b[k] = complex(d.ramos[k][4])

    y = 1./z  # inverso de cada elemento

    # Elementos fora da diagonal da matriz Y
    for k in range(0, d.nr):
        Y[d.bini[k]-1,d.bfim[k]-1] = Y[d.bini[k]-1,d.bfim[k]-1] - y[k]/a[k]
        Y[d.bfim[k]-1,d.bini[k]-1] = Y[d.bini[k]-1,d.bfim[k]-1]

    # formação dos elementos da diagonal da matriz Y 
    for m in range(0,d.nb):
        for n in range(0,d.nr):
            if d.bini[n]-1 == m:
                Y[m,m] = Y[m,m] + y[n]/(np.square(a[n])) + b[n]
            elif d.bfim[n]-1 == m:
                Y[m,m] = Y[m,m] + y[n] + b[n]

    #matrizes G e B
    G = np.real(Y)  
    B = np.imag(Y)


    # ------- PROCESSO ITERATIVO DO MÉTODO NEWTON RAPSHON --------------------
    ## inicializações
    iter = 0
    tol = 1 # tolerância
    conv = False

    #inicio do processo iterativo
    while (iter < 500 and conv == False):
        
        # Calcular P and Q
        P = np.zeros((d.nb,1))
        Q = np.zeros((d.nb,1))
        for i in range(0, d.nb):
            for k in range(0, d.nb):
                P[i] = P[i] + d.vb[i][1] * d.vb[k][1] * (G[i,k]*np.cos(d.ab[i][1]-d.ab[k][1]) + B[i,k]*np.sin(d.ab[i][1]-d.ab[k][1]))
                Q[i] = Q[i] + d.vb[i][1] * d.vb[k][1] * (G[i,k]*np.sin(d.ab[i][1]-d.ab[k][1]) - B[i,k]*np.cos(d.ab[i][1]-d.ab[k][1]))
    
        
        #checar limites de reativo
        if iter <= 7 and iter > 2:    # checar após a sétima iteração
            for n in range(0,d.nb):
                if d.tipo_barras[n] == 1:
                    QG = Q[n] + (d.ql[n,1]/d.sbase)
                    if QG < (d.qmin[n]/d.sbase):
                        d.vb[n,1] = d.vb[n,1] + 0.01
                    elif QG > (d.qmax[n]/d.sbase):
                        d.vb[n,1] = d.vb[n,1] - 0.01

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

        tol = np.max(np.absolute(mismatch))

        if tol <= 1e-4:
            conv = True
        
        iter+=1

    #convertendo tensão de polar para retangular
    vm = np.zeros((d.nb,1),dtype=complex)
    for k in range(0,d.nb):
        vm[k] = complex(d.vb[k][1]*np.cos(d.ab[k][1]),d.vb[k][1]*np.sin(d.ab[k][1]))

    # angulos em graus
    v_ang = np.zeros((d.nb,1))
    for k in range(0,d.nb):
        v_ang[k] = (180 * d.ab[k][1])/np.pi  


    Iij = np.zeros((d.nb,d.nb), dtype=complex)
    Iijm = np.zeros((d.nb,d.nb), dtype=complex)
    Iija = np.zeros((d.nb,d.nb), dtype=complex)
    Sij = np.zeros((d.nb,d.nb), dtype=complex)
    Pij = np.zeros((d.nb,d.nb), dtype=complex)
    Qij = np.zeros((d.nb,d.nb), dtype=complex)
    Si = np.zeros((d.nb,1), dtype=complex)

    # Injeções de correntes nas barras
    I = np.matmul(Y, vm)
    Im = np.abs(I)
    Ia = np.angle(I)

    # fluxo de corrente nas linhas
    for m in range(0,d.nr):
        p = d.bini[m]
        q = d.bfim[m]
        p -=1
        q -=1
        Iij[p,q] = -(vm[p] - vm[q]) * Y[p,q] # Y(m,n) = -y(m,n)..
        Iij[q,p] = -Iij[p,q]

    Iijm = np.abs(Iij)
    Iija = np.angle(Iij)

    # fluxo de potência nas Linhas
    for m in range(0, d.nb):
            for n in range(0,d.nb):
                if m != n:
                    Sij[m,n] = vm[m] * np.conj(Iij[m,n]) * d.sbase

    Pij = np.real(Sij)
    Qij = np.imag(Sij)

    # Perdas nas linhas
    Lij = np.zeros((d.nr,1), dtype = complex)
    for m in range(0, d.nr):
        p = d.bini[m]
        q = d.bfim[m]
        p -=1
        q -=1
        Lij[m] = Sij[p,q] + Sij[q,p]

    Lpij = np.real(Lij)
    Lqij = np.imag(Lij)

    # Injeção de potência nas barras
    for i in range(0,d.nb):
        for k in range(0,d.nb):
            Si[i] = Si[i] + np.conj(vm[i]) * vm[k] * Y[i,k] * d.sbase

    Pi = np.real(Si)
    Qi = -np.imag(Si)
    Pg = Pi + d.pl[:,1].reshape(-1,1)
    Qg = Qi + d.ql[:,1].reshape(-1,1)

    if conv:
        print('Sistema convergiu em {0} iterações'.format(iter))
    else:
        print('Sistema não convergente')
        print('{} iterações'.format(iter))


    return d.vb, d.ab, J, Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Y, d , npvpq, npq