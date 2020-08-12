import numpy as np

class FuzzyMath:
    f = None
    pertf = None

    '''Classe que implementa a artimética fuzzy para números fuzzy com funções de pertinência triangulares'''
    #  http://computacaointeligente.com.br/artigos/operacoes-aritmeticas-entre-numeros-fuzzy/
    
    def __init__(self, *args):
        nargs = len(args)
        if nargs == 0:
            self.f = np.array([0,0,0],dtype=float)
            self.pertf = 0
        elif nargs == 1:
            self.f = args[0]
            self.pertf = 1
        elif nargs == 2:
            self.f = args[0]
            self.pertf = args[1]
        else:
            print('Máximo número de argumentos 2')
            raise ValueError

    def __str__ (self):
        return '< %s; %s>' % (self.f, self.pertf)

    def __add__ (self, mf):
        s = FuzzyMath()
        s.f[0] = self.f[0] + mf.f[0]
        s.f[1] = self.f[1] + mf.f[1]
        s.f[2] = self.f[2] + mf.f[2]
        s.pertf = (self.pertf + mf.pertf) - (self.pertf * mf.pertf)
        #s.f = np.sort(s.f)
        return s
    
    def  __sub__ (self, mf):
        s = FuzzyMath()
        s.f[0] = self.f[0] - mf.f[2]
        s.f[1] = self.f[1] - mf.f[1]
        s.f[2] = self.f[2] - mf.f[0]
        s.pertf = (self.pertf - mf.pertf) + (self.pertf * mf.pertf)
        return s
    '''
    def __sub__ (self, mf):
        s = FuzzyMath()
        s.f[0] = np.min([self.f[0] - mf.f[0], self.f[2] - mf.f[2]])
        s.f[2] = np.max([self.f[0] - mf.f[0], self.f[2] - mf.f[2]])
        s.f[1] = (s.f[0] + s.f[2])/2
        s.pertf = (self.pertf - mf.pertf) + (self.pertf * mf.pertf)
        return s
    '''

    def __mul__ (self, p):
        s = FuzzyMath()
        if (type(p) == float) or (type(p) == int):
            #if p < 0:
            #    print("Erro : Valor necessita ser maior que zero")
            #    raise ValueError
            #else:
            s.f[0] = self.f[0] * p
            s.f[1] = self.f[1] * p
            s.f[2] = self.f[2] * p
            #s.pertf = 1 - (1-self.pertf)**p
        else:
            aX = np.absolute((self.f[2]-self.f[0])/2)
            aY = np.absolute((p.f[2]-p.f[0])/2)
            mX = self.f[1]
            mY = p.f[1]
            mA = mX*mY
            aOrd =np.sort([np.absolute(mA - ((mX - aX )*(mY - aY))),np.absolute(mA - ((mX - aX )*(mY + aY ))),np.absolute(mA - ((mX + aX )*(mY - aY ))),np.absolute(mA - ((mX + aX )*(mY + aY )))])
            s.f[1] = mA
            s.f[0] = mA - aOrd[1]
            s.f[2] = mA + aOrd[1]
            s.pertf = self.pertf * p.pertf
        return s

    def __truediv__ (self, mf):
        s = FuzzyMath()
        s.f[0] = self.f[0] / mf.f[0]
        s.f[1] = self.f[1] / mf.f[1]
        s.f[2] = self.f[2] / mf.f[2]
        s.pertf = self.pertf / mf.pertf
        return s

    def __pow__ (self, lamb):
        s = FuzzyMath()
        s.f = self.f ** lamb
        s.pertf = self.pertf ** lamb
        return s
    
    def cos(self):
        s = FuzzyMath()
        m = self.f[1]
        a = np.absolute((self.f[2]-self.f[0])/2)
        mcos   =  (np.cos(m - a) + np.cos(m + a)) / 2
        acos   = (np.absolute((np.absolute(mcos)-np.absolute(np.cos(m - a)))) + np.absolute((np.absolute(mcos) + np.absolute(np.cos(m + a)))))/2
        s.f[1] = mcos
        s.f[0] = mcos - acos
        s.f[2] = mcos + acos
        #s.pertf = self.pertf * p.pertf
        return s

    def sen(self):
        s = FuzzyMath()
        m = self.f[1]
        a = np.absolute((self.f[2]-self.f[0])/2)
        msin   =  (np.sin(m - a) + np.sin(m + a)) / 2
        asin   = (np.absolute((np.absolute(msin)-np.absolute(np.sin(m - a)))) + np.absolute((np.absolute(msin) + np.absolute(np.sin(m + a)))))/2
        s.f[1] = msin
        s.f[0] = msin - asin
        s.f[2] = msin + asin
        #s.pertf = self.pertf * p.pertf
        return s
    
    def convToArray(self):
        s = np.zeros((len(self),3))
        for k in range(len(s)):
            s[k,0] = self[k].f[0]
            s[k,1] = self[k].f[1]
            s[k,2] = self[k].f[2]  
        return s