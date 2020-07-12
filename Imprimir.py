#coding: utf-8
import numpy as np
class Imprimir:
    
    file = None
    path_file = None

    
    def __init__(self, path_file, d):
        self.path_file  = path_file
        self.d = d

    def write_results(self,Pij, Qij, Pi, Qi, Pg, Qg, Lpij, Lqij, Pik, Pki,Qik, Qki,PerdasPik):
        self.file = open(self.path_file,'w', encoding='utf-8') #abre arquivo
        self.file.write('{0:*^120s}'.format('VALORES DETERMINÍSTICOS') + '\n')
        self.file.write('\n')
        self.file.write('{0:*^120s}'.format('Dados de Barra') + '\n')                                                                                                                
        self.file.write('{0:^5}{1:^7}{2:^7}{3:^7}{4:^7}{5:^7}{6:^7}{7:^7}{8:^7}{9:^7}{10:^7}{11:^7}{12:^7}{13:^7}{14:^7}{15:^7}{16:^7}'.format('Barra','|','V','|','Ang','|','Pinj','|','Qinj','|','Pg','|','Qg','|','Pl','|','Ql') + '\n')
        for k in range(self.d.nb):
            self.file.write('{0:>5d}{1:^7}{2:>7.3f}{3:^7}{4:>7.2f}{5:^7}{6:>7.2f}{7:^7}{8:>7.2f}{9:^7}{10:>7.2f}{11:^7}{12:>7.2f}{13:^7}{14:>7.2f}{15:^7}{16:>7.2f}'.format(
                k+1,'|', float(self.d.vb[k,1]),'|', 180*float(self.d.ab[k,1])/np.pi, '|', 
                float(Pi[k]),'|', float(Qi[k]),'|', float(Pg[k]),'|', float(Qg[k]),'|',
                float(self.d.pl[k,1]),'|', float(self.d.ql[k,1])) + '\n')
        
        #self.file.write('********************{0:^60}********************'.format('Dados de Barra') + '\n')                                                                                                                    ==' + '\n')
        #self.file.write('{0:^10}{1:^10}{2:^10}{3:^10}{4:^10}{5:^10}{6:^10}{7:^10}{8:^10}'.format('Barra','|','va','vb','vc','|','AngA','AngB','AngC') + '\n')
        #for k in range(self.d.nb):
        #    self.file.write('{0:^10d}{1:^10}{2:^10.3f}{3:^10.3f}{4:^10.3f}{5:^10}{6:^10.2f}{7:^10.2f}{8:^10.2f}'.format(k+1,'|',float(self.d.vb[k,0]),float(self.d.vb[k,1]),float(self.d.vb[k,2]),'|',float(self.d.ab[k,0]),float(self.d.ab[k,1]),float(self.d.ab[k,2])) + '\n')

        self.file.write('{:-^120s}'.format('-') + '\n')
        self.file.write('\n')
        self.file.write('{0:*^120}'.format('  Dados de Ramos  ') + '\n')
        self.file.write('\n')
        self.file.write('{0:*^120s}'.format('  Potência ativa e Reativa [MW] e [Mvar]  ') + '\n')
        self.file.write('{:-^120s}'.format('-') + '\n')
        self.file.write('{0:^4}{1:^2}{2:^4}{3:>10}{4:^10}{5:^4}{6:^2}{7:^4}{8:>10}{9:^10}{10:^4}{11:^2}{12:^4}{13:>10}{14:^10}{15:^4}{16:^2}{17:^4}{18:>10.2}'.format('DE','-','PARA','P','|','PARA','-','DE','P','|','DE','-','PARA','Q','|','PARA','-','DE','Q')+'\n')
        for m in range(self.d.nr):
            i = self.d.bini[m]
            j = self.d.bfim[m]
            i -= 1
            j -= 1
            self.file.write('{0:^4d}{1:^2}{2:^4d}{3:>10.2f}{4:^10}{5:^4d}{6:^2}{7:^4d}{8:>10.2f}{9:^10}{10:^4d}{11:^2}{12:^4d}{13:>10.2f}{14:^10}{15:^4d}{16:^2}{17:^4d}{18:>10.2f}'.format(i+1,'-',j+1, Pij[i,j],'|',j+1,'-',i+1,Pij[j,i],'|',i+1,'-',j+1, Qij[i,j],'|',j+1,'-',i+1,Qij[j,i])+'\n')
        self.file.write(''+'\n')
        self.file.write('{:-^120s}'.format('-') + '\n')
        self.file.write(''+'\n')

        self.file.write('{0:*^120s}'.format('  Perdas ativa [MW]') + '\n')
        self.file.write('{0:^4}{1:^2}{2:^4}{3:>10}{4:^10}{5:>10}'.format('DE','-','PARA','PLoss','|','QLoss')+'\n')
        for m in range(self.d.nr):
            i = self.d.bini[m]
            j = self.d.bfim[m]
            self.file.write('{0:^4d}{1:^2}{2:^4d}{3:>10.4f}{4:^10}{5:>10.4f}'.format(i,'-',j, Lpij[m,0],'|',Lqij[m,0])+'\n')
        self.file.write(''+'\n')
        self.file.write('{:-^120s}'.format('-') + '\n')
        self.file.write('\n')
        self.file.write('{0:*^120s}'.format('VALORES FUZZY') + '\n')
        self.file.write('\n')
        self.file.write('{0:*^120s}'.format('Dados de Barra') + '\n')                                                                                                                
        self.file.write('{0:^5}{1:^7}{2:^7}{3:^7}{4:^7}{5:^7}{6:^7}{7:^7}{8:^7}{9:^7}{10:^7}{11:^7}{12:<7}'.format(
            'Barra','|','V1','|','V2','|','V3','|','Ang1','|','Ang2','|','Ang3') + '\n')
        for k in range(self.d.nb):
            self.file.write('{0:>5d}{1:^7}{2:>7.3f}{3:^7}{4:>7.3f}{5:^7}{6:>7.3f}{7:^7}{8:>7.3f}{9:^7}{10:>7.3f}{11:^7}{12:>7.3f}'.format(
                k+1,'|', float(self.d.vb[k,0]),'|', float(self.d.vb[k,1]), '|', float(self.d.vb[k,2]),'|',
                180*float(self.d.ab[k,0])/np.pi,'|', 180*float(self.d.ab[k,1])/np.pi,'|', 180*float(self.d.ab[k,2])/np.pi) + '\n')

        self.file.write('{:-^120s}'.format('-') + '\n')
        self.file.write('\n')

        self.file.write('{0:*^120s}'.format('  Potência ativa e Reativa [MW] e [Mvar]  ') + '\n')
        self.file.write('{:-^120s}'.format('-') + '\n')
        self.file.write('{0:^4}{1:^2}{2:^4}{3:>10}{4:>10}{5:>10}{6:^10}{7:^4}{8:^2}{9:^4}{10:>10}{11:>10}{12:>10}{13:^10}'.format('DE','-','PARA','P1','P2','P3','|','PARA','-','DE','P1','P2','P3','|')+'\n')
        for m in range(self.d.nr):
            i = self.d.bini[m]
            j = self.d.bfim[m]
            self.file.write('{0:^4d}{1:^2}{2:^4d}{3:>10.2f}{4:>10.2f}{5:>10.2f}{6:^10}{7:^4d}{8:^2}{9:^4d}{10:>10.2f}{11:>10.2f}{12:>10.2f}{13:^10}'.format(i,'-',j, Pik[m,0], Pik[m,1], Pik[m,2],'|',j,'-',i,Pki[m,0], Pki[m,1], Pki[m,2],'|') +'\n')
        self.file.write(''+'\n')
        self.file.write('{:-^120s}'.format('-') + '\n')

        self.file.write('{0:^4}{1:^2}{2:^4}{3:>10}{4:>10}{5:>10}{6:^10}{7:^4}{8:^2}{9:^4}{10:>10.2}{11:>10.2}{12:>10.2}'.format('DE','-','PARA','Q1','Q2','Q3','|','PARA','-','DE','Q1','Q2','Q3')+'\n')
        for m in range(self.d.nr):
            i = self.d.bini[m]
            j = self.d.bfim[m]
            self.file.write('{0:^4d}{1:^2}{2:^4d}{3:>10.2f}{4:>10.2f}{5:>10.2f}{6:^10}{7:^4d}{8:^2}{9:^4d}{10:>10.2f}{11:>10.2f}{12:>10.2f}'.format(i,'-',j,Qik[m,0],Qik[m,1],Qik[m,2],'|', j,'-', i, Qki[m,0], Qki[m,1], Qki[m,2]) +'\n')
        self.file.write(''+'\n')
        self.file.write('{:-^120s}'.format('-') + '\n')

        self.file.write('{0:*^120s}'.format('  Perdas Ativas  ') + '\n')
        self.file.write('{:-^120s}'.format('-') + '\n')
        self.file.write('{0:^4}{1:^2}{2:^4}{3:>10}{4:>10}{5:>10}{6:^10}{7:^4}{8:^2}{9:^4}{10:>10}{11:>10}{12:>10}{13:^10}'.format('DE','-','PARA','P1','P2','P3','|','PARA','-','DE','P1','P2','P3','|')+'\n')
        for m in range(self.d.nr):
            i = self.d.bini[m]
            j = self.d.bfim[m]
            self.file.write('{0:^4d}{1:^2}{2:^4d}{3:>10.5f}{4:>10.5f}{5:>10.5f}{6:^10}{7:^4d}{8:^2}{9:^4d}{10:>10.5f}{11:>10.5f}{12:>10.5f}{13:^10}'.format(i,'-',j, PerdasPik[m,0], PerdasPik[m,1], PerdasPik[m,2],'|',j,'-',i,PerdasPik[m,0], PerdasPik[m,1], PerdasPik[m,2],'|') +'\n')
        self.file.write(''+'\n')
        self.file.write('{:-^120s}'.format('-') + '\n')


        self.file.close() # fecha arquivo