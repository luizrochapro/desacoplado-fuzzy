from FuzzyInfSystem import *
import matplotlib.pyplot as plt



f = FuzzyInfSystem(10,1,0.1,0.01)

x = f.pert_funcs_df()
y = f.pert_funcs_dx()

dp = fuzz.trimf(f.uni_dis_F, [0, 5, 10])

act_mfs = f.activate_mfs(dp, x)

mfs_saida = f.calc_mfs_saida(act_mfs, y)

mfs_agregadas = f.agregar_mfs_saida(mfs_saida)

dx = f.calc_centroide(mfs_agregadas)

print(dx)