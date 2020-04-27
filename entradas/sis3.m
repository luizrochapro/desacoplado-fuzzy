Sbase = 100. ; % em MVA

% bar  Tip  V ang  Pg    Qg    Pc    Qc    Qsh
% Pg,Pc em MW   ;   Qg,Qc,Qshunt em Mvar, a,b,c,d,fuzzy
barras = [
 1  2  1.000  0.0    0.0  0.0    0.0     0.0    0.0  0 
 2  1  1.000  0.0   50.0 10.0    0.0     0.0    0.0  0  
 3  0  1.000  0.0    0.0  0.0  200.0     0.0    0.0  1 
];

% b_inicial  b_final  R  X  Q_sh  tap
% R,X em %  ;   b_sh em Mvar
ramos = [
 1   2   0.2   2.0    40.0  1.00
 2   3   0.2   2.0    40.0  1.00
 1   3   0.2   2.0    40.0  1.00
];

% Geração e Carga fuzzy - trapezoidal
% PG QG PL QL
fuzzy = [
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
3 0 0 0 0 0 0 0 0 190.0 195.0 205.0 210.0 0 0 0 0 
];