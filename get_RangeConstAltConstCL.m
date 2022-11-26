function [Range2] = get_RangeConstAltConstCL(BFR)

H = 360*30.48;
Cd0 = 0.0221;
K = 0.0259;
T0 = 288.15;
g0 = 9.8065;
Wmtow_N = 217000*g0;
Wi = Wmtow_N;
S = 363.1;
gamma = 1.4;
R = 287.05287; 
% Mach = 0.85; %----- needs to be calculated for ctp
We = Wi - (BFR*Wi);

[~, Ts, rhos, ~] = AtmosProp(H);

CL = sqrt(Cd0/(3*K));

CD = Cd0 + K*CL^2;

q = Wi/(S*CL);

Mach = sqrt((2*q)/(rhos*gamma*R*Ts));

[C_tp] = get_TSFC(Mach, Ts/T0);

Range2 = (2/(g0*C_tp))*(((2*Wi)/(rhos*S*CL))^0.5)*(CL/CD)*(1-(Wi/We)^-0.5);







end