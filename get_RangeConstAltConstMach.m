function [Range3] = get_RangeConstAltConstMach(BFR)


Mach = 0.85;
CD0 = 0.021;
K = 0.0259;
H = 360*30.48;
g0 = 9.8065;
Wmtow_N = 217000*g0;
Wi = Wmtow_N;
We = Wi - (BFR*Wi);
S = 363.1;
T0 = 288.15;
gamma = 1.4;


[ps, Ts, ~, as] = AtmosProp(H);

[C_tp] = get_TSFC(Mach, Ts/T0);

q = 0.5*ps*gamma*Mach^2;

Range3 = (as*Mach/(g0*C_tp))*sqrt(1/(K*CD0))*(atan((Wi/(q*S))*sqrt(K/CD0))-atan(We/(q*S)*sqrt(K/CD0)));



end 


