function [Range1] = get_RangeConstMachConstCL(BFR)

Cd0 = 0.0221;
K = 0.0259;
Mach = 0.85;
H = 360*30.48;
T0 = 288.15;
g0 = 9.8065;


[~, Ts, ~, as] = AtmosProp(H);

[C_tp] = get_TSFC(Mach, Ts/T0);

V = Mach*as;

CL = sqrt(Cd0/(3*K));

CD = Cd0 + K*CL^2;


Range1 = (V/(C_tp*g0))*(CL/CD)*log(1/(1-BFR));

end 