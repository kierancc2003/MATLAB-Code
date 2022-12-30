function [VNE_ms, VME_ms, VstallTO_ms, VstallLD_ms] = get_VNEspeeds(H)

V_NE = 400;
knots_2_ms = 0.51444;
CAS_ms = V_NE * knots_2_ms;
gamma = 1.4;
g0 = 9.8065;
Wmtow_N = 283720*g0;
S = 338.9;
CL_Max_TO = 2.33;
CL_Max_LD = 2.86;
Mach_NE = 0.92;
R = 287.05287;


[p0, ~, ~, a0] = AtmosProp(0);
[ps, Ts, rhos, as] = AtmosProp(H);


qc = p0 * (( 1 + (gamma-1)/2 * (CAS_ms/a0)^2   )^(gamma/(gamma-1))-1);


VNE_ms= sqrt ( ( ( (qc/ps) +1) ^ ( (gamma-1) /gamma) -1)*(2*as^2/ (gamma-1)));

VstallTO_ms = sqrt( (2/rhos) * (Wmtow_N/S) * (1/CL_Max_TO));
VstallLD_ms= sqrt( (2/rhos) * (Wmtow_N/S) * (1/CL_Max_LD));
VME_ms = Mach_NE * sqrt(gamma*Ts*R);



end