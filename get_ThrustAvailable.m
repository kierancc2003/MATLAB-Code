function [TA0_N] = get_ThrustAvailable(H)

rho0 = 1.225;
T_A0N = 274000; %lb
% lbf2kgf = 0.454;
% T_A0_kgf = T_A0 * lbf2kgf;
% g0 = 9.8065;
% T_A0N = T_A0_kgf * g0;

Thr_Laps_Rate = 0.221;


[~,~,rho350,~] = AtmosProp(350*30.48);
[~,~,rhos,~] = AtmosProp(H);


x = rho350/rho0;

m = log(Thr_Laps_Rate)/log(x);

TA0_N = T_A0N * (rhos/rho0)^m;
end