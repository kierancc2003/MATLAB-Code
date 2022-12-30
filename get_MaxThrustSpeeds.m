function [VlowTA_ms,VhighTA_ms] = get_MaxThrustSpeeds(H)

K = 0.05;
CD0 = 0.027;
g0 = 9.8065;
Wmtow_N = 283720*g0;
S = 338.9;
Engine_No = 3;

[TA0_N] = get_ThrustAvailable(H);
[~,~,rhos,~] = AtmosProp(H);

A = 1/(rhos*CD0);
B = (Engine_No*TA0_N)/Wmtow_N;
C = Wmtow_N/S;
D = CD0 * K;



VhighTA_ms = (A * (B*C + C*sqrt(B^2 - 4*D)))^0.5;

VlowTA_ms = (A * (B*C - C*sqrt(B^2 - 4*D)))^0.5;




end