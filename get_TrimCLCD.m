function [CL_trim,CD_trim] = get_TrimCLCD(FlightPathAngle_Deg,W_N,Mach,H)

S = 338.9;
C_D0 = 0.027;
K = 0.05;
gamma = 1.4;

 %retrieve atmos properties
[p,~,~,~] = AtmosProp(H);

%determine dynamic pressure q
q = 0.5*gamma*p*Mach.^2;

%determine lift required
Lreq = W_N*cos(FlightPathAngle_Deg./(180/pi));

%determine lift coefficient required
CL_trim = Lreq/(q*S);

%determine frag coefficient required
CD_trim = C_D0 + K*CL_trim.^2;



end