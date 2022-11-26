function [SAR, SE] = getSAR_SE(W_N, H, Mach, FlightPathAngle_Deg)

T0 = 288.15; %given sea level constant 

[~,Ts,~,a] = AtmosProp(H); %use atmos props to get temperature and speed of sound ot later convert to velocity

Relative_Temp = Ts/T0; %used for get tsfc function

[C_tp] = get_TSFC(Mach,Relative_Temp); %calculating tsfc 

[CL_trim,CD_trim] = get_TrimCLCD(FlightPathAngle_Deg,W_N,Mach,H); %calculating Cl and CD trim

SE = (1/C_tp)*(1/W_N)*(CL_trim/CD_trim); %calculating SE 

SAR = ((a*Mach)/C_tp)*(CL_trim/CD_trim)*(1/W_N); %equation in notes gives L/D in the equation but is equal to Cl/CD from function


end