function [C_tp] = get_TSFC(Mach,Relative_Temp)

% Declare known constants
n2 = 0.508;     % TSFC exponent, -
ct2 = 0.725;  % TSFC constant, s-1
g0 = 9.8065; % Gravitational accelerational at sea-level

% Determine TSFC at this relative temp, Mach, kg/s/N
C_tp = (1/(3600*g0))*ct2*sqrt(Relative_Temp)*Mach.^n2;

end