function [mdotf_kgs] = Fuel_Flow_Rate(H, Mach, Treq_N)


T0 = 288.15;
 

[~, T, ~, ~] = AtmosProp(H);

Relative_Temp = T/T0;

[C_tp] = get_TSFC(Mach,Relative_Temp);

mdotf_kgs = C_tp * Treq_N;


end