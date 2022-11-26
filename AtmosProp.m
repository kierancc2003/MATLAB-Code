function [ps, Ts, rhos, as] = AtmosProp(H)



g = 9.8065;
R = 287.05287;
gamma = 1.4;
L0 = -0.0065;
L11 = 0;
T0 = 288.15;
T11 = 216.65;
P0 = 101325;
P11 = 22632.559;
L20 = 0.001;
H11 = 11000;
H20 = 20000;



if H <= H11
        Ts = T0 + L0*(H);
        as = sqrt(gamma*R*Ts);
        ps = P0*((1+(H*(L0/T0)))^-(g/(R*L0)));
        rhos = ps/(R*Ts);
    elseif H > H11 && H <= H20 
        Ts = T11 + L11*(H-H11);
        as = sqrt(gamma*R*Ts);
        ps = P11*exp(-(g/(R*T11))*(H-H11));
        rhos = ps/(R*Ts);
    else
        Ts = T11 + L20*(H-H20);
        as = sqrt(gamma*R*Ts);
        P20 = P11*exp(-(g/(R*T11))*(H20-H11));
        ps = P20*exp(-(g/(R*T11))*(H-H20));
        rhos = ps/(R*Ts);
end


end
