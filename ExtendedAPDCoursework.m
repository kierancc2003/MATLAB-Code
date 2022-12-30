clc;
clear;
close all;


%==================== Section of Constants ======================
 
gamma = 1.4;   % Ratio of specific heats
P0 = 101325;   % Sea level pressure
T0 = 288.15;   % Sea level temperature
rho0 = 1.225;   % Sea level air density
R = 287.05287;   % Gas constant
l1 = 2/(gamma-1); 
l2 = (gamma-1)/gamma;
C_t2 = 0.611;   % TSFC Lapse Constant                           
n2 = 0.432;   % TSFC Lapse Exponent   
g = 9.8065;   % Graviational acceleration
Wmtow_N = 217000*g;   % Aircraft max. take-off weight
Wmfuel_N = 0.358*Wmtow_N;   % Aircraft max. fuel weight
gamma = 1.4;   % Ratio of specific heats
S = 363.10;   % Wing Area, m2
rho_fuel = 840;   % Fuel Density 
PAX = 290;   % Number of Passengers
L2gallon = 4.55;   % Litres to gallon conversion rate
Km2Mile = 0.621;   % Km to Mile conversion rate
CD0 = 0.0221;
K = 0.0259;



Mach_crit = 0.85;
H_SAR = [0 100 360]'.*30.48;
W_SAR = [215 195 175]'.*1000.*g;


Norm_Op_Empty_W = 0.545;   %OEW/MTOW
Norm_Max_Struct_Pay = 0.233;   %MSP/MTOW
Norm_Max_Fuel_Cap = 0.358;   %MFC/MTOW
Fres = 0.1;   %t all times
OEW = Norm_Op_Empty_W * Wmtow_N;
MSP = Norm_Max_Struct_Pay * Wmtow_N;
MFC = Norm_Max_Fuel_Cap * Wmtow_N;
MRF = Fres * Wmtow_N;

CL_Max_TO = 2.51;
CL_Max_LD = 2.73;
V_NE = 330; %Knots
Knots_MS_Conversion = 0.51444;
V_NE_Ms = 330 * Knots_MS_Conversion;
Thr_Laps_Rate = 0.299; %@FL350 and M = 0.8
T_AO = 67500; %lb 
lb_F_Conversion = 0.454;
T_AO_Kg = T_AO * lb_F_Conversion;


Trim_Altitudes = [0:1:15000]';

%====================== Coursework 1 ========================



%==================== All Vector Creations======================= 



run DaikoScript.m

ft(:,1) = DAIKOtrajectoryv2(:, "TimeElapseds");
alt(:,1) = DAIKOtrajectoryv2(:, "Geoaltitudem");
ip(:,1) = DAIKOtrajectoryv2(:, "ImpactPressurePa");
Measure_Temp(:,1) = DAIKOtrajectoryv2(:, "MeasStaticTempK");
H = alt.Geoaltitudem(:,1);
ImpactPressure = ip.ImpactPressurePa(:,1);
FlightTime = ft.TimeElapseds(:,1);
MeasuredTemperature = Measure_Temp.MeasStaticTempK;


Ts = zeros(height(H),1); 
ps = zeros(height(H),1); 
rhos = zeros(height(H),1); 
as = zeros(height(H),1);
Mach = zeros(height(H),1);
Total_Pressure = zeros(height(H),1);
Total_Over_Static = zeros(height(H),1);
Relative_Pressure = zeros(height(H),1);
Relative_Temperature = zeros(height(H),1);
Relative_Density = zeros(height(H),1);
Relative_Measure_Temp = zeros(height(H),1);
True_AS_MeasuredT = zeros(height(H),1);
True_AS_PredictedT = zeros(height(H),1);



DeltaTime = zeros(height(H),1);
DeltaALT = zeros(height(H),1);
ABSDistance = zeros(height(H),1);



W_N = zeros(length(H(:,1)),1) + Wmtow_N; %------------------- constant for now, will derive later ---------------------


% Pre-allocate arrays. All rows represent a time step
FlightPathAngle_Deg = zeros(length(H(:,1))-1,1);    % Flight path angle between successive traj points (deg.)
mdotf_kgs = zeros(length(H(:,1))-1,1);      % Ave. fuel flow rate (kg/s)
del_Wf_kg = zeros(length(H(:,1))-1,1);      % Fuel burn per time step (kg)
Treq_N = zeros(length(H(:,1))-1,1);         % Thrust required (N)

CL_trim = zeros(height(H),1);   % Coefficent of trim lift 
CD_trim = zeros(height(H),1);   % Coefficent of trim drag 
q = zeros(height(H),1);   % Dynamic pressure
C_tp = zeros(height(H),1);   % Thrust specific fuel consumption


SAR = zeros(length(H(:,1)),1);
SE = zeros(length(H(:,1)),1);
mpg = zeros(length(H(:,1)),1);
pmpg = zeros(length(H(:,1)),1);



Mach_SAR = linspace(0.1, Mach_crit, 4537)';
SAR_3Alts = zeros(height(Mach_SAR), height(H_SAR));
SE_3Alts = zeros(height(Mach_SAR), height(H_SAR));


BFR = linspace(0,Norm_Max_Fuel_Cap,4537)';

TOW_VarBFR = zeros(height(BFR),1);
RangeVarMach_km = zeros(height(BFR),1);
RangeVarCL_km = zeros(height(BFR),1);
RangeVarH_km = zeros(height(BFR),1);

Range_VH = zeros(height(BFR),1);
Range_VM = zeros(height(BFR),1);
Range_VCL = zeros(height(BFR),1);



VNE_ms = zeros(height(Trim_Altitudes),1);
VME_ms = zeros(height(Trim_Altitudes),1);
VstallTO_ms = zeros(height(Trim_Altitudes),1);
VstallLD_ms = zeros(height(Trim_Altitudes),1);
VlowTA_ms = zeros(height(Trim_Altitudes),1);
VhighTA_ms = zeros(height(Trim_Altitudes),1);

%=============== Calculating Atmospheric Properties ==================


for i = 1:height(H)
    
    [ps(i,1 ),Ts(i,1),rhos(i,1),as(i,1)] = AtmosProp(H(i,1));
    Total_Pressure(i,:) = ImpactPressure(i,:) + ps(i,:);
    Total_Over_Static(i,:) = Total_Pressure(i,:)/ps(i,:);
    Mach(i,:) = sqrt(l1*((Total_Over_Static(i,:)^l2 )- 1));

%------------------------------------------------------------------
%Calculating the relative ratio values to be able to plot our graph
   
    Relative_Pressure(i,:) = ps(i,:)/P0;
    Relative_Temperature(i,:) = Ts(i,:)/T0;
    Relative_Density(i,:) = rhos(i,:)/rho0;
    Relative_Measure_Temp(i,:) = MeasuredTemperature(i,:)/T0;

%------------------------------------------------------------------
%Calculating the TAS from predicted and measured temperatures for given
%altitudes

    True_AS_MeasuredT(i,:) = Mach(i,:)*sqrt(gamma*R*MeasuredTemperature(i,:));
    True_AS_PredictedT(i,:) = Mach(i,:)*sqrt(gamma*R*Ts(i,:));

end


%=================== Coursework 1 Plots ======================


Time = FlightTime/60;                   %Flight Time in Minutes

figure(1);
hold on;
plot(Time,Relative_Pressure, ':', 'LineWidth',2);
plot(Time, Relative_Density, '-.', 'LineWidth',1.5);
plot(Time,Relative_Temperature, 'b', 'LineWidth',1.5);
plot(Time, Relative_Measure_Temp, 'c--', 'LineWidth',1.5);
plot(Time, Mach, 'k', 'LineWidth',1.5);
xlim([0 310]);
ylim([0 1.05]);
legend('ISA Relative Pressure \delta', 'ISA Relative Density \sigma', ...
    'ISA Relative Temperature \theta', 'Relative Measure Temperature \theta_m_e_a_s', ...
    'Mach Number {\it M}');
legend('Location', 'southoutside');
xlabel('Flight Time Elapsed {\it t_m_i_n_u_t_e_s}');
ylabel('Relative Property Value or Mach Number /\theta /\delta /\sigma /{\it M} (-)');
title('Variation of Relative Atmospheric Properties and Mach Number along Flight Path');
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'off')


figure(2);
hold on;
plot(Time,True_AS_PredictedT, 'b', 'LineWidth',1.0);
plot(Time,True_AS_MeasuredT, 'k', 'LineWidth',1.0); 
xlim([0 310]);
legend('TAS Derived From Predicted Static Temperature, {\it TAS_I_S_A}', ...
    'TAS Derived From Measured Static Temperature, {\it TAS_m_e_a_s}');
legend('Location', 'south');
title('Comparison of True Airspeed, as calculated based on ISA-Derived Temperature ', ...
    'vs. aircraft-measured static temperature, along Flight Path');
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'off')
axis(gca,'square')




%======================== Coursework 2 =====================




% NOTE: to compute things like flight path angle and fuel burn, we need to
% think about what happens BETWEEN adjacent trajectory points. Therefore,
% the new variables to consider will have 1 fewer rows than the trajectory
% point variables
%  x--------------x---------------x     i.e if this is the trajectory...
%   ^              ^               ^  
%   i             i+1             i+2
%    |___________|   |___________|
%          ^               ^
%          n              n+1

% Pre-allocate arrays. All rows represent a trajectory point


for i = 1:height(H)-1

    DeltaALT(i,:) = H(i+1,:) - H(i,:);

    DeltaTime(i,:) = FlightTime(i+1,:) - FlightTime(i,:);

    ABSDistance(i,:) = 0.5*(True_AS_PredictedT(i,:) + True_AS_PredictedT(i+1,:))*DeltaTime(i,:);

    FlightPathAngle_Deg(i,:) = (asin(DeltaALT(i,:)/ABSDistance(i,:)))*(180/pi);


end

FlightPathAngle_Deg(end+1) = FlightPathAngle_Deg(end);

%====================Calculating CL_trim and CD_trim======================

for i = 1:height(H)-1

 [CL_trim(i,:),CD_trim(i,:)] = get_TrimCLCD(FlightPathAngle_Deg(i,:),W_N(i,:),Mach(i,:),H(i,:));

end


%====================Calculating Thrust Required=========================

for i = 1:height(H)

    q(i,:) = 0.5*ps(i,:)*gamma*Mach(i,:)^2;
    Treq_N(i,:) = q(i,:)*S*CD_trim(i,:);

end

%====================Calculating Fuel Flow Rate=========================

for i = 1:height(mdotf_kgs)

  [mdotf_kgs(i,:)] = Fuel_Flow_Rate(H(i,:), Mach(i,:), Treq_N(i,:));
  del_Wf_kg(i,:) = mdotf_kgs(i,:)*DeltaTime(i,:);
  W_N(i+1,:) = W_N(i,:) - (del_Wf_kg(i,:)*g);

end

Treq_N(1,1) = Treq_N(2,1);
mdotf_kgs(1,1) = mdotf_kgs(1,1);
mdotf_kgs(end+1) = mdotf_kgs(end);
del_Wf_kg(end+1) = del_Wf_kg(end);

run run_C2_plots.m

%% C.3


for i = 1:height(H)-1

[SAR(i,:), SE(i,:)] = getSAR_SE(W_N(i,:), H(i,:), Mach(i,:), FlightPathAngle_Deg(i,:));

mpg(i,:) = (((SAR(i,:)/1000)*Km2Mile)/rho_fuel)*1000*L2gallon;

pmpg(i,:) = mpg(i,:)*PAX;

end

SAR(end,1) = SAR(end-1,1);
pmpg(end,1) = pmpg(end-1,1);


figure(4);
hold on;
plot(Time, SAR, 'k', 'LineWidth',1.0); 
plot(Time, pmpg, 'b--', 'LineWidth',1.0);
xlim([0 300]);
ylim([0 200]);
legend('Specific Air Range {\it SAR (m/kg)} ', 'Passenger MPG, {\it p - mpg} (miles/gal) ')
xlabel('Flight Time Elapsed {\it t_m_i_n_u_t_e_s}');
ylabel('Specific Air Range / Passenger Miles per Gallon')
title('Specific Air Range and Passenger MPG along Flight Path');
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'off')
axis(gca,'square')


%embeded for loop calculates SAR and SE for a varying mach number at a set
%altitude, for all given altitudes and overwrites vectors into collums of
%corresponding SAR and SE for each mach at a given altitude
for i = 1:height(H_SAR)
    for j = 1:height(Mach_SAR)

        [SAR_3Alts(j,i), SE_3Alts(j,i)] = getSAR_SE(195000*9.81, H_SAR(i), Mach_SAR(j),0);

    end
end


SAR_3Weights = zeros(height(Mach_SAR), height(H_SAR));
SE_3Weights = zeros(height(Mach_SAR), height(H_SAR));


%embeded for loop calculates SAR and SE for a varying mach number at a set
%altitude, for all given altitudes and overwrites vectors into collums of
%corresponding SAR and SE for each mach at a given altitude
for i = 1:height(W_SAR)
    for j = 1:height(Mach_SAR)

        [SAR_3Weights(j,i), SE_3Weights(j,i)] = getSAR_SE(W_SAR(i), 360*30.48, Mach_SAR(j),0);

    end
end


figure(5)
plot(Mach_SAR, SAR_3Alts, 'LineWidth',1.5);
legend('SAR @ Sea Level', 'SAR @ 10000ft', 'SAR @ 36000ft')
xlabel('Mach')
ylabel('SAR using getSAR-SE Function')
title('SAR as function of Mach, 3 flight levels, constant W = 195 kgf')
ylim([0 200])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'off')
axis(gca,'square')



figure(6)
plot(Mach_SAR, SAR_3Weights, 'LineWidth',1.5);
legend('SAR @ 215000kgf', 'SAR @ 195000kgf', 'SAR @ 175000kgf')
xlabel('Mach')
ylabel('SAR using getSAR-SE Function')
title('SAR as function of Mach, 1 Flight Level, Varying Weight')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'off')
axis(gca,'square')

figure(7)
plot(Mach_SAR, SE_3Alts, 'LineWidth',1.5);
legend('SE @ Sea Level', 'SE @ 10000ft', 'SE @ 36000ft')
xlabel('Mach')
ylabel('SE using getSAR-SE Function')
title('SE as function of Mach, 3 Flight Levels, Constant Weight')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'off')
axis(gca,'square')



figure(8)
plot(Mach_SAR, SE_3Weights, 'LineWidth',1.5);
legend('SE @ 215000kgf', 'SE @ 195000kgf', 'SE @ 175000kgf')
xlabel('Mach')
ylabel('SE using getSAR-SE Function')
title('SE as function of Mach, 1 Flight Level, Varying Weight')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'off')
axis(gca,'square')



%% C4
for i = 1:height(BFR)
    [Range_VH(i,:)] = get_RangeConstMachConstCL(BFR(i,:));
    [Range_VM(i,:)] = get_RangeConstAltConstCL(BFR(i,:));
    [Range_VCL(i,:)] = get_RangeConstAltConstMach(BFR(i,:));
      
end

figure(9)
hold on
plot(BFR,Range_VH/1000, 'k',  'LineWidth', 1.5);
plot(BFR,Range_VM/1000, '--b', 'LineWidth', 1.5);
plot(BFR,Range_VCL/1000, ':', 'LineWidth', 2, 'color', '#EDB120');
legend('Range Under Variable Altitude Program','Range Under Variable Mach Program', 'Range Under Variable CL Program')
legend('Location', 'northwest')
xlabel('Block Fuel Ratio \zeta')
ylabel('Range, (Km)')
title('Aircraft Rang as a function of Block Fuel Ratio (W = MTOW)')
yticks([2500 5000 7500 10000 12500 15000])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
axis(gca,'square')

Zeta_MPR = (Wmtow_N - OEW  - MSP)/Wmtow_N;
Zeta_MER = (MFC)/Wmtow_N;
Zeta_MFR = (MFC)/(MFC + OEW);

[~, ~, rhos_360, a_360] = AtmosProp(360*30.48);
V_BR_Cruise = Mach_crit * a_360;


CL_MPR = (2/rhos_360) * (Wmtow_N/S) * (1/V_BR_Cruise^2);
CL_MER =(2/rhos_360) * (Wmtow_N/S) * (1/V_BR_Cruise^2);
CL_MFR = (2/rhos_360) * ((MFC + OEW)/S) * (1/V_BR_Cruise^2);

CD_MPR = CD0 + K*CL_MPR^2;
CD_MER = CD0 + K*CL_MER^2;
CD_MFR = CD0 + K*CL_MFR^2;


[~, T_360, ~, a_360] = AtmosProp(360*30.48);

[C_tp_360] = get_TSFC(Mach_crit, T_360/T0);

MPR = (V_BR_Cruise/(C_tp_360*g))*(CL_MPR/CD_MPR)*log(1/(1-Zeta_MPR));
MER = (V_BR_Cruise/(C_tp_360*g))*(CL_MER/CD_MER)*log(1/(1-Zeta_MER));
MFR = (V_BR_Cruise/(C_tp_360*g))*(CL_MFR/CD_MFR)*log(1/(1-Zeta_MFR));


%================================= Task 4==============================
BFR_BR = linspace(0,Zeta_MFR,10000)';
Range_Array = linspace(0,MFR,10000)';
PAY = zeros(height(Range_Array),1);
TOW = zeros(height(Range_Array),1);
FUEL = zeros(height(Range_Array),1); 



for i = 1:height(Range_Array)
    if BFR_BR(i,:) <= Zeta_MPR
        PAY(i,:) = MSP + OEW;
        TOW(i,:) = (PAY(i,:))/(1-BFR_BR(i,:));
        FUEL(i,:) = TOW(i,:) - PAY(i,:);
    elseif BFR_BR(i,:) > Zeta_MPR && BFR_BR(i,:) <= Zeta_MER
        TOW(i,:) = Wmtow_N;
        FUEL(i,:) = TOW(i,:) * BFR_BR(i,:);
        PAY(i,:) = TOW(i,:) - FUEL(i,:);
    elseif BFR_BR(i,:) > Zeta_MER
        FUEL(i,:) = MFC;
        TOW(i,:) = FUEL(i,:)/BFR_BR(i,:);
        PAY(i,:) = TOW(i,:) - FUEL(i,:);
    end
end




figure (10)
hold on
plot(Range_Array/(10^3),TOW/g, 'k',  'LineWidth', 1.5, 'color', "#EDB120");
plot(Range_Array/(10^3),PAY/g,  'LineWidth', 1.5, 'color', "#D95319")
yline(OEW/g,  'LineWidth', 1.5, 'color', "#0072BD")
yline(Wmtow_N/g,  'LineWidth', 2, 'LineStyle',':','color','b')
legend('OEW + PAY + FUEL = TOW', 'OEW + PAY', 'OEW', 'MTOW')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'on')
axis(gca,'square')
xline(MER/10^3)
xline(MPR/10^3)
xline(MFR/10^3)
%% Sheet 5

for i = 1:height(Trim_Altitudes)
    [VNE_ms(i,:), VME_ms(i,:), VstallTO_ms(i,:), VstallLD_ms(i,:)] = get_VNEspeeds(Trim_Altitudes(i,:));
    [VlowTA_ms(i,:),VhighTA_ms(i,:)] = get_MaxThrustSpeeds(Trim_Altitudes(i,:));
end


figure (11)

set(gcf,'Color',[1,1,1])
hold on
plot(VNE_ms,Trim_Altitudes/30.48,'g-','LineWidth',2)
plot(VME_ms,Trim_Altitudes/30.48,'k-','LineWidth',2)
plot(VstallTO_ms,Trim_Altitudes/30.48, 'r-','LineWidth',2)
plot(VstallLD_ms,Trim_Altitudes/30.48,'r--','LineWidth',2)
plot(VlowTA_ms(1:12300,1),Trim_Altitudes(1:12300,1)/30.48, 'LineWidth', 2, 'Color', 'c')
plot(VhighTA_ms(1:12300,1),Trim_Altitudes(1:12300,1)/30.48, 'LineWidth', 2, 'Color', 'b')
legendnames = [{'Never Exceed (Low Alt.)'},{'Never Exceed (High Alt.)'},{'Stall, T/O Config'},{'Stall, L/D Config'},{'Max. Thrust (Low Speed)'},{'Max. Thrust (High Speed) '}]; 
legend(gca,legendnames,'Interpreter','latex') % instruct Matlab to read late
set(legend,'Location','NorthWest') % set legend location
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.025 0.025])
set(gca,'LineWidth',1)
grid(gca,'off')  
box(gca,'on')
axis(gca,'square')
ylim([0 600]);
title("A330 - 300 Flight Envelope")
xlabel('True Airspeed, m/s'); % SET X-AXIS LABEL
ylabel('Altitude, FL');








