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
plot(Time,Relative_Pressure, ':', 'LineWidth',1.5);
plot(Time, Relative_Density, '-.', 'LineWidth',1.0);
plot(Time,Relative_Temperature, 'b', 'LineWidth',1.0);
plot(Time, Relative_Measure_Temp, 'c--', 'LineWidth',1.0);
plot(Time, Mach, 'k', 'LineWidth',1.0);
xlim([0 310]);
ylim([0 1.05]);
grid on;
legend('ISA Relative Pressure \delta', 'ISA Relative Density \sigma', ...
    'ISA Relative Temperature \theta', 'Relative Measure Temperature \theta_m_e_a_s', ...
    'Mach Number {\it M}');
legend('Location', 'southoutside');
xlabel('Flight Time Elapsed {\it t_m_i_n_u_t_e_s}');
ylabel('Relative Property Value or Mach Number /\theta /\delta /\sigma /{\it M} (-)');
title('Variation of Relative Atmospheric Properties and Mach Number along Flight Path');


figure(2);
hold on;
grid on;
plot(Time,True_AS_PredictedT, 'b', 'LineWidth',1.0);
plot(Time,True_AS_MeasuredT, 'k', 'LineWidth',1.0); 
xlim([0 310]);
legend('TAS Derived From Predicted Static Temperature, {\it TAS_I_S_A}', ...
    'TAS Derived From Measured Static Temperature, {\it TAS_m_e_a_s}');
legend('Location', 'south');
title('Comparison of True Airspeed, as calculated based on ISA-Derived Temperature ', ...
    'vs. aircraft-measured static temperature, along Flight Path');




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
grid on;
legend('Specific Air Range {\it SAR (m/kg)} ', 'Passenger MPG, {\it p - mpg} (miles/gal) ')
xlabel('Flight Time Elapsed {\it t_m_i_n_u_t_e_s}');
ylabel('Specific Air Range / Passenger Miles per Gallon')
title('Specific Air Range and Passenger MPG along Flight Path');




%% C.3 Task 3

%Create and array of altitudes weights and mach numbers


%assigning constants and variables to calculate from and converting them
%into the relevant units for the functions
Mach_crit = 0.85;
H_SAR = [0 100 360]'.*30.48;
W_SAR = [215 195 175]'.*1000.*g;


%pre allocating vectors to be overwritten in the for loops below
Mach_SAR = linspace(0.1, Mach_crit, 4537)';
SAR_3Alts = zeros(height(Mach_SAR), height(H_SAR));
SE_3Alts = zeros(height(Mach_SAR), height(H_SAR));


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
grid on 
ylim([0 200])



figure(6)
plot(Mach_SAR, SAR_3Weights, 'LineWidth',1.5);
legend('SAR @ 215000kgf', 'SAR @ 195000kgf', 'SAR @ 175000kgf')
xlabel('Mach')
ylabel('SAR using getSAR-SE Function')
title('SAR as function of Mach, 1 Flight Level, Varying Weight')
grid on


figure(7)
plot(Mach_SAR, SE_3Alts, 'LineWidth',1.5);
legend('SE @ Sea Level', 'SE @ 10000ft', 'SE @ 36000ft')
xlabel('Mach')
ylabel('SE using getSAR-SE Function')
title('SE as function of Mach, 3 Flight Levels, Constant Weight')
grid on 



figure(8)
plot(Mach_SAR, SE_3Weights, 'LineWidth',1.5);
legend('SE @ 215000kgf', 'SE @ 195000kgf', 'SE @ 175000kgf')
xlabel('Mach')
ylabel('SE using getSAR-SE Function')
title('SE as function of Mach, 1 Flight Level, Varying Weight')
grid on 


%% C4

% constants given in the task sheet

Norm_Op_Empty_W = 0.545;   %OEW/MTOW
Norm_Max_Struct_Pay = 0.233;   %MSP/MTOW
Norm_Max_Fuel_Cap = 0.358;   %MFC/MTOW
Fres = 0.1;   %t all times

BFR = linspace(0,Norm_Max_Fuel_Cap,4537)';

TOW_VarBFR = zeros(height(BFR),1);
RangeVarMach_km = zeros(height(BFR),1);
RangeVarCL_km = zeros(height(BFR),1);
RangeVarH_km = zeros(height(BFR),1);

Range_VH = zeros(height(BFR),1);
Range_VM = zeros(height(BFR),1);
Range_VCL = zeros(height(BFR),1);

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
grid on
legend('Location', 'northwest')
xlabel('Block Fuel Ratio \zeta')
ylabel('Range, (Km)')
title('Aircraft Rang as a function of Block Fuel Ratio (W = MTOW)')
yticks([2500 5000 7500 10000 12500 15000])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
axis(gca,'square')




