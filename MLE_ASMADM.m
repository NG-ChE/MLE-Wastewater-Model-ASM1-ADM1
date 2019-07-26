tic
%% TO-DO

% **Figure out why ammonia is not nitrifying**

% Implement FOG for AD
    % No chemical characeristics given, not sure what to do or if even
    % implement
% Implement primary bypass
% Implement dewatering(after AD)/dryer(after AD)
    % Thickener implemented
% Implement denitrification filter
% Implement proper secondary clarifier model
    % Primary clarifier is done but is not removing enough TSS than
    % whats actually happening at the plant

% Create PID for oxygen transfer in aeration zone
% Create GUI for kinetics/ recycle ratios/ simulation time/ flow variables

%% Issues
% Digester going crazy after ~300 days, believed to be from huge TSS load
% in stream 15 at day 234
% Carbon balance on AD is wrong and needs to be adjusted for new ASM/ADM
% conversion implementation
    % Stick to ADM variables to calculate carbon flow
% ***** CHECK UNITS *****

clc

%% Recycle ratios for MLE system
Var.Rir = 1.88; % Specify internal recycle ratio MLR/PC_effluent 
Var.Rr = 0.71; % Specify return activated sludge recycle ratio RAS/PC_effluent 
Var.fsc = 0.974; % Secondary Clarifier underflow separation (has to be less than or equal to 1)
if Var.fsc > 1
    Var.fsc = 1;
elseif Var.fsc < 0
    Var.fsc = 0;
else
end
Var.ft = 0.1976; % Flow fraction of concentrated TSS stream with respect to inflow

%% biological parameters and volumes
% Check parameter values with other reported values
% Make them a function of temperature
% Correct for units
Var.param = [0.67 0.24 0.08 0.08 ...
0.06 4 10 0.2 ...
0.5 0.8 0.8 1 ...
0.4 0.3 0.05 0.05 ...
3 0.1 0.8 100 ...
10 ...
416.83 ... %1782.92895 ... 416.83 ... % Primary clarifier volume m3 (470,000 gal) taken from B&V
5522.15 ... %999.348711 ... 5522.15 ... % Anoxic tank volume m3 (264,000 gal) taken from B&V
36339.95 .... %999.348711 ... 36339.95 .... % Aeration tank volume m3 (264,000 gal) taken from B&V
4803.84 ... %19229.891863 ... 4803.84 ... % Secondary clarifier volume m3 (5,080,000 gal) taken from B&V
0 ...
6056.65888]'; % Digester total volume in m3 (taken as 1.6 MG)

%% AD parameters
% ADJUST VALUES TO PLANT SPECIFIC DIMENSIONS
l = IndataADM1_v2;

%% Import data for varying influent flow

[status,sheets] = xlsfinfo('simuPlantData.xlsx');
sumData = [];
for s = 1:numel(sheets)
    [data,titles] = xlsread('simuPlantData.xlsx',s);
    sumData = [sumData;data];
end
% Alkalinity set to 316 mg/L of calcium carbonate
% Value received from Albert 
fixData = sumData;
fixData(:,7) = 316;
InfluentData = fillmissing(fixData,'movmedian',100); 
% plot(sumData(:,1),InfluentData(:,2:7),'r.-',sumData(:,1),sumData(:,2:7),'b.-') 
% legend('Filled Missing Data','Original Data')

%% Convert typical units coming into plant to ASM1 variables
% Concentration for each component taken as an average for initial
% condition of influent, however, it gets overwritten in ODE
% Units g/m3 = mg/L
infChar.TSS = InfluentData(:,2); % mg/L
% If no VSS data, assume VSS/TSS ratio of 0.69 or 0.75
infChar.tvRatio = 0.69;
infChar.VSS = infChar.TSS.*infChar.tvRatio;
infChar.FSS = infChar.TSS - infChar.VSS; % fixed suspended solids
% USING CBOD FROM DATA
infChar.cBOD = InfluentData(:,3); % mg/L
infChar.TKN = InfluentData(:,4); % mg-N/L
infChar.NH3 = InfluentData(:,5); % mg-N/L
infChar.NO3 = InfluentData(:,6); % mg-N/L
% Convert cBOD to total BOD5
% Source from https://www.wastewaterelearning.com/elearning/pluginfile.php/69/mod_resource/content/4/What%20is%20the%20Difference%20in%20BOD5-CBOD.pdf
infChar.BOD = 1.5*infChar.cBOD + 4.6*infChar.TKN;
% Alkalinity not given, will have Thursday from Albert
infChar.ALK = InfluentData(:,7); % mg/L
% % Conversion to ASM1 variables using Biological Wastewater Treatment by
% % Grady, Daigger, Love and Filipe, from Chapter 9.6
infChar.CODt = 2.1.*infChar.BOD; % mgCOD/L, Converts BOD5 to total COD, if not available
infChar.CODbo = 1.71.*infChar.BOD; % mgCOD/L, biodegradable COD
infChar.CODio = infChar.CODt - infChar.CODbo; % mgCOD/L, inert COD
% ASM.Xio = 0.375*1.5.*infChar.VSS; % mgCOD/L, particulate inert COD
% ASM.Sio = infChar.CODio - ASM.Xio; % mgCOD/L, soluble inert COD
% ASM.f_readily = 0.43; % fraction of biodegradable COD that is readily biodegradable
% ASM.Sso = infChar.CODbo.*ASM.f_readily; % mgCOD/L, readily biodegradable substrate
% ASM.Xso = infChar.CODbo - ASM.Sso; % mgCOD/L, slowly biodegradable substrate
% ASM.ONtotal = infChar.TKN - infChar.NH3; % mg-N/L, total organic nitrogen
% ASM.Snio = 1.5; % mg-N/L, soluble inert organic nitrogen
% ASM.in_xd = 0.06; % mass of nitrogen per mass of COD in biomass
% ASM.Xnio = ASM.in_xd.*ASM.Xio; % mg-N/L, 
% ASM.Snso_Xnso = ASM.ONtotal - ASM.Snio - ASM.Xnio; % mg-N/L, biodegradable organic nitrogen
% ASM.Sndo = ASM.Snso_Xnso.*(ASM.Sso./(ASM.Sso + ASM.Xso)); % mg-N/L, soluble biodegradable nitrogen
% ASM.Xndo = ASM.Snso_Xnso - ASM.Sndo; % mg-N/L, particulate biodegradable nitrogen
ASM.Xbho = repelem(0.000000001,length(InfluentData))'; % mgCOD/L, heterotrophic active biomass -> cant be zero, but very close to it
ASM.Xbao = repelem(0.000000001,length(InfluentData))'; % mgCOD/L, autrophic active biomass -> cant be zero, but very close to it
ASM.Soo = repelem(0,length(InfluentData))'; % mgO2/L, oxygen concentration
ASM.Xpo = repelem(0,length(InfluentData))'; % mgCOD/L, biomass debris 
ASM.Salko = infChar.ALK./100; % mol/L, alkalinity
ASM.Snho = infChar.NH3; % Initial ammonia
ASM.Snoo = infChar.NO3; % Initial nitrite/nitrate
% alt method (pL github)
ASM.Sio = 0.13.*infChar.CODt;
ASM.Xio = infChar.CODio - ASM.Sio;
ASM.Xso = 1.6.*infChar.VSS - ASM.Xio;
ASM.Sso = infChar.CODbo - ASM.Xso;
ASM.nb_TKN = infChar.TKN.*0.03;
ASM.sol_bio_orgN_ratio = ASM.Sso./(ASM.Sso + ASM.Xso);
ASM.Sndo = (infChar.TKN - ASM.Snho - ASM.nb_TKN).*(ASM.sol_bio_orgN_ratio);
ASM.Xndo = (infChar.TKN - ASM.Snho - ASM.nb_TKN).*(1 - ASM.sol_bio_orgN_ratio);

%% Intial conditions for system
MLE_influent = [ASM.Sio ASM.Sso ASM.Xio ASM.Xso ASM.Xbho ASM.Xbao ASM.Xpo ...
    ASM.Soo ASM.Snoo ASM.Snho ASM.Sndo ASM.Xndo ASM.Salko];
MLE_int = mean(MLE_influent); % Create vector of initial values for the MLE system based on total average
AD_int = [0.009;... % S_su
0.0009;...          % S_aa
0.0009;...          % S_fa 
0.0009;...          % S_va
0.0009;...          % S_bu
0.0009;...          % S_pro
0.0009;...          % S_ac
2.3594e-9;...       % S_h2
2.3594e-6;...       % S_ch4
0.039;...           % S_IC
0.13023;...         % S_IN
0.009;...           % S_I
0.30870;...         % X_c
0.02795;...         % X_ch
0.10260;...         % X_pr
0.02948; ...        % X_li
0.42016;...         % X_su
1.17917;...         % X_aa
0.24303;...         % X_fa
0.43192;...         % X_c4
0.13730;...         % X_pro
0.76056;...         % X_ac
0.31702;...         % X_h2
25.61739;...        % X_I
0.04;...            % S_cat
0.02;...            % S_an
0.0116;...          % S_vam
0.01322;...         % S_bum
0.01574;...         % S_prom
0.19724;...         % S_acm
0.14278;...         % S_hco3m
0.00409;...         % S_nh3
1.023e-5;...        % h2
1.62125;...         % ch4
0.01411]';          % co2
% "start up"
sys_int = [MLE_int AD_int]; % Combine initial conditions
x = sys_int(:)*ones(1,17); % Format to an array of [components,streams]
% steady state simulation
%x = initialValuesAll(MLE_int);

% Manipulate data for ODE input
Var.Qt = InfluentData(:,1); % Time, days
Var.Ct = Var.Qt; % Time, days
Var.Qflow = InfluentData(:,11); % Flow at time, Qt
Var.C = MLE_influent; % Conc at time, Ct

% Correct for duplicates in data by adding a small increment to each element
k = 1;
while k < (length(Var.Qflow)+1)
    if k == 1
       Var.Qt(k) = 1; % Set time to exactly 1 to avoid interpolation errors on initial points for Qflow and C in ODE
       Var.Ct(k) = 1;
    else
    Var.Qt(k) = Var.Qt(k) + (rand*rand)*1E-5 + 1; 
    Var.Ct(k) = Var.Qt(k);
    end
Var.Qflow(k) = Var.Qflow(k) + (rand*rand)*1E-3;
Var.C(k,:) = Var.C(k,:) + (rand*rand)*1E-7;
k = k + 1;
end

%% Simulation time span (days)
% Increase timespan interval for more data points in result section, can
% signicantly increase time if interval is very small
%t = 1:0.1:fixData(end,1);
t = 1:0.1:50;

% Sample rate
    % Decrease sample rate for better DO control, but ODE takes longer
sp = 1/5;
Var.tsample = t + sp*t;
Var.timespan = t;

%% ODE
odefunc = @(t,x) MLE(t,x,Var,l);
opts = odeset('MStateDependence','JPattern');
ODE_sol = ode15s(odefunc,t,x,opts);
toc
ODEToc = toc;

%% Results
tic
disp('Manipulating large data set, be patient')
% Get the data for each stream from the column of Concentration, then
% manipulate each array to be used for graphing
% Set first 13 rows to Concentration (ASM1 variables) 
% Set the next 35 rows to Conc_AD (ADM1 variables)
% Reloop
arrayManip = Array.mleArray;
int_len = length(arrayManip);
compASM = 13; % Number of components in ASM1
compADM = 35; % Number of components in ADM1
loop_len = int_len/(compADM + compASM);
Concentration = [];
Conc_AD = [];
i = 1;
na_int = 1; % Start of ASM1 variables
na_end = na_int + 12; % End of ASM1 variables
oa_int = na_end + 1; % Start of ADM1 variables
oa_end = na_end + 35; % End of ADM1 variables
while i < (loop_len + 1)
    Concentration = [Concentration;arrayManip(na_int:na_end,:)];
    Conc_AD = [Conc_AD;arrayManip(oa_int:oa_end,15:17)];
    na_int = oa_end + 1;
    na_end = na_int + 12;
    oa_int = na_end + 1;
    oa_end = na_end + 35;
    i = i + 1;
end

% Set variables
time = Array.tArray;
[~,col] = size(Concentration);

%% Plotting separate results
%% MLE + AD Results
% Create ASM Plant Stream Structure (Based on ASM1 Variables)
% Stream 14 is excluded as it is only a gas stream of 3 AD Components
ASMstream.One = reshape(Concentration(:,1),[compASM,length(time)])';
ASMstream.Two = reshape(Concentration(:,2),[compASM,length(time)])';
ASMstream.Three = reshape(Concentration(:,3),[compASM,length(time)])';
ASMstream.Four = reshape(Concentration(:,4),[compASM,length(time)])';
ASMstream.Five = reshape(Concentration(:,5),[compASM,length(time)])';
ASMstream.Six = reshape(Concentration(:,6),[compASM,length(time)])';
ASMstream.Seven = reshape(Concentration(:,7),[compASM,length(time)])';
ASMstream.Eight = reshape(Concentration(:,8),[compASM,length(time)])';
ASMstream.Nine = reshape(Concentration(:,9),[compASM,length(time)])';
ASMstream.Ten = reshape(Concentration(:,10),[compASM,length(time)])';
ASMstream.Eleven = reshape(Concentration(:,11),[compASM,length(time)])';
ASMstream.Twelve = reshape(Concentration(:,12),[compASM,length(time)])';
ASMstream.Thirteen = reshape(Concentration(:,13),[compASM,length(time)])';
ASMstream.Fourteen = reshape(Concentration(:,14),[compASM,length(time)])';
ASMstream.Fifteen = reshape(Concentration(:,15),[compASM,length(time)])';
% ASM1 variables that have been converted from ADM1
ASMstream.Seventeen = reshape(Concentration(:,17),[compASM,length(time)])';

figure(1)
subplot(4,5,1)
plot(time,ASMstream.One)
title('Plant Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,2)
plot(time,ASMstream.Two)
title('Primary Clarifier Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,3)
plot(time,ASMstream.Three)
title('Primary Clarifier WAS')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,4)
plot(time,ASMstream.Four)
title('Mixing Point Before Anoxic Tank')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,5)
plot(time,ASMstream.Five)
title('Anoxic Tank Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,6)
plot(time,ASMstream.Six)
title('Aeration Tank Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,7)
plot(time,ASMstream.Seven)
title('Aeration Tank Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,8)
plot(time,ASMstream.Eight)
title('Internal Recycle')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,9)
plot(time,ASMstream.Nine)
title('Secondary Clarifier Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,10)
plot(time,ASMstream.Ten)
title('Secondary Clarifier Underflow')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,11)
plot(time,ASMstream.Eleven)
title('WAS Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,12)
plot(time,ASMstream.Twelve)
title('RAS Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,13)
plot(time,ASMstream.Thirteen)
title('Mixing of WAS and PS')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,14)
plot(time,ASMstream.Fourteen)
title('Thickener Centrate Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,15)
plot(time,ASMstream.Fifteen)
title('Thickener Sludge Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,16)
plot(time,ASMstream.Seventeen)
title('Digester Sludge Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,17)
plot(time,Array.Qarray)
title('Plant Influent Flowrate')
ylabel('Volumetric Flow, m3/day')
xlabel('Time, days')
subplot(4,5,18)
plot(time,ASMstream.Nine(:,10))
title('Plant Effluent Ammonia')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,19)
plot(time,ASMstream.Nine(:,9))
title('Plant Effluent Nitrate/Nitrite')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,20)
plot(time,ASMstream.Nine(:,8))
title('Plant Effluent DO')
ylabel('Concentration, mg/L')
xlabel('Time, days')

%% Anaerobic Digester Results
% Create Array for AD Components
ADMstream.One = reshape(Conc_AD(:,1),[compADM,length(time)])';
ADMstream.Two = reshape(Conc_AD(:,2),[compADM,length(time)])';
ADMstream.Three = reshape(Conc_AD(:,3),[compADM,length(time)])';

% Only components in ADMstream.Two should be the Gas components
% Delete dummy components in the array
ADMstream.Two(:,1:32) = [];
ADMstream.S_gas_h2 = ADMstream.Two(:,1);
ADMstream.S_gas_ch4 = ADMstream.Two(:,2);
ADMstream.S_gas_co2 = ADMstream.Two(:,3);

% Determine partial pressures for each gas component, converting CH4 and H2 from kgCOD to
% kmole C (S_gas_h2/S_gas_ch4)
ADMstream.P_gas_h2 = ADMstream.S_gas_h2*l.R*l.T_op/16;
ADMstream.P_gas_ch4 = ADMstream.S_gas_ch4*l.R*l.T_op/64;
ADMstream.P_gas_co2 = ADMstream.S_gas_co2*l.R*l.T_op;

% Look into changing gas flow equation
%q_gas = ((l.R*l.T_op)/(P_atm - l.p_gas_h2o))*l.V_liq*((rho_T_8/16) + (rho_T_9/64) + rho_T_10);
% Gas flow
ADMstream.P_gas = ADMstream.P_gas_h2 + ADMstream.P_gas_ch4 + ...
    ADMstream.P_gas_co2 + l.p_gas_h2o; % Total gas pressure
ADMstream.q_gas = l.k_p*(ADMstream.P_gas - l.P_atm).* ...
    ADMstream.P_gas/l.P_atm; % Total gas flow

% Calculate gas component gas flow
ADMstream.q_gas_h2 = ADMstream.P_gas_h2./ADMstream.P_gas.*ADMstream.q_gas;
ADMstream.q_gas_ch4 = ADMstream.P_gas_ch4./ADMstream.P_gas.*ADMstream.q_gas;
ADMstream.q_gas_co2 = ADMstream.P_gas_co2./ADMstream.P_gas.*ADMstream.q_gas;
ADMstream.q_gas_h2o = l.p_gas_h2o./ADMstream.P_gas.*ADMstream.q_gas;

% Mole fraction of gas components
ADMstream.PyH2 = ADMstream.P_gas_h2./ADMstream.P_gas;
ADMstream.PyCO2 = ADMstream.P_gas_co2./ADMstream.P_gas;
ADMstream.PyCH4 = ADMstream.P_gas_ch4./ADMstream.P_gas;

figure(2)
subplot(2,2,1);
plot(time,ADMstream.One);
title('Anaerobic Digester Influent')
ylabel('Concentration, kg/m3')
xlabel('Time, days')
subplot(2,2,2);
plot(time,ADMstream.Two);
title('Anaerobic Digester Gas Effluent')
ylabel('Concentration, kg/m3')
xlabel('Time, days')
legend('H2','CH4','CO2')
subplot(2,2,3);
plot(time,ADMstream.Three);
title('Anaerobic Digester Liquid Effluent')
ylabel('Concentration, kg/m3')
xlabel('Time, days')
subplot(2,2,4);
plot(time,ADMstream.PyCH4,time,ADMstream.PyCO2,time,ADMstream.PyH2);
title('Anaerobic Digester Partial Pressure')
ylabel('Partial Pressure')
xlabel('Time, days')
legend('CH4','CO2','H2')

%% Carbon (grams) in system
% Carbon content fractions
% Fractions taken from B&V study
CF.C_S_I = 0.32; % C content of soluble inert material (Si) [gC/gCOD]
CF.C_I_P = 0.32; % C content of inert particulate material (Xi) [gC/gCOD]
CF.C_SB_S = 0.32; % C content of slowly biodegradable substrate (Xs) [gC/gCOD]
CF.C_S_S = 0.32; % C content of soluble substrate (Ss) [gC/gCOD]
CF.C_C = 0.32; % C content of colloidal material [gC/gCOD]
CF.C_A = 0.375; % C content of acetate [gC/gCOD]
CF.C_P = 0.321; % C content of propionate [gC/gCOD]
CF.C_M = 0.25; % C content of methanol [gC/gCOD]
CF.C_S_M = 0.188; % C content of soluble methane [gC/gCOD]
CF.C_CO2 =  0.273; % C content of carbon dioxide [gC/gCO2]
CF.C_CC = 0.12; % C content of calcium carbonate (Salk) [gC/gCaCO3]
CF.C_HB = 0.366; % C content of heterotrophic biomass (Xbh) [gC/gCOD]
CF.C_AOB = 0.366; % C content of ammonia-oxidizing biomass [gC/gCOD]
CF.C_NOB = 0.366; % C content of nitrite-oxidizing biomass [gC/gCOD]
% AOB/NOB are autotrophic biomass
CF.C_AB = 0.366; % C content of autotrophic biomass (Xba) [gC/gCOD]
CF.C_UB = 0.366; % C content of unbiodegradable residue (Xp) [gC/gCOD]

% Create arrays for carbon
% Each carbon component relates to a component in the ASM1 Concentration
% data set, unless otherwise specified, which will later be multiplied by the carbon content fraction
loop.CSI = 1; Carbon.CSI = [];
loop.CSS = 2; Carbon.CSS = [];
loop.CIP = 3; Carbon.CIP = [];
loop.CSBS = 4; Carbon.CSBS = [];
loop.CHB = 5; Carbon.CHB = [];
loop.CAB = 6; Carbon.CAB = [];
loop.CUB = 7; Carbon.CUB = [];
loop.CCC = 13; Carbon.CCC = [];
% Non-ASM1 variables
loop.CSM = 34; Carbon.CSM = [];
loop.CCO2 = 35; Carbon.CCO2 = [];
p = 1;
while p < (length(time) + 1)
    Carbon.CSI = [Carbon.CSI;Concentration(loop.CSI,:)];
    Carbon.CSS = [Carbon.CSS;Concentration(loop.CSS,:)];
    Carbon.CIP = [Carbon.CIP;Concentration(loop.CIP,:)];
    Carbon.CSBS = [Carbon.CSBS;Concentration(loop.CSBS,:)];
    Carbon.CHB= [Carbon.CHB;Concentration(loop.CHB,:)];
    Carbon.CAB = [Carbon.CAB;Concentration(loop.CAB,:)];
    Carbon.CUB = [Carbon.CUB;Concentration(loop.CUB,:)];
    Carbon.CCC = [Carbon.CCC;Concentration(loop.CCC,:)];
    Carbon.CSM = [Carbon.CSM;Conc_AD(loop.CSM,:)];
    Carbon.CCO2 = [Carbon.CCO2;Conc_AD(loop.CCO2,:)];
    loop.CSI = loop.CCC + 1;
    loop.CSS = loop.CCC + 2;
    loop.CIP = loop.CCC + 3;
    loop.CSBS = loop.CCC + 4;
    loop.CHB = loop.CCC + 5;
    loop.CAB = loop.CCC + 6;
    loop.CUB = loop.CCC + 7;
    loop.CCC = loop.CCC + 13;
    loop.CSM = loop.CCO2 + 34;
    loop.CCO2 = loop.CCO2 + 35;
    p = p + 1;
end 
% Convert concentration to mass/time 
Carbon.CSIgrams = Carbon.CSI.*Array.Qarray;
Carbon.CSIgrams = Carbon.CSI.*Array.Qarray;
Carbon.CSSgrams = Carbon.CSI.*Array.Qarray;
Carbon.CIPgrams = Carbon.CSI.*Array.Qarray;
Carbon.CSBSgrams = Carbon.CSI.*Array.Qarray;
Carbon.CHBgrams = Carbon.CSI.*Array.Qarray;
Carbon.CABgrams = Carbon.CSI.*Array.Qarray;
Carbon.CUBgrams = Carbon.CSI.*Array.Qarray;
Carbon.CCCgrams = Carbon.CSI.*Array.Qarray;
Carbon.CSMkgram = Carbon.CSM.*ADMstream.q_gas.*1000; % Convert to grams
Carbon.CCO2kmoleC = Carbon.CCO2.*ADMstream.q_gas.*1000; % Convert to mol

% Determine total carbon mass flow rate [gC/day] in stream m
% Multiply carbon fraction by the corresponding carbon component
% Methanol still needs to be implemented
m = 1;
Cflow = ones(length(time),col); % Preallocate array
while m < (col + 1)
    Cflow(:,m) = Carbon.CSIgrams(:,m).*CF.C_S_I + Carbon.CIPgrams(:,m).*CF.C_I_P + Carbon.CSBSgrams(:,m).*CF.C_SB_S +... 
    Carbon.CSSgrams(:,m).*CF.C_S_S + Carbon.CCCgrams(:,m).*CF.C_CC+ Carbon.CHBgrams(:,m).*CF.C_HB + Carbon.CABgrams(:,m).*CF.C_AB +...
    Carbon.CUBgrams(:,m).*CF.C_UB;
    m = m + 1;
end
% Adjust for gas stream
% Convert from COD to grams and mol to grams
Cflow(:,16) = Carbon.CSMkgram(:,2).*CF.C_S_M + Carbon.CCO2kmoleC(:,2).*12.01; 
% Convert grams to lbs 
Cflow = 0.00220462.*Cflow;
% Mass balance check
PerError = [];
PerError(:,1) = 100.*abs(Cflow(:,1) - (Cflow(:,2) + Cflow(:,3)))./Cflow(:,1);
PerError(:,2) = 100.*abs(Cflow(:,4) - (Cflow(:,12) + Cflow(:,8)+ Cflow(:,2)))./Cflow(:,4);
PerError(:,3) = 100.*abs(Cflow(:,4) - Cflow(:,5))./Cflow(:,4);
PerError(:,4) = 100.*abs(Cflow(:,5) - Cflow(:,6))./Cflow(:,5);
PerError(:,5) = 100.*abs(Cflow(:,6) - (Cflow(:,8) + Cflow(:,7)))./Cflow(:,6);
PerError(:,6) = 100.*abs(Cflow(:,7) - (Cflow(:,9) + Cflow(:,10)))./Cflow(:,7);
PerError(:,7) = 100.*abs(Cflow(:,10) - (Cflow(:,12) + Cflow(:,11)))./Cflow(:,10);
PerError(:,8) = 100.*abs(Cflow(:,1) - (Cflow(:,11) + Cflow(:,9) + Cflow(:,3)))./Cflow(:,1);
PerError(:,9) = 100.*abs(Cflow(:,13) - (Cflow(:,11) + Cflow(:,3)))./Cflow(:,13);
PerError(:,10) = 100.*abs(Cflow(:,13) - Cflow(:,15) - Cflow(:,14))./Cflow(:,13);
% Ignore Carbon balance on AD for now
%PerError(:,11) = 100.*abs(C(:,15) - C(:,17) - C(:,16))./C(:,13);
% Plot percent error
figure(3)
plot(time,PerError)
title('Mass Balance Percent Error')
ylabel('Error, %')
xlabel('Time, days')

% Plot carbon flow for each stream separately
figure(4)
subplot(4,5,1)
plot(time,Cflow(:,1))
title('Plant Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,2)
plot(time,Cflow(:,2))
title('Primary Clarifier Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,3)
plot(time,Cflow(:,3))
title('Primary Clarifier WAS')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,4)
plot(time,Cflow(:,4))
title('Anoxic Tank Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,5)
plot(time,Cflow(:,5))
title('Anoxic Tank Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,6)
plot(time,Cflow(:,6))
title('Aeration Tank Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,7)
plot(time,Cflow(:,7))
title('SC Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,8)
plot(time,Cflow(:,8))
title('Internal Recycle')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,9)
plot(time,Cflow(:,9))
title('SC Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,10)
plot(time,Cflow(:,10))
title('SC Underflow')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,11)
plot(time,Cflow(:,11))
title('WAS Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,12)
plot(time,Cflow(:,12))
title('RAS Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
% Ignore carbon balance on AD for now
% subplot(4,5,13)
% plot(time,C(:,13))
% title('Mixing of WAS and PS')
% ylabel('Carbon Flow, lb/day')
% xlabel('Time, days')
% subplot(4,5,16)
% plot(time,C(:,16))
% title('Anaerobic Digester Gas Stream')
% ylabel('Carbon Flow, lb/day')
% xlabel('Time, days')
% subplot(4,5,17)
% plot(time,C(:,17))
% title('Anaerobic Digester Sludge Stream')
% ylabel('Carbon Flow, lb/day')
% xlabel('Time, days')

toc
resultsToc = toc;
totalTime = ODEToc + resultsToc;
fprintf('The total elapsed time is %.2f minutes',(totalTime/60))

function [Conc] = MLE(t,dCdt,Var,l)
persistent KLa
persistent Array
persistent ft
CompASM = 13;
CompADM = 35;
%% Dynamic flow
% Constant Plant flow - > testing average data -> using gal/min going into
% NT flow converted to m3/day
Qplant = interp1(Var.Qt,Var.Qflow,t); % Interpolate data set of volumetric flow at specified time

%% Solve flow balance
if t == 1
ft = Var.ft;
else
end
Q = zeros(1,17);
% Q(1) = Plant influent
% Q(2) = Primary clarifier effluent
% Q(3) = Primary clarifier waste sludge
% Q(4) = Anoxic influent
% Q(5) = Anoxic effluent
% Q(6) = Aeration Effluent
% Q(7) = Secondary clarifier influent
% Q(8) = Internal recycle (MLR)
% Q(9) = Secondary clarifier effluent
% Q(10) = Secondary clarifier underflow
% Q(11) = WAS
% Q(12) = RAS
% Q(13) = mixing of WAS and PC waste sludge sent to AD
% Q(14) = gas flow -> this is generated -> not in mass balance
% Q(15) = liquid flow
Q(1) = Qplant;
Q(3) = 189.27; 
Q(2) = Q(1) - Q(3);
Q(8) = Var.Rir*Q(2);
Q(12) = Var.Rr*Q(2);
Q(10) = Q(12)/Var.fsc;
Q(4) = Q(2) + Q(8) + Q(12);
Q(5) = Q(4);
Q(6) = Q(5);
Q(7) = Q(6) - Q(8);
Q(9) = Q(7) - Q(10);
Q(11) = Q(10) - Q(12);
Q(13) = Q(3) + Q(11);


% Fixes structure of intial values from vector to an array
dCdt = reshape(dCdt,[length(dCdt)/length(Q),length(Q)]);

%% Assigning Parameters
% Stoichiometric
Yh = Var.param(1); % gCOD(biomassformed)/gCOD(substrate removed)
Ya = Var.param(2); % gCOD(biomassformed)/gN(N oxidized)
ixb = Var.param(3); % gN/gCOD
fp = Var.param(4);
ixp = Var.param(5); % gN/gCOD

% Kinetic parameters
muh = Var.param(6); % 1/day
Ks = Var.param(7); % gCOD/m3
Koh = Var.param(8); % gN/m3
Kno = Var.param(9); % 0.5 gO2/m3
ng = Var.param(10); 
mua = Var.param(11); % 1/day
Knh = Var.param(12); % gN/m3
Koa = Var.param(13); % gO2/m3
bh = Var.param(14); % 1/day
ba = Var.param(15); % 1/day
ka = Var.param(16); % m3/gCOD/day
kh = Var.param(17); % 1/day
Kx = Var.param(18); 
nh = Var.param(19);

% O2 transference
So_sat = Var.param(21); % gO2/m3

% Equipment volumes
Vol1 = Var.param(22); % primary clarifier m3, taken L,W,B of one NT clarifier
Vol2 = Var.param(23); % anoxic m3, taken as total volume of NT anoxic tanks
Vol3 = Var.param(24); % aeration m3, taken as total volume of NT aeration tanks
Vol4 = Var.param(25); % secondary clarifier m3, taken as volume of a cylinder, and only one NT clarifier
%Vol5 = Var.param(26); % denit filters m3
Vol6 = Var.param(27); % AD volume, taken as two digesters in parallel as per B&V study


%% Component Identification
%Si  = dCdt(1,i); % Soluble inert organic matter
%Ss  = dCdt(2,i); % Readily biodegradable substrate
%Xi  = dCdt(3,i); % Particulate inert organic matter
%Xs  = dCdt(4,i); % Slowly biodegradable substrate
%Xbh = dCdt(5,i); % Active heterotrophic biomass
%Xba = dCdt(6,i); % Active autotrophic biomass
%Xp  = dCdt(7,i); % Particulate products arising from biomass decay
%So  = dCdt(8,i); % Oxygen
%Sno = dCdt(9,i); % Nitrate and nitrite nitrogen
%Snh = dCdt(10,i); % NH 4+ + NH 3 nitrogen
%Snd = dCdt(11,i); % Soluble biodegradable organic nitrogen
%Xnd = dCdt(12,i); % Particulate biodegradable organic nitrogen
%Salk = dCdt(13,i); % Alkalinity

%% Stoichiometric matrix
K =[0 -1/Yh  0  0          1   0      0   -(1-Yh)/Yh       0                     -ixb          0     0                -ixb/14 ;...
    0 -1/Yh  0  0          1   0      0   0                -(1-Yh)/(2.86*Yh)     -ixb          0     0                (1-Yh)/(14*2.86*Yh)-ixb/14 ;...
    0  0     0  0          0   1      0   -(4.56-Ya)/Ya    1/Ya                  -(ixb+1/Ya)   0     0                -ixb/14-1/(7*Ya) ;...
    0  0     0  (1-fp)    -1   0     fp   0                0                     0             0     (ixb-fp*ixp)     0 ;...
    0  0     0  (1-fp)     0  -1     fp   0                0                     0             0     (ixb-fp*ixp)     0 ;...
    0  0     0  0          0   0      0   0                0                     1            -1     0                1/14 ;...
    0  1     0  -1         0   0      0   0                0                     0             0     0                0;...
    0  0     0  0          0   0      0   0                0                     0             1    -1                0];

%% Reaction rates vector

% Anoxic tank dCdt(i,5)
theta1 = [muh*(dCdt(2,5)/(Ks+dCdt(2,5)))*(dCdt(8,5)/(Koh+dCdt(8,5)))*dCdt(5,5);...
    muh*(dCdt(2,5)/(Ks+dCdt(2,5)))*(Koh/(Koh+dCdt(8,5)))*(dCdt(9,5)/(Kno+dCdt(9,5)))*ng*dCdt(5,5);...
    mua*(dCdt(10,5)/(Knh+dCdt(10,5)))*(dCdt(8,5)/(Koa+dCdt(8,5)))*dCdt(6,5);...
    bh*dCdt(5,5);...
    ba*dCdt(6,5);...
    ka*dCdt(11,5)*dCdt(5,5);...
    kh*dCdt(4,5)/dCdt(5,5)/(Kx+dCdt(4,5)/dCdt(5,5))*((dCdt(8,5)/(Koh+dCdt(8,5)))+nh*Koh/(Koh+dCdt(8,5))*dCdt(9,5)/(Kno+dCdt(9,5)))*dCdt(5,5);...
    kh*dCdt(4,5)/dCdt(5,5)/(Kx+dCdt(4,5)/dCdt(5,5))*((dCdt(8,5)/(Koh+dCdt(8,5)))+nh*Koh/(Koh+dCdt(8,5))*dCdt(9,5)/(Kno+dCdt(9,5)))*dCdt(5,5)*dCdt(12,5)/dCdt(4,5)];

% Aeration tank dCdt(i,6)
theta2 = [muh*(dCdt(2,6)/(Ks+dCdt(2,6)))*(dCdt(8,6)/(Koh+dCdt(8,6)))*dCdt(5,6);...
    muh*(dCdt(2,6)/(Ks+dCdt(2,6)))*(Koh/(Koh+dCdt(8,6)))*(dCdt(9,6)/(Kno+dCdt(9,6)))*ng*dCdt(5,6);...
    mua*(dCdt(10,6)/(Knh+dCdt(10,6)))*(dCdt(8,6)/(Koa+dCdt(8,6)))*dCdt(6,6);...
    bh*dCdt(5,6);...
    ba*dCdt(6,6);...
    ka*dCdt(11,6)*dCdt(5,6);...
    kh*dCdt(4,6)/dCdt(5,6)/(Kx+dCdt(4,6)/dCdt(5,6))*((dCdt(8,6)/(Koh+dCdt(8,6)))+nh*Koh/(Koh+dCdt(8,6))*dCdt(9,6)/(Kno+dCdt(9,6)))*dCdt(5,6);...
    kh*dCdt(4,6)/dCdt(5,6)/(Kx+dCdt(4,6)/dCdt(5,6))*((dCdt(8,6)/(Koh+dCdt(8,6)))+nh*Koh/(Koh+dCdt(8,6))*dCdt(9,6)/(Kno+dCdt(9,6)))*dCdt(5,6)*dCdt(12,6)/dCdt(4,6)];

%% MLE mass balance
% Q(1) - Q(3) - Q(2); % primary clarifier
% Q(2) + Q(8) + Q(12) - Q(4); % inflow to anox 
% Q(5) - Q(4); % anox reactor
% Q(6) - Q(5); % aeration reactor
% Q(7) + Q(8) - Q(6); % recycle split
% Q(9) + Q(10) - Q(7); % secondary clarifier
% Q(12) + Q(11) - Q(10); % RAS/WAS split
% Q(3) + Q(11) - Q(13); % Mixing of WAS and primary sludge
% Q(14) + Q(15) - Q(13); % Thickener balance
% Q(15) - Q(17); % Digester balance

%% ODE for MLE system
% Concentration of component i in stream 1-12
% K(1,i)*theta1(1) + K(2,i)*theta1(2) + K(3,i)*theta1(3) + K(4,i)*theta1(4) + K(5,i)*theta1(5) + K(6,i)*theta1(6) + K(7,i)*theta1(7) + K(8,i)*theta1(8) 
i = 1;
Conc = zeros(length(dCdt),numel(Q));
while i < (CompASM + 1)
    % Model won't be used for my plant simulation
    %% Modeling primary clarifier - Otterpohl and Freund 1992
%     hrt = Vol1/Q(1) % Hydraulic residence time
%     n_COD = 2.7*log(hrt*hrt + 9)/100; % Removal efficiency
%     XCOD1 = dCdt(3,1) + dCdt(4,1) + dCdt(5,1) + dCdt(6,1) + dCdt(7,1); % Particulate COD in influent
%     CODin = dCdt(1,1) + dCdt(2,1) + XCOD1; % Total COD in influent
%     n_x = (n_COD*CODin)/XCOD1;
%     if n_x > 0.95
%        n_x = 0.95;
%     elseif n_x < 0.05
%        n_x = 0.05;
%     else
%        n_x = n_x;
%     end
%     n_x = (1-n_x)
    %% TSS removal
    n_x = 0.533; % Fraction of TSS left in effluent, taken as average from ST and NT from Appendix B GPS-X files from B&V
    % Determine which components are separated
    % Comment out the two lines below to run constant influent data
    Dyn_conc = interp1(Var.Ct,Var.C,t); % Interpolate data set of concentration at specified time
    dCdt(i,1) = Dyn_conc(:,i);
    
    if i < 3
        dCdt(i,2) = dCdt(i,1);
    elseif (2 < i) && (i < 8)
        dCdt(i,2) = n_x*dCdt(i,1)*Q(1)/Q(2);
    elseif (7 < i) && (i < 12)
        dCdt(i,2) = dCdt(i,1); 
    elseif (11 < i) && (i < 13)
        dCdt(i,2) = n_x*dCdt(i,1)*Q(1)/Q(2);
    else
        dCdt(i,2) = dCdt(i,1);
    end
    if i < 3
        dCdt(i,3) = dCdt(i,1);
    elseif (2 < i) && (i < 8)
        dCdt(i,3) = (1 - n_x)*dCdt(i,1)*Q(1)/Q(3);
    elseif (7 < i) && (i < 12)
        dCdt(i,3) = dCdt(i,1);
    elseif (11 < i) && (i < 13)
        dCdt(i,3) = (1 - n_x)*dCdt(i,1)*Q(1)/Q(3);
    else
        dCdt(i,3) = dCdt(i,1);
    end
    dCdt(i,2) = (dCdt(i,1)*Q(1) - dCdt(i,3)*Q(3))/Q(2); % Mass balance for flow into/out of Primary Clarifier
    
        %% Anox/Aer
    for intRec = 1:150
        if intRec == 1
        GC8(i,1) = dCdt(i,8);
        GC12(i,1) = dCdt(i,12);
        else
        end
    dCdt(i,4) = (Q(2)*dCdt(i,2) + Q(8)*GC8(i,intRec) + Q(12)*GC12(i,intRec))/Q(4); % mixing point
    if i == 8
    if t == 1
    KLa = Var.param(20); % Oxygen transfer coefficient
    else
    end
%     maxNH3 = 5;
%     if dCdt(10,6) > maxNH3
%         KLa = Var.param(20);
%     else
%         KLa = 0;
%     end
    else
    end
    Conc(i,5) = 1/Vol2*(Q(4)*dCdt(i,4) - Q(5)*dCdt(i,5)) + K(1,i)*theta1(1) + K(2,i)*theta1(2) + K(3,i)*theta1(3) + K(4,i)*theta1(4) + K(5,i)*theta1(5) + K(6,i)*theta1(6) + K(7,i)*theta1(7) + K(8,i)*theta1(8); % Anoxic balance
    Conc(i,6) = 1/Vol3*(Q(5)*dCdt(i,5) - Q(6)*dCdt(i,6)) + K(1,i)*theta2(1) + K(2,i)*theta2(2) + K(3,i)*theta2(3) + K(4,i)*theta2(4) + K(5,i)*theta2(5) + K(6,i)*theta2(6) + K(7,i)*theta2(7) + K(8,i)*theta2(8); % Aeration general balance
    if i == 8
        % Rough control of DO in aeration zone
%         if t == 1
%             KLa = Var.param(20); % Oxygen transfer coefficient
%         else
%         end
%         adjT = round(t*100)/100;
%         findT = find(adjT == Var.timespan);
%         if numel(findT) > 0
%             if dCdt(10,6) > 5
%                 KLa = 100;
%             else
%                 KLa = 0;
%             end
%         else
%         end
        Conc(8,6) = Conc(8,6) + KLa*(So_sat - dCdt(8,6)); % Effect of aeration on the Oxygen concentration
    else 
        Conc(i,6) = Conc(i,6);
    end
    % Calculate MCRT
    TSSas = 0.75*(dCdt(5,6) + dCdt(6,6))*Vol3; % Particulate COD in aeration tank
    TSSsc = 0.75*(dCdt(5,7) + dCdt(6,7))*Vol4; % Particulate COD in secondary clarifier
    SLe = TSSas*Q(9);
    SLw = TSSsc*Q(11);
    SRT = (TSSas + TSSsc)/(SLe + SLw);
    
    %% Recycle split, same concentration

    dCdt(i,8) = dCdt(i,6);
    dCdt(i,7) = dCdt(i,6);
    
    % Model not worth implementing, values taken from B&V will be used from
    % App. B from GPS-X pdf file
    %% Secondary clarifier model
    c_x = 0.0015; % Fraction of TSS left in effluent, taken as average from ST and NT of B&V App. B file
    XCOD7 = dCdt(3,7) + dCdt(4,7) + dCdt(5,7) + dCdt(6,7) + dCdt(7,7); % Particulate COD in influent
    if i < 3
        dCdt(i,9) = dCdt(i,7);
    elseif (2 < i) && (i < 8)
        dCdt(i,9) = c_x*dCdt(i,7)*Q(7)/Q(9);
    elseif (7 < i) && (i < 12)
        dCdt(i,9) = dCdt(i,7);
    elseif (11 < i) && (i < 13)
        dCdt(i,9) = c_x*dCdt(i,7)*Q(7)/Q(9);
    else
        dCdt(i,9) = dCdt(i,7);
    end
    if i < 3
        dCdt(i,10) = dCdt(i,7);
    elseif (2 < i) && (i < 8)
        dCdt(i,10) = (1 - c_x)*dCdt(i,7)*Q(7)/Q(10);
    elseif (7 < i) && (i < 12)
        dCdt(i,10) = dCdt(i,7);
    elseif (11 < i) && (i < 13)
        dCdt(i,10) = (1 - c_x)*dCdt(i,7)*Q(7)/Q(10);
    else
        dCdt(i,10) = dCdt(i,7);
    end
    
    dCdt(i,7) = (dCdt(i,9)*Q(9) + dCdt(i,10)*Q(10))/Q(7); % Mass balance for flow into/out of Primary Clarifier
    XCOD9 = dCdt(3,9) + dCdt(4,9) + dCdt(5,9) + dCdt(6,9) + dCdt(7,9); % Particulate COD in effluent
    XCOD10 = dCdt(3,10) + dCdt(4,10) + dCdt(5,10) + dCdt(6,10) + dCdt(7,10); % Particulate COD in waste stream
    
    %% Waste split

    dCdt(i,11) = dCdt(i,10);
    dCdt(i,12) = dCdt(i,10);
    
    if dCdt(i,8) <= 0 
        GC8(i,intRec) = 0;
        err1 = 0;
    else
        err1 = 100.*((GC8(i,intRec) - dCdt(i,8))./dCdt(i,8));
    end
    if dCdt(i,12) <= 0
        GC12(i,intRec) = 0;
        err2 = 0;
    else
        err2 = 100.*((GC12(i,intRec) - dCdt(i,12))./dCdt(i,12));
    end
    format long g
    table(GC8(i),dCdt(i,8))
    
    errT = abs(err1) + abs(err2);
    tol = 1;
    if errT > tol
        if err1 > 0
            q1 = 0.5;
            GC8(i,intRec + 1) = q1.*GC8(i,intRec) + (1 - q1)*dCdt(i,8);
        else
            q1 = -0.25;
            GC8(i,intRec + 1) = q1.*GC8(i,intRec) + (1 - q1)*dCdt(i,8);
        end
        if err2 > 0
            q2 = 0.5;
            GC12(i,intRec + 1) = q2.*GC12(i,intRec) + (1 - q2)*dCdt(i,12);
        else
            q2 = -0.25;
            GC12(i,intRec + 1) = q2.*GC12(i,intRec) + (1 - q2)*dCdt(i,12);
        end
    else
        break
    end
    end


%% Conversion from ASM1 to ADM1
    % Using paper: Benchmark Simulation Model No.2 (BSM2)
    %% Waste sludge mixing
    dCdt(i,13) = (dCdt(i,11)*Q(11) + dCdt(i,3)*Q(3))/Q(13);
    
    %% Thickener 
    TSS15 = dCdt(3,15) + dCdt(4,15) + dCdt(5,15) + dCdt(6,15) + dCdt(7,15);
    % Need to possible change internval for better control stability
    maxTSS = 46000;
    minTSS = 44000;
    setTSS = 45000;
    if TSS15 > maxTSS || TSS15 < minTSS
        ft = (dCdt(3,13) + dCdt(4,13) + dCdt(5,13) + dCdt(6,13) + dCdt(7,13))/(setTSS);
            if ft > 1
                ft = 1;
            elseif ft < 0
                ft = 0.001;
            else
            end
    else
    end
%     % Verify
%     yyaxis left
%     plot(Array.tArray,Array.Qarray(:,15));
%     yyaxis right
%     plot(Array.tArray,Array.TSS15)
    Q(15) = ft*Q(13);
    Q(14) = Q(13) - Q(15);
    Q(17) = Q(15);
    Q(16) = 0;
   
    % Removal efficiency of 0.8 of TSS
    % dCdt(i,15) going to digester
    TSSx = 1 - 0.8;
    if i < 3
        dCdt(i,14) = dCdt(i,13);
    elseif (2 < i) && (i < 8)
        dCdt(i,14) = TSSx*dCdt(i,13)*Q(13)/Q(14);
    elseif (7 < i) && (i < 12)
        dCdt(i,14) = dCdt(i,1); 
    elseif (11 < i) && (i < 13)
        dCdt(i,14) = TSSx*dCdt(i,13)*Q(13)/Q(14);
    else
        dCdt(i,14) = dCdt(i,13);
    end
    if i < 3
        dCdt(i,15) = dCdt(i,13);
    elseif (2 < i) && (i < 8)
        dCdt(i,15) = (1 - TSSx)*dCdt(i,13)*Q(13)/Q(15);
    elseif (7 < i) && (i < 12)
        dCdt(i,15) = dCdt(i,13);
    elseif (11 < i) && (i < 13)
        dCdt(i,15) = (1 - TSSx)*dCdt(i,13)*Q(13)/Q(15);
    else
        dCdt(i,15) = dCdt(i,13);
    end
    dCdt(i,13) = (dCdt(i,14)*Q(14) + dCdt(i,15)*Q(15))/Q(13);

    %% Reducing total incoming COD for Ss,Xs,Xbh,Xba in that specific order
    % Maybe optimize for loop
    for lenComp = 1:13
        xtemp(lenComp) = dCdt(lenComp,15);
    end
    for lenADM = 1:24
        adm(lenADM) = 0;
    end
    
%% Parameters
 tfac = 1/298.15 - 1/(273.15 + 35);
 bigr = 0.08314;
% CODequiv is the conversion factor for COD demand of nitrate
% exact value of ASM1 2.86
 CODequiv = 40/14;
% fraction of N in amino acids and Xpr as in ADM1 report
 fnaa = 0.098;
% N content of biomass based on BSM1, same in AS and AD
 fnbac = 0.08;
% N content of composite material based on BSM2
 fnxc = 0.0376;
% N content of inerts XI and XP, same in AS and AD
 fxni = 0.06;
% N content of SI, zero in ASM1 and BSM1
 fsni = 0;
% N content of SI in the AD system
 fsni_adm = 0.06;
% fnbac, fxni and fsni are adjusted to fit the benchmark values of iXB=0.08
% and
% iXP=0.06 and iSI=0.
% i.e 8% N content mgCOD/l <-> mgN/l = iXB, in ADM1 8.75%
 nbac = fnbac/14*14000;
% i.e. 3.76% N content mgCOD/l <-> mgN/l
 nxc = fnxc/14*14000;
% i.e. 9.8% N content mgCOD/l <-> mgN/l
 naa = fnaa/14*14000;
% i.e. 6% N content mgCOD/l <-> mgN/l = iXP = iXI
 xni = fxni/14*14000;
% i.e. 0% N content mgCOD/l <-> mgN/l = iSI
 sni = fsni/14*14000;
 sni_adm = fsni_adm/14*14000;
% lipid fraction of non-nitrogenous XS in BSM2
 frlixs = 0.7;
% anaerobically degradable fraction of biomass in BSM2
 frxs = 0.68;
% lipid fraction of non-nitrogenous biomass in BSM2
 frlixb = 0.4;
% amount of XI and XP degradable in AD, zero in BSM2
 fdegrade = 0;
% Let CODdemand be the COD demand of available electron
% acceptors prior to the anaerobic digester, i.e. oxygen and nitrate
 CODdemand = dCdt(8,15) + CODequiv*dCdt(9,15);
     if CODdemand > (dCdt(2,15) + dCdt(4,15) + dCdt(5,15) + dCdt(6,15))
        disp('Warning! Influent characterization may need to be evaluated, not enough COD available')
     else
     end
 xtemp(8) = 0;
 xtemp(9) = 0;
    if CODdemand > dCdt(2,15)
        remaina = CODdemand - dCdt(2,15);
        xtemp(2) = 0;
        if remaina > dCdt(4,15)
            remainb = remaina - dCdt(4,15);
            xtemp(4) = 0;
                if remainb > dCdt(5,15)
                    remainc = remainb - dCdt(5,15);
                    xtemp(10) = xtemp(10) + dCdt(5,15)*fnbac;
                    xtemp(5) = 0;
                        if remainc > dCdt(6,15)
                            remaind = remainc - dCdt(6,15);
                            xtemp(10) = xtemp(10) + dCdt(6,15)*fnbac;
                            xtemp(6) = 0;
                            xtemp(8) = remaind;
                            disp('ERROR: COD shortage when removing inital oxygen and nitrate')
                        else
                            xtemp(6) = dCdt(6,15) - remainc;
                            xtemp(10) = xtemp(10) + remainc*fnbac;
                        end
                else
                    xtemp(5) = dCdt(5,15) - remainb;
                    xtemp(10) = xtemp(10) + remainb*fnbac;
                end
        else
            xtemp(4) = dCdt(4,15) - remaina;
        end
    else
        xtemp(2) = dCdt(2,15) - CODdemand;
    end
    
    sorgn = dCdt(11,15)/fnaa;
    if sorgn >= xtemp(2)
        adm(2) = xtemp(2);
        xtemp(11) = xtemp(11) - xtemp(2)*fnaa;
        xtemp(2) = 0;
    else
        adm(2) = sorgn;
        xtemp(2) = xtemp(2) - sorgn;
        xtemp(11) = 0;
    end
    
    xorgn = dCdt(12,15)/fnaa;
    if xorgn >= xtemp(4)
        xprtemp = xtemp(4);
        xtemp(12) = xtemp(12) - xtemp(4)*fnaa;
        xtemp(4) = 0;
        xlitemp = 0;
        xchtemp = 0;
    else
        xprtemp = xorgn;
        xlitemp = frlixs*(xtemp(4) - xorgn);
        xchtemp = (1 - frlixs)*(xtemp(4) - xorgn);
        xtemp(4) = 0;
        xtemp(12) = 0;
    end
    
    biomass = xtemp(5) + xtemp(6);
    biomass_nobio = biomass*(1 - frxs);
    biomass_bioN = biomass*fnbac - biomass_nobio*fxni;
    if biomass_bioN < 0
        disp('ERROR: not enough biomass N to map the requested inert part')
    end
    if ((biomass_bioN/fnaa) <= (biomass - biomass_nobio))
        xprtemp2 = biomass_bioN/fnaa;
        remainCOD = biomass - biomass_nobio - xprtemp2;
        if ((xtemp(12)/fnaa) > remainCOD)
            xprtemp2 = xprtemp2 + remainCOD;
            xtemp(12) = xtemp(12) - remainCOD*fnaa;
            remainCOD = 0;
            xtemp(5) = 0;
            xtemp(6) = 0;
        else
            xprtemp2 = xprtemp2 + xtemp(12)/fnaa;
            remainCOD = remainCOD - xtemp(12)/fnaa;
            xtemp(12) = 0;
        end
        xlitemp2 = frlixb*remainCOD;
        xchtemp2 = (1 - frlixb)*remainCOD;
    else
        xprtemp2 = biomass - biomass_nobio;
        xtemp(12) = xtemp(12) + biomass*fnbac - biomass_nobio*fxni - xprtemp2*fnaa;
    end
    xtemp(5) = 0;
    xtemp(6) = 0;
    
    inertX = (1 - fdegrade)*(xtemp(3) + xtemp(7));
    xc = 0;
    xlitemp3 = 0;
    xchtemp3 = 0;
    if fdegrade > 0
        noninertX = fdegrade*(xtemp(3) + xtemp(7));
        if ((noninertX*fxni) < (noninertX*fnxc))
            xc = noninertX*fxni/fnxc;
            noninertX = noninertX - noninertX*fxni/fnxc;
            if xtemp(12) < noninertX*fnxc
                xc = xc + xtemp(12)/fnxc;
                noninertX = noninertX - xtemp(12)/fnxc;
                xtemp(12) = 0;
                if (xtemp(11) < noninertX*fnxc)
                    xc = xc + temp(11)/fnxc;
                    noninertX = nonintertX - xtemp(11)/fnxc;
                    xtemp(11) = 0;
                    if (xtemp(10) < noninertX*fnxc)
                        xc = xc + xtemp(10)/fnxc;
                        noninertX = noninertX - xtemp(10)/fnxc;
                        xtemp(10) = 0;
                        disp('ERROR: Nitrogen shortage when converting biodegradable XI&XP')
                        disp('Putting remaining XI&XP as lipids (50c) and carbohydrates (50%)') 
                        xlitemp3 = 0.5*noninertX;
                        xchtemp3 = 0.5*noninertX;
                        noninertX = 0;
                    else
                        xc = xc + noninertX;
                        xtemp(10) = xtemp(10) - noninertX*fnxc;
                        noninertX = 0;
                    end
                else
                    xc = xc + noninertX;
                    xtemp(11) = xtemp(11) - noninertX*fnxc;
                    noninertX = 0;
                end
            else
                xc = xc + noninertX;
                xtemp(12) = xtemp(12) - noninertX*fnxc;
                noninertX = 0;
            end
        else
            xc = xc + noninertX;
            xtemp(12) = xtemp(12) + noninertX*(fxni - fnxc);
            noninertX = 0;
        end
    end
    
    inertS = 0;
    if ((xtemp(1)*fsni) < (xtemp(1)*fsni_adm))
        inertS = xtemp(1)*fsni/fsni_adm;
        xtemp(1) = xtemp(1) - xtemp(1)*fsni/fsni_adm;
        if (xtemp(11) < (xtemp(1)*fsni_adm))
            inertS = inertS + xtemp(11)/fsni_adm;
            xtemp(1) = xtemp(1) - xtemp(11)/fsni_adm;
            xtemp(11) = 0;
            if (xtemp(12) < (xtemp(1)*fsni_adm))
                inertS = inertS + xtemp(12)/fsni_adm;
                xtemp(1) = xtemp(1) - xtemp(12)/fsni_adm;
                xtemp(12) = 0;
                if (xtemp(10) < (xtemp(1)*fsni_adm))
                    inertS = inertS + xtemp(10)/fsni_adm;
                    xtemp(1) = xtemp(1) - xtemp(10)/fsni_adm;
                    xtemp(10) = 0;
                    disp('ERROR: Nitrogen shortage when converting SI')
                    disp('Putting remaining SI as monosacharides') 
                    xtemp(2) = xtemp(2) + xtemp(1);
                    xtemp(1) = 0;
                else
                    inertS = inertS + xtemp(1);
                    xtemp(10) = xtemp(10) - xtemp(1)*fsni_adm;
                    xtemp(1) = 0;
                end
            else
                inertS = inertS + xtemp(1);
                xtemp(12) = xtemp(12) - xtemp(1)*fsni_adm;
                xtemp(1) = 0;
            end
        else
            inertS = inertS + xtemp(1);
            xtemp(11) = xtemp(11) - xtemp(1)*fsni_adm;
            xtemp(1) = 0;
        end
    else
        inertS = inertS + xtemp(1);
        xtemp(11) = xtemp(11) + xtemp(1)*(fsni - fsni_adm);
        xtemp(1) = 0;
    end
    
    adm(1) = xtemp(2)/1000;
    adm(2) = adm(2)/1000;
    adm(10) = xtemp(13)/1000;
    adm(11) = (xtemp(10) + xtemp(11) + xtemp(12))/14000;
    adm(12) = inertS/1000;
    adm(13) = xc/1000;
    adm(14) = (xchtemp + xchtemp2 + xchtemp3)/1000;
    adm(15) = (xprtemp + xprtemp2)/1000;
    adm(16) = (xlitemp + xlitemp2 + xlitemp3)/1000;
    adm(24) = (biomass_nobio + inertX)/1000;
% Calculation of adm(10) (Sic)
    % Take average pH from before digester
    ph = 7;
    alfachac = (-1/64)/(1 + 10^(4.76 - ph));
    alfachpro = (-1/112)/(1 + 10^(4.88 - ph));
    alfachbu = (-1/160)/(1 + 10^(4.82 - ph));
    alfachva = (-1/208)/(1 + 10^(4.86 - ph));
    pkk = 9.25 - log10(exp(51965/bigr/100*tfac));
    alfachin = (10^(pkk - ph))/(1 + 10^(pkk - ph));
    pkk = 6.35 - log10(exp(7646/bigr/100*tfac));
    alfachic = -1/(1 + 10^(pkk - ph));
    % modifie / PV
    alfachnh = 1/14;
    alfachno = -1/14;
    alfachalk = -1;
    chargeasm1 = (dCdt(13,15)*alfachalk + dCdt(10,15)*alfachnh +...
        dCdt(9,15)*alfachno)/1000;
    chargeadm1 = adm(7)*alfachac + adm(6)*alfachpro +...
        adm(5)*alfachbu + adm(4)*alfachva + adm(11)*alfachin;
    adm(10) = (chargeasm1 - chargeadm1)/alfachic;
    pkk = 14 - log10(exp(55900/bigr/100*tfac));
    ancat = chargeadm1 + adm(10)*alfachic - 10^(-ph) + 10^(-pkk + ph);
    
    % Check mass balances
    totCODin = 0;
    for vec = 1:7
        totCODin = totCODin + dCdt(vec,15);
    end
    
    totNin = dCdt(9,15) + dCdt(10,15) + dCdt(11,15) + dCdt(12,15)...
        + (nbac/1000)*(dCdt(5,15) + dCdt(6,15)) +  (sni/1000)*dCdt(1,15)...
        + (xni/1000)*(dCdt(3,15) + dCdt(7,15));
    
    totCODout = 0;
    for vec1 = 1:9
        totCODout = totCODout + adm(vec1)*1000;
    end
    for vec2 = 12:24
        totCODout = totCODout + adm(vec2)*1000;
    end
    
    totNout = nbac*(adm(17) + adm(18) + adm(19) + adm(20)...
        + adm(21) + adm(22) + adm(23)) + naa*(adm(2) + adm(15))...
        + adm(11)*14000 + sni_adm*adm(12) + nxc*adm(13) + xni*adm(24);
    MBCODin = totCODin - totCODout;
    MBNin = totNin - totNout;
    
    % Set adm vector to influent components
    % Units of [kgCOD/m3]
    S_su_in = adm(1);
    S_aa_in = adm(2);
    S_fa_in = adm(3);
    S_va_in = adm(4);
    S_bu_in = adm(5);
    S_pro_in = adm(6);
    S_ac_in = adm(7);
    S_h2_in = adm(8);
    S_ch4_in = adm(9);
    S_IC_in = adm(10);
    S_IN_in = adm(11);
    S_I_in = adm(12);
    X_c_in = adm(13);
    X_ch_in = adm(14);
    X_pr_in = adm(15);
    X_li_in = adm(16);
    X_su_in = adm(17);
    X_aa_in = adm(18);
    X_fa_in = adm(19);
    X_c4_in = adm(20);
    X_pro_in = adm(21);
    X_ac_in = adm(22);
    X_h2_in = adm(23);
    X_I_in = adm(24);
    % The loop below is based off of the Matlab code on page 57 of
    % Benchmark Simulation Model no. 2 (BSM2) paper
    S_an_in = ancat;
    if S_an_in < 0
        S_an_in = -1*S_an_in;
        S_cat_in = 0;
    else
        S_cat_in = S_an_in;
        S_an_in = 0;
    end
    % Save to concentration
    dCdt(14,15) = S_su_in;
    dCdt(15,15) = S_aa_in;
    dCdt(16,15) = S_fa_in;
    dCdt(17,15) = S_va_in;
    dCdt(18,15) = S_bu_in;
    dCdt(19,15) = S_pro_in;
    dCdt(20,15) = S_ac_in;
    dCdt(21,15) = S_h2_in;
    dCdt(22,15) = S_ch4_in;
    dCdt(23,15) = S_IC_in;
    dCdt(24,15) = S_IN_in;
    dCdt(25,15) = S_I_in;
    dCdt(26,15) = X_c_in;
    dCdt(27,15) = X_ch_in;
    dCdt(28,15) = X_pr_in;
    dCdt(29,15) = X_li_in;
    dCdt(30,15) = X_su_in;
    dCdt(31,15) = X_aa_in;
    dCdt(32,15) = X_fa_in;
    dCdt(33,15) = X_c4_in;
    dCdt(34,15) = X_pro_in;
    dCdt(35,15) = X_ac_in;
    dCdt(36,15) = X_h2_in;
    dCdt(37,15) = X_I_in;
    dCdt(38,15) = S_cat_in;
    dCdt(39,15) = S_an_in;
    % Variable change to check against AD paper
    %l.V_gas = 19.98697; % m3 headspace volume taken from B&V study
    %l.V_liq = Vol6 - l.V_gas; % m3
    
%% AD differential equations
% Data/Equations pulled from Aspects on ADM1 implementation within the BSM2 framework
% CHECK PARAMETERS -> CHANGE VOLUMES (headspace/liquid)
    % Intial conditions inside reactor
    S_su = dCdt(14,17);
    S_aa = dCdt(15,17);
    S_fa = dCdt(16,17);
    S_va = dCdt(17,17);
    S_bu = dCdt(18,17);
    S_pro = dCdt(19,17);
    S_ac = dCdt(20,17);
    S_h2 = dCdt(21,17);
    S_ch4 = dCdt(22,17);
    S_IC = dCdt(23,17);
    S_IN = dCdt(24,17);
    S_I = dCdt(25,17);
    X_c = dCdt(26,17);
    X_ch = dCdt(27,17);
    X_pr = dCdt(28,17);
    X_li = dCdt(29,17);
    X_su = dCdt(30,17);
    X_aa = dCdt(31,17);
    X_fa = dCdt(32,17);
    X_c4 = dCdt(33,17);
    X_pro = dCdt(34,17);
    X_ac = dCdt(35,17);
    X_h2 = dCdt(36,17);
    X_I = dCdt(37,17);
    S_cat = dCdt(38,17);
    S_an = dCdt(39,17);
    S_vam = dCdt(40,17);
    S_bum = dCdt(41,17);
    S_prom = dCdt(42,17);
    S_acm = dCdt(43,17);
    S_hco3m = dCdt(44,17);
    S_nh3 = dCdt(45,17);
    S_gas_h2 = dCdt(46,16);
    S_gas_ch4 = dCdt(47,16);
    S_gas_co2 = dCdt(48,16);
    
    % Gas pressure
    P_gas_h2 = S_gas_h2*l.R*l.T_op/16;
    P_gas_ch4 = S_gas_ch4*l.R*l.T_op/64;
    P_gas_co2 = S_gas_co2*l.R*l.T_op;
    P_gas = P_gas_h2 + P_gas_ch4 + P_gas_co2 + l.p_gas_h2o;
    q_gas = l.k_p*(P_gas - l.P_atm)*P_gas/l.P_atm;
    
    % Calculation of potential Sh
    S_nh4 = S_IN - S_nh3;
    S_co2 = S_IC - S_hco3m;
    Theta = S_cat + S_nh4 - S_hco3m - (S_acm/64) - (S_prom/112) - (S_bum/160)...
    - (S_vam/208) - S_an;
    Sh = -Theta/2 + sqrt(Theta^2 + 4*l.K_w)/2;
    if(Sh <= 0)
    Sh = 1e-12;
    end
    pH = -log10(Sh);
    
    % pH inhibition
    if(pH < l.pH_UL_aa)
    I_pH_aa = exp(-3*((pH - l.pH_UL_aa)/(l.pH_UL_aa - l.pH_LL_aa))^2);
    else
    I_pH_aa = 1;
    end
    if(pH < l.pH_UL_ac)
    I_pH_ac = exp(-3*((pH - l.pH_UL_ac)/(l.pH_UL_ac - l.pH_LL_ac))^2);
    else
    I_pH_ac = 1;
    end
    if(pH < l.pH_UL_h2)
    I_pH_h2 = exp(-3*((pH - l.pH_UL_h2)/(l.pH_UL_h2 - l.pH_LL_h2))^2);
    else
    I_pH_h2 = 1;
    end
    
    % Process inhibition
    I_IN_lim = 1/(1 + l.K_S_IN/S_IN);
    I_h2_fa = 1/(1 + S_h2/l.K_I_h2_fa);
    I_h2_c4 = 1/(1 + S_h2/l.K_I_h2_c4);
    I_h2_pro = 1/(1 + S_h2/l.K_I_h2_pro);
    I_nh3 = 1/(1 + S_nh3/l.K_I_nh3);
    I_5 = I_pH_aa*I_IN_lim;
    I_6 = I_5;
    I_7 = I_pH_aa*I_IN_lim*I_h2_fa;
    I_8 = I_pH_aa*I_IN_lim*I_h2_c4;
    I_9 = I_8;
    I_10 = I_pH_aa*I_IN_lim*I_h2_pro;
    I_11 = I_pH_ac*I_IN_lim*I_nh3;
    I_12 = I_pH_h2*I_IN_lim;
    
    % Biochemical process rates
    rho_1 = l.k_dis*X_c;
    rho_2 = l.k_hyd_ch*X_ch;
    rho_3 = l.k_hyd_pr*X_pr;
    rho_4 = l.k_hyd_li*X_li;
    rho_5 = l.k_m_su*(S_su/(l.K_S_su + S_su))*X_su*I_5;
    rho_6 = l.k_m_aa*(S_aa/(l.K_S_aa + S_aa))*X_aa*I_6;
    rho_7 = l.k_m_fa*(S_fa/(l.K_S_fa + S_fa))*X_fa*I_7;
    rho_8 = l.k_m_c4*(S_va/(l.K_S_c4 + S_va))*X_c4*(S_va/(S_bu + S_va + 1e-6))*I_8;
    rho_9 = l.k_m_c4*(S_bu/(l.K_S_c4 + S_bu))*X_c4*(S_bu/(S_va + S_bu + 1e-6))*I_9;
    rho_10 = l.k_m_pro*(S_pro/(l.K_S_pro + S_pro))*X_pro*I_10;
    rho_11 = l.k_m_ac*(S_ac/(l.K_S_ac + S_ac))*X_ac*I_11;
    rho_12 = l.k_m_h2*(S_h2/(l.K_S_h2 + S_h2))*X_h2*I_12;
    rho_13 = l.k_dec_X_su*X_su;
    rho_14 = l.k_dec_X_aa*X_aa;
    rho_15 = l.k_dec_X_fa*X_fa;
    rho_16 = l.k_dec_X_c4*X_c4;
    rho_17 = l.k_dec_X_pro*X_pro;
    rho_18 = l.k_dec_X_ac*X_ac;
    rho_19 = l.k_dec_X_h2*X_h2;
    
    % Acid/base rates
    rho_A_4 = l.k_A_B_va*(S_vam*(l.K_a_va + Sh) - l.K_a_va*S_va);
    rho_A_5 = l.k_A_B_bu*(S_bum*(l.K_a_bu + Sh) - l.K_a_bu*S_bu);
    rho_A_6 = l.k_A_B_pro*(S_prom*(l.K_a_pro + Sh) - l.K_a_pro*S_pro);
    rho_A_7 = l.k_A_B_ac*(S_acm*(l.K_a_ac + Sh) - l.K_a_ac*S_ac);
    rho_A_10 = l.k_A_B_co2*(S_hco3m*(l.K_a_co2 + Sh) - l.K_a_co2*S_IC);
    rho_A_11 = l.k_A_B_IN*(S_nh3*(l.K_a_IN + Sh) - l.K_a_IN*S_IN);
    
    % Gas transfer rates
    rho_T_8 = l.k_L_a*(S_h2 - 16*l.K_H_h2*P_gas_h2);
    rho_T_9 = l.k_L_a*(S_ch4 - 64*l.K_H_ch4*P_gas_ch4);
    rho_T_10 = l.k_L_a*(S_co2 - l.K_H_co2*P_gas_co2);
    
    % Processes (carbon)
    s_1 = -l.C_xc + l.f_sI_xc*l.C_sI + l.f_ch_xc*l.C_ch + l.f_pr_xc*l.C_pr + l.f_li_xc*l.C_li + l.f_xI_xc*l.C_xI;
    s_2 = -l.C_ch + l.C_su;
    s_3 = -l.C_pr + l.C_aa;
    s_4 = -l.C_li + (1 - l.f_fa_li)*l.C_su + l.f_fa_li*l.C_fa;
    s_5 = -l.C_su + (1 - l.Y_su)*(l.f_bu_su*l.C_bu + l.f_pro_su*l.C_pro + l.f_ac_su*l.C_ac) + l.Y_su*l.C_bac;
    s_6 = -l.C_aa + (1 - l.Y_aa)*(l.f_va_aa*l.C_va + l.f_bu_aa*l.C_bu + l.f_pro_aa*l.C_pro + l.f_ac_aa*l.C_ac) + ...
    l.Y_aa*l.C_bac;
    s_7 = -l.C_fa + (1 - l.Y_fa)*0.7*l.C_ac + l.Y_fa*l.C_bac;
    s_8 = -l.C_va + (1 - l.Y_c4)*0.54*l.C_pro + (1 - l.Y_c4)*0.31*l.C_ac + l.Y_c4*l.C_bac;
    s_9 = -l.C_bu + (1 - l.Y_c4)*0.8*l.C_ac + l.Y_c4*l.C_bac;
    s_10 = -l.C_pro + (1 - l.Y_pro)*0.57*l.C_ac + l.Y_pro*l.C_bac;
    s_11 = -l.C_ac + (1 - l.Y_ac)*l.C_ch4 + l.Y_ac*l.C_bac;
    s_12 = (1 - l.Y_h2)*l.C_ch4 + l.Y_h2*l.C_bac;
    s_13 = -l.C_bac + l.C_xc;
    
    % Differential equations 1-4, soluble matter
                                                                                                % State No.
    Conc(14,17) = Q(15)/l.V_liq*(S_su_in - S_su) + rho_2 + (1 - l.f_fa_li)*rho_4 - rho_5;       % 1
    Conc(15,17) = Q(15)/l.V_liq*(S_aa_in - S_aa) + rho_3 - rho_6;                               % 2
    Conc(16,17) = Q(15)/l.V_liq*(S_fa_in - S_fa) + l.f_fa_li*rho_4 - rho_7;                     % 3
    Conc(17,17) = Q(15)/l.V_liq*(S_va_in - S_va) + (1 - l.Y_aa)*l.f_va_aa*rho_6 - rho_8;        % 4
    
    % Differential equations 5-8, soluble matter
    Conc(18,17) = Q(15)/l.V_liq*(S_bu_in - S_bu) + (1 - l.Y_su)*l.f_bu_su*rho_5 + ...           % 5
    (1 - l.Y_aa)*l.f_bu_aa*rho_6 - rho_9; 
    Conc(19,17) = Q(15)/l.V_liq*(S_pro_in - S_pro) + (1 - l.Y_su)*l.f_pro_su*rho_5 + ...          % 6 
    (1 - l.Y_aa)*l.f_pro_aa*rho_6 + (1-l.Y_c4)*0.54*rho_8 - rho_10;
    Conc(20,17) = Q(15)/l.V_liq*(S_ac_in - S_ac) + (1 - l.Y_su)*l.f_ac_su*rho_5 + ...             % 7
    (1 - l.Y_aa)*l.f_ac_aa*rho_6 + (1 - l.Y_fa)*0.7*rho_7 + ... 
    (1 - l.Y_c4)*0.31*rho_8 + (1 - l.Y_c4)*0.8*rho_9 + ... 
    (1 - l.Y_pro)*0.57*rho_10 - rho_11;
    Conc(21,17) = Q(15)/l.V_liq*(S_h2_in - S_h2) + (1 - l.Y_su)*l.f_h2_su*rho_5 + ...           % 8
    (1 - l.Y_aa)*l.f_h2_aa*rho_6 + (1 - l.Y_fa)*0.3*rho_7 + ... 
    (1 - l.Y_c4)*0.15*rho_8 + (1 - l.Y_c4)*0.2*rho_9 + ... 
    (1 - l.Y_pro)*0.43*rho_10 - rho_12 - rho_T_8;

    % Differential equations 9-12, soluble matter
    Conc(22,17) = Q(15)/l.V_liq*(S_ch4_in - S_ch4) + (1 - l.Y_ac)*rho_11 + ...                  % 9
    (1 - l.Y_h2)*rho_12 - rho_T_9; 
    Conc(23,17) = Q(15)/l.V_liq*(S_IC_in - S_IC) - (s_1*rho_1 + s_2*rho_2 + ...                 % 10
    s_3*rho_3 + s_4*rho_4 + s_5*rho_5 + s_6*rho_6 + s_7*rho_7 + ... 
    s_8*rho_8 + s_9*rho_9 + s_10*rho_10 + s_11*rho_11 + s_12*rho_12 + ... 
    s_13*(rho_13 + rho_14 + rho_15 + rho_16... 
    + rho_17 + rho_18 + rho_19)) - rho_T_10; 
    Conc(24,17) = Q(15)/l.V_liq*(S_IN_in - S_IN) - l.Y_su*l.N_bac*rho_5 + ...                   % 11
    (l.N_aa - l.Y_aa*l.N_bac)*rho_6 - l.Y_fa*l.N_bac*rho_7 - ... 
    l.Y_c4*l.N_bac*rho_8 - l.Y_c4*l.N_bac*rho_9 - ... 
    l.Y_pro*l.N_bac*rho_10 - l.Y_ac*l.N_bac*rho_11 - ... 
    l.Y_h2*l.N_bac*rho_12 + (l.N_bac - l.N_xc)*... 
    (rho_13 + rho_14 + rho_15 + rho_16 + rho_17 + rho_18 +... 
    rho_19) + (l.N_xc - l.f_xI_xc*l.N_I - l.f_sI_xc*l.N_I - l.f_pr_xc*l.N_aa)*rho_1; 
    Conc(25,17) = Q(15)/l.V_liq*(S_I_in - S_I) + l.f_sI_xc*rho_1;                               % 12
    
    % Differential equations 13-16, particulate matter
    Conc(26,17) = Q(15)/l.V_liq*(X_c_in - X_c) - rho_1 + ...                                   % 13
    rho_13 + rho_14 + rho_15 + rho_16 + rho_17 + rho_18 + rho_19; 
    Conc(27,17) = Q(15)/l.V_liq*(X_ch_in - X_ch) + l.f_ch_xc*rho_1 - rho_2;                     % 14
    Conc(28,17) = Q(15)/l.V_liq*(X_pr_in - X_pr) + l.f_pr_xc*rho_1 - rho_3;                     % 15
    Conc(29,17) = Q(15)/l.V_liq*(X_li_in - X_li) + l.f_li_xc*rho_1 - rho_4;                     % 16
    
    % Differential equations 17-20, particulate matter
    Conc(30,17) = Q(15)/l.V_liq*(X_su_in - X_su) + l.Y_su*rho_5 - rho_13;                       % 17
    Conc(31,17) = Q(15)/l.V_liq*(X_aa_in - X_aa) + l.Y_aa*rho_6 - rho_14;                       % 18
    Conc(32,17) = Q(15)/l.V_liq*(X_fa_in - X_fa) + l.Y_fa*rho_7 - rho_15;                       % 19
    Conc(33,17)  = Q(15)/l.V_liq*(X_c4_in - X_c4) + l.Y_c4*rho_8 + l.Y_c4*rho_9 - rho_16;       % 20
    
    % Differential equations 21-24, particulate matter
    Conc(34,17) = Q(15)/l.V_liq*(X_pro_in - X_pro) + l.Y_pro*rho_10 - rho_17;                   % 21
    Conc(35,17) = Q(15)/l.V_liq*(X_ac_in - X_ac) + l.Y_ac*rho_11 - rho_18;                      % 22
    Conc(36,17) = Q(15)/l.V_liq*(X_h2_in - X_h2) + l.Y_h2*rho_12 - rho_19;                      % 23
    Conc(37,17) = Q(15)/l.V_liq*(X_I_in - X_I) + l.f_xI_xc*rho_1;                               % 24
    
    % Differential equations 25-26, cations and anions
    Conc(38,17) = Q(15)/l.V_liq*(S_cat_in - S_cat);                                             % 25
    Conc(39,17) = Q(15)/l.V_liq*(S_an_in - S_an);                                               % 26
    
    % Differential equations 27-32, ion states
    Conc(40,17) = -rho_A_4;                                                                     % 27
    Conc(41,17) = -rho_A_5;                                                                     % 28
    Conc(42,17) = -rho_A_6;                                                                     % 29
    Conc(43,17) = -rho_A_7;                                                                     % 30
    Conc(44,17) = -rho_A_10;                                                                    % 31
    Conc(45,17) = -rho_A_11;                                                                    % 32
    
    % Differential equations 33-35, gas phase equations
    Conc(46,16) = -S_gas_h2*q_gas/l.V_gas + rho_T_8*l.V_liq/l.V_gas;                            % 33
    Conc(47,16) = -S_gas_ch4*q_gas/l.V_gas + rho_T_9*l.V_liq/l.V_gas;                           % 34
    Conc(48,16) = -S_gas_co2*q_gas/l.V_gas + rho_T_10*l.V_liq/l.V_gas;                          % 35

%% Conversion from ADM1 to ASM1
    for vec3 = 1:13
        asmm(vec3) = 0;
    end
    for vec4 = 1:24
        xtemp(vec4) = dCdt(13 + vec4,17);
    end
% Set parameter values
    fnaa = 0.098;
    fnxc = 0.0376;
    fnbac = 0.08;
    fsni_adm = 0.06;
    fxni = 0.06;
    fsni = 0;
    nbac = fnbac/14*14000;
    nxc = fnxc/14*14000;
    naa = fnaa/14*14000;
    sni_adm = fsni_adm/14*14000;
    xni = fxni/14*14000;
    sni = fsni/14*14000;
    frxs_AS = 0.79;
    fdegrade_AS = 0;

    biomass = 0;
    for vec5 = 17:23
        biomass = biomass + xtemp(vec5)*1000;
    end
    biomass_nobio = biomass*(1 - frxs_AS);
    biomass_bioN = biomass*fnbac - biomass_nobio*fxni;
    remainCOD = 0;
    if biomass_bioN < 0
        disp('Warning 1')
        xptemp = biomass*fnbac/fxni;
        % added from matlab code of BSM2
        biomass_nobio = xptemp;
        biomass_bioN = 0;
    else
        xptemp = biomass_nobio;
    end
    if ((biomass_bioN/fnxc) <= (biomass - biomass_nobio))
        xstemp = biomass_bioN/fnxc;
        remainCOD = biomass - biomass_nobio - xstemp;
        % added >= instead of > from matlab code of BSM2 
        if ((xtemp(11)*14000/fnxc) >= remainCOD)
            xstemp = xstemp + remainCOD;
        else
            disp('Warning 2')
            disp('System failure')
        end
    else
        xstemp = biomass - biomass_nobio;
    end
    xtemp(11) = xtemp(11) + biomass*fnbac/14000 - xptemp*fxni/14000 ...
        - xstemp*fnxc/14000;
        
    asmm(4) = xstemp;
    for vec6 = 13:16
        asmm(4) = asmm(4) + xtemp(vec6)*1000;
    end
    
    asmm(7) = xptemp;
    
    inertX = (1 - fdegrade_AS)*xtemp(24)*1000;
    xstemp2 = 0;
    noninertX = 0;
    if (fdegrade_AS > 0)
        noninertX = fdegrade_AS*xtemp(24)*1000;
        if fxni < fnxc
            xstemp2 = noninertX*fxni/fnxc;
            noninertX = noninertX - noninertX*fnxi/fnxc;
            if ((xtemp(11)*14000) < (noninertX*fnxc))
                xstemp2 = xstemp2 + xtemp(11)*14000/fnxc;
                noninertX = noninertX - xtemp(11)*14000/fnxc;
                xtemp(11) = 0;
                disp('Warning 3')
                inertX = inertX + noninertX;
            else
                xstemp2 = xstemp2 + noninertX;
                xtemp(11) = xtemp(11) - noninertX*fnxc/14000;
                noninertX = 0;
            end
        else
            xstemp2 = xstemp2 + noninertX;
            xtemp(11) = xtemp(11) + noninertX*(fnxi - fnxc)/14000;
            noninertX = 0;
        end
    end
    
    asmm(3) = inertX;
    asmm(4) = asmm(4) + xstemp2;
    
    inertS = 0;
    if (fsni_adm < fsni)
        inertS = xtemp(12)*fsni_adm/fsni;
        % added from matlab code of BSM2
        xtemp(12) = xtemp(12) - xtemp(12)*fsni_adm/fsni;
        if ((xtemp(11)*14) < (xtemp(12)*fsni))
            inertS = inertS + xtemp(11)*14/fsni;
            xtemp(12) = xtemp(12) - xtemp(11)*14/fsni;
            xtemp(11) = 0;
            disp('Warning 5')
        else
            inertS = inertS + xtemp(12);
            xtemp(11) = xtemp(11) - xtemp(12)*fsni/14;
            xtemp(12) = 0;
        end
    else
        inertS = inertS + xtemp(12);
        xtemp(11) = xtemp(11) + xtemp(12)*(fsni_adm - fsni)/14;
        xtemp(12) = 0;
    end
    
    asmm(1) = inertS*1000;
    
    asmm(12) = fnxc*(xstemp + xstemp2) + nxc*xtemp(13) + naa*xtemp(15);
    
    for vec7 = 1:7
        asmm(2) = asmm(2) + xtemp(vec7)*1000;
    end
    
    asmm(11) = naa*xtemp(2);
    
    asmm(10) = xtemp(11)*14000;
    
    %% VERIFY pH IS CORRECT DURING ENTIRE SIMULATION
% Calculation of Salk
    alfachac = (-1/64)/(1 + 10^(4.76 - pH));
    alfachpro = (-1/112)/(1 + 10^(4.88 - pH));
    alfachbu = (-1/160)/(1 + 10^(4.82 - pH));
    alfachva = (-1/208)/(1 + 10^(4.86 - pH));
    pkk = 9.25 - log10(exp(51965/bigr/100*tfac));
    alfachin = (10^(pkk - pH))/(1 + 10^(pkk - pH));
    pkk = 6.35 - log10(exp(7646/bigr/100*tfac));
    alfachic = -1/(1 + 10^(pkk - pH));
% modifie / PV
    alfachnh = 1/14;
    alfachno = -1/14;
    alfachalk = -1;
    chargeasm1 = asmm(10)*alfachnh + asmm(9)*alfachno;
    chargeadm1 = (adm(7)*alfachac + adm(6)*alfachpro + adm(11)*alfachin...
        + adm(5)*alfachbu + adm(4)*alfachva + adm(10)*alfachic)*1000;
    asmm(13) = chargeasm1 - chargeadm1;
% Check mass balances
    totCODin = 0;
    for vec8 = 1:7
        totCODin = totCODin + adm(vec8);
    end
    for vec9 = 12:24
        totCODin = totCODin + adm(vec9)*1000;
    end
    totTKNin = 0;
    for vec10 = 17:23
        totTKNin = totTKNin + nbac*adm(vec10);
    end
    totTKNin = totTKNin + nxc*adm(13) + naa*adm(15) + naa*adm(2) +...
        adm(11)*14000 + sni_adm*adm(12) + xni*adm(24);
    totCODout = 0;
    for vec11 = 1:7
        totCODout = totCODout + asmm(vec11);
    end
    % SI_N not included here below
    totTKNout = asmm(10) + asmm(11) + asmm(12) + fsni*asmm(1) +...
        fnbac*(asmm(5) + asmm(6)) + fxni*(asmm(3) + asmm(7));
    MBCODout = totCODin - totCODout;
    MBTNKout = totTKNin - totTKNout;

    % CONVERT UNITS CORRECTLY BACK TO ASM1
    % CODconserved = CODt_anaerobic - Sh2 - Sch4;
    % COD Conversions
    dCdt(1,17) = asmm(1);
    dCdt(2,17) = asmm(2);
    dCdt(3,17) = asmm(3);
    dCdt(4,17) = asmm(4);
    dCdt(5,17) = asmm(5);
    dCdt(6,17) = asmm(6);
    dCdt(7,17) = asmm(7);
    dCdt(8,17) = asmm(8);
    dCdt(9,17) = asmm(9);
    dCdt(10,17) = asmm(10);
    dCdt(11,17) = asmm(11);
    dCdt(12,17) = asmm(12);
    dCdt(13,17) = asmm(13);
    i = i + 1;
end
vec_len = numel(Conc);
Conc = reshape(Conc,[vec_len,1]);
% Due to Conc only changing in the reactors, dCdt variable needs to be
% saved and sent to the workspace as it used throughout the code to keep track of components
% At time = 1, dont store variable into structure
if t == 1
    Array.mleArray = dCdt;
    Array.tArray = t;
    Array.Qarray = Q;
    Array.TSS15 = TSS15;
    Array.ft = ft;
    Array.pH = pH;
elseif t >= 1.01
    % Avoid large dataset and quicken solver time by only aquiring data at the specified time span
%     findT = ismembertol(t,Var.timespan,0.00001);
%     T_vector = find(findT == 1);
    AdjT = round(t*100)/100;
    FindT = find(AdjT == Var.timespan);
    if numel(FindT) > 0
        Array.mleArray = [Array.mleArray;dCdt];
        Array.tArray = [Array.tArray;t];
        Array.Qarray = [Array.Qarray;Q];
        Array.TSS15 = [Array.TSS15;TSS15];
        Array.ft = [Array.ft;ft];
        Array.pH = [Array.pH;pH];
        assignin('base','Array',Array);
    else
    end
else
end
fprintf('Current simulation time is %6.6f days\n',t)
end