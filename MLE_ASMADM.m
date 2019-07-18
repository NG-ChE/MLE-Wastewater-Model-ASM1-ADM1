tic
%% TO-DO
% GO OVER FORTRAN CODE (pg 65 out of 99) OF BENCHMARK SIMULATION MODEL NO.2 FOR
% ASM/ADM CONVERSION

% Adjust initial values for each stream throughout the entire plant
% Convert ASM values to ADM
    % Needs to be double checked for errors (units,etc.)
% Implement FOG for AD    
% Implement thickener(before AD)/dewatering(after AD)/dryer(after AD)
% Implement denitrification filter

% Create PID for oxygen transfer in aeration zone
% Create GUI for kinetics/ recycle ratios/ simulation time/ flow variables

%% Issues
% AD isn't working properly
% Carbon balance on AD is flat out wrong
% ***** CHECK UNITS *****


%% Changes from last version
clc
%clear all

%% Simple simulation of MLE system with no digester.

%% Recycle ratios for MLE system
Var.Rir = 1.88; % Specify internal recycle ratio MLR/PC_effluent 
Var.Rr = 0.71; % Specify return activated sludge recycle ratio RAS/PC_effluent 
Var.fpc = 0.017554; % Primary Clarifier flow separation (has to be less than or equal to 1)
if Var.fpc > 1
    Var.fpc = 1;
elseif Var.fpc < 0
   Var.fpc = 0;
else
end
Var.fsc = 0.974; % Secondary Clarifier underflow separation (has to be less than or equal to 1)
if Var.fsc > 1
    Var.fsc = 1;
elseif Var.fsc < 0
    Var.fsc = 0;
else
end

%% biological parameters and volumes
% Check parameter values with other reported values
% Make them a function of temperature
% Correct for units
Var.param = [0.67 0.24 0.08 0.08 ...
0.06 4 10 0.2 ...
0.5 0.8 0.8 1 ...
0.4 0.3 0.05 0.05 ...
3 0.1 0.8 40 ...
10 416.83 5522.15 36339.95 ....
4803.84 0 3028.329]';

%% AD parameters
% ADJUST VALUES
% ADJUST INITIAL VALUES INSIDE DIFF EQs
%l = SE_IndataADM1_v2;
l = IndataADM1_v2;

%% simulation time span (days)
% Increase timespan interval for more data points in result section, can
% signicantly reduce time if interval is very small
t = 1:0.1:29;
% Sample rate
    % Decrease sample rate for better DO control, but ODE takes longer
sp = 1/5;
Var.tsample = t + sp*t;
Var.timespan = t;


%% Convert typical units coming into plant to ASM1 variables
% Units g/m3 = mg/L
% TSS = 204.85; % mg/L
% % If no VSS data, assume VSS/TSS ratio of 0.69 or 0.75
% tv = 0.69;
% VSS = TSS*tv;
% FSS = TSS - VSS; % fixed suspended solids
% % USING CBOD5 FROM DATA
% BOD5 = 150.76; % mg/L
% TKN = 44.19; % mg-N/L
% NH3 = 29.1; % mg-N/L
% NO3 = 0.02; % mg-N/L
% % Alkalinity not given
% ALK = 200; % mg-N/L
% % Conversion to ASM1 Variables
% CODt = 2.1*BOD5; % mgCOD/L, Converts BOD5 to total COD, if not available
% CODbo = 1.71*BOD5; % mgCOD/L, biodegradable COD
% CODio = CODt - CODbo; % mgCOD/L, inert COD
% % Xio = 0.375*1.5*VSS; % mgCOD/L, particulate inert COD
% % Sio = CODio - Xio; % mgCOD/L, soluble inert COD
% % f_readily = 0.43; % fraction of biodegradable COD that is readily biodegradable
% % Sso = CODbo*f_readily; % mgCOD/L, readily biodegradable substrate
% % Xso = CODbo - Sso; % mgCOD/L, slowly biodegradable substrate
% % ONtotal = TKN - NH3; % mg-N/L, total organic nitrogen
% % Snio = 1.5; % mg-N/L, soluble inert organic nitrogen
% % in_xd = 0.06; % mass of nitrogen per mass of COD in biomass
% % Xnio = in_xd*Xio; % mg-N/L, 
% % Snso_Xnso = ONtotal - Snio - Xnio; % mg-N/L, biodegradable organic nitrogen
% % Sndo = Snso_Xnso*(Sso/(Sso+Xso)); % mg-N/L, soluble biodegradable nitrogen
% % Xndo = Snso_Xnso - Sndo; % mg-N/L, particulate biodegradable nitrogen
% Xbho = 0.000000001; %0; % mgCOD/L, heterotrophic active biomass -> cant be zero, but very close to it
% Xbao = 0.000000001; %0; % mgCOD/L, autrophic active biomass -> cant be zero, but very close to it
% Soo = 0; % mgO2/L, oxygen concentration
% Xpo = 0; % mgCOD/L, biomass debris 
% Salko = ALK/100; % mM/L, alkalinity
% Snho = NH3; % Initial ammonia
% Snoo = 0; % Initial nitrite/nitrate
% % alt method (pL github)
% Sio = 0.13*CODt;
% Xio = CODio - Sio;
% Xso = 1.6*VSS - Xio;
% Sso = CODbo - Xso;
% nb_TKN = TKN*0.03;
% sol_bio_orgN_ratio = Sso/(Sso+ Xso);
% Sndo = (TKN - Snho - nb_TKN)*(sol_bio_orgN_ratio);
% Xndo = (TKN - Snho - nb_TKN)*(1 - sol_bio_orgN_ratio);

%% Intial conditions for system
% assume the initial conditions into plant are those of the reactors
% Since we arent dealing with startup, look to adjust initial conditions for each stream 
Sio  = 30; % Soluble inert organic matter
Sso  = 69.5; % Readily biodegradable substrate
Xio  = 51.2; % Particulate inert organic matter
Xso  = 202.32; % Slowly biodegradable substrate
Xbho = 28.17; % Active heterotrophic biomass
Xbao = 25; % Active autotrophic biomass
Xpo  = 0; % Particulate products arising from biomass decay
Soo  = 0; % Oxygen
Snoo = 0; % Nitrate and nitrite nitrogen
Snho = 5; % NH 4+ + NH 3 nitrogen
Sndo = 6.95; % Soluble biodegradable organic nitrogen
Xndo = 5; % Particulate biodegradable organic nitrogen
Salko = 7; % Alkalinity
MLE_int = [Sio Sso Xio Xso Xbho Xbao Xpo Soo Snoo Snho Sndo Xndo Salko]; % Create vector of initial values for the MLE system
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
sys_int = [MLE_int AD_int]; % Combine initial conditions
x = sys_int(:)*ones(1,15); % Format to an array of [components,streams]

%% Import data for varying influent flow
dat_dry = importdata('datos/Inf_dry_2006.txt','\t',1);
dat_rain = importdata('datos/Inf_rain_2006.txt','\t',1);
dat_strm = importdata('datos/Inf_strm_2006.txt','\t',1);

% Append data
drydata = dat_dry.data;
stormdata = dat_strm.data;
stormdata(:,1) = stormdata(:,1) + drydata(end,1);
InfluentData = [drydata;stormdata];
DataTime = InfluentData(:,1);

% Manipulate data for ODE input
Var.Qt = InfluentData(:,1); % Time, days
Var.Ct = Var.Qt; % Time, days
InfluentData(:,1) = []; % Remove column of time data
Var.Qflow = InfluentData(:,14); % Flow at time, Qt
InfluentData(:,14) = []; % Remove column of influent flow data
Var.C = InfluentData; % Conc at time, Ct

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
    Conc_AD = [Conc_AD;arrayManip(oa_int:oa_end,:)];
    na_int = oa_end + 1;
    na_end = na_int + 12;
    oa_int = na_end + 1;
    oa_end = na_end + 35;
    i = i + 1;
end
% Remove dummy AD columns
Conc_AD(:,1:12) = [];
% Set variables
time = Array.tArray;
[row,col] = size(Concentration);

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
% ASM1 variables that have been converted from ADM1
ASMstream.Fifteen = reshape(Concentration(:,15),[compASM,length(time)])';

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
plot(time,ASMstream.Fifteen)
title('Digester Sludge Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,15)
plot(time,Array.Qarray)
title('Plant Influent Flowrate')
ylabel('Volumetric Flow, m3/day')
xlabel('Time, days')
subplot(4,5,16)
plot(time,ASMstream.Nine(:,10))
title('Plant Effluent Ammonia')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,17)
plot(time,ASMstream.Nine(:,9))
title('Plant Effluent Nitrate/Nitrite')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,5,18)
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
S_gas_h2 = ADMstream.Two(:,1);
S_gas_ch4 = ADMstream.Two(:,2);
S_gas_co2 = ADMstream.Two(:,3);
% Determine partial pressures for each gas component, converting CH4 and H2 from kgCOD to
% kmole C (S_gas_h2/S_gas_ch4)
P_gas_h2 = S_gas_h2*l.R*l.T_op/16;
P_gas_ch4 = S_gas_ch4*l.R*l.T_op/64;
P_gas_co2 = S_gas_co2*l.R*l.T_op;

% Gas flow
P_gas = P_gas_h2 + P_gas_ch4 + P_gas_co2 + l.p_gas_h2o; % Total gas pressure
q_gas = l.k_p*(P_gas - l.P_atm).*P_gas/l.P_atm; % Total gas flow

% Calculate gas component gas flow
q_gas_h2 = P_gas_h2./P_gas.*q_gas;
q_gas_ch4 = P_gas_ch4./P_gas.*q_gas;
q_gas_co2 = P_gas_co2./P_gas.*q_gas;
q_gas_h2o = l.p_gas_h2o./P_gas.*q_gas;

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
plot(time,P_gas_ch4/P_gas,time,P_gas_co2/P_gas,time,P_gas_h2/P_gas);
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
CSI = 1; Carbon.CSI = [];
CSS = 2; Carbon.CSS = [];
CIP = 3; Carbon.CIP = [];
CSBS = 4; Carbon.CSBS = [];
CHB = 5; Carbon.CHB = [];
CAB = 6; Carbon.CAB = [];
CUB = 7; Carbon.CUB = [];
CCC = 13; Carbon.CCC = [];
% Non-ASM1 variables
CSM = 34; Carbon.CSM = [];
CCO2 = 35; Carbon.CCO2 = [];
p = 1;
while p < (length(time) + 1)
    Carbon.CSI = [Carbon.CSI;Concentration(CSI,:)];
    Carbon.CSS = [Carbon.CSS;Concentration(CSS,:)];
    Carbon.CIP = [Carbon.CIP;Concentration(CIP,:)];
    Carbon.CSBS = [Carbon.CSBS;Concentration(CSBS,:)];
    Carbon.CHB= [Carbon.CHB;Concentration(CHB,:)];
    Carbon.CAB = [Carbon.CAB;Concentration(CAB,:)];
    Carbon.CUB = [Carbon.CUB;Concentration(CUB,:)];
    Carbon.CCC = [Carbon.CCC;Concentration(CCC,:)];
    Carbon.CSM = [Carbon.CSM;Conc_AD(CSM,:)];
    Carbon.CCO2 = [Carbon.CCO2;Conc_AD(CCO2,:)];
    CSI = CCC + 1;
    CSS = CCC + 2;
    CIP = CCC + 3;
    CSBS = CCC + 4;
    CHB = CCC + 5;
    CAB = CCC + 6;
    CUB = CCC + 7;
    CCC = CCC + 13;
    CSM = CCO2 + 34;
    CCO2 = CCO2 + 35;
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
Carbon.CSMkgram = Carbon.CSM.*q_gas.*1000; % Convert to grams
Carbon.CCO2kmoleC = Carbon.CCO2.*q_gas.*1000; % Convert to mol

% Determine total carbon mass flow rate [gC/day] in stream m
% Multiply carbon fraction by the corresponding carbon component
% Methanol still needs to be implemented
m = 1;
C = ones(length(time),col); % Preallocate array
while m < (col + 1)
    C(:,m) = Carbon.CSIgrams(:,m).*CF.C_S_I + Carbon.CIPgrams(:,m).*CF.C_I_P + Carbon.CSBSgrams(:,m).*CF.C_SB_S +... 
    Carbon.CSSgrams(:,m).*CF.C_S_S + Carbon.CCCgrams(:,m).*CF.C_CC+ Carbon.CHBgrams(:,m).*CF.C_HB + Carbon.CABgrams(:,m).*CF.C_AB +...
    Carbon.CUBgrams(:,m).*CF.C_UB;
    m = m + 1;
end
% Adjust for gas stream
% Convert from COD to grams and mol to grams
C(:,14) = Carbon.CSMkgram(:,2).*CF.C_S_M + Carbon.CCO2kmoleC(:,2).*12.01; 
% Convert grams to lbs 
C = 0.00220462.*C;
% Mass balance check
PerError = [];
PerError(:,1) = 100.*abs(C(:,1) - (C(:,2) + C(:,3)))./C(:,1);
PerError(:,2) = 100.*abs(C(:,4) - (C(:,12) + C(:,8)+ C(:,2)))./C(:,4);
PerError(:,3) = 100.*abs(C(:,4) - C(:,5))./C(:,4);
PerError(:,4) = 100.*abs(C(:,5) - C(:,6))./C(:,5);
PerError(:,5) = 100.*abs(C(:,6) - (C(:,8) + C(:,7)))./C(:,6);
PerError(:,6) = 100.*abs(C(:,7) - (C(:,9) + C(:,10)))./C(:,7);
PerError(:,7) = 100.*abs(C(:,10) - (C(:,12) + C(:,11)))./C(:,10);
PerError(:,8) = 100.*abs(C(:,1) - (C(:,11) + C(:,9) + C(:,3)))./C(:,1);
PerError(:,9) = 100.*abs(C(:,13) - (C(:,11) + C(:,3)))./C(:,13);
% Ignore Carbon balance on AD for now
%PerError(:,10) = 100.*abs(C(:,13) - C(:,15) - C(:,14))./C(:,13);
% Plot percent error
figure(3)
plot(time,PerError)
title('Mass Balance Percent Error')
ylabel('Error, %')
xlabel('Time, days')

% Plot carbon flow for each stream separately
figure(4)
subplot(4,5,1)
plot(time,C(:,1))
title('Plant Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,2)
plot(time,C(:,2))
title('Primary Clarifier Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,3)
plot(time,C(:,3))
title('Primary Clarifier WAS')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,4)
plot(time,C(:,4))
title('Anoxic Tank Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,5)
plot(time,C(:,5))
title('Anoxic Tank Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,6)
plot(time,C(:,6))
title('Aeration Tank Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,7)
plot(time,C(:,7))
title('SC Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,8)
plot(time,C(:,8))
title('Internal Recycle')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,9)
plot(time,C(:,9))
title('SC Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,10)
plot(time,C(:,10))
title('SC Underflow')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,11)
plot(time,C(:,11))
title('WAS Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,5,12)
plot(time,C(:,12))
title('RAS Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
% Ignore carbon balance on AD for now
% subplot(4,5,13)
% plot(time,C(:,13))
% title('Mixing of WAS and PS')
% ylabel('Carbon Flow, lb/day')
% xlabel('Time, days')
% subplot(4,5,14)
% plot(time,C(:,14))
% title('Anaerobic Digester Gas Stream')
% ylabel('Carbon Flow, lb/day')
% xlabel('Time, days')
% subplot(4,5,15)
% plot(time,C(:,15))
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
CompASM = 13;
CompADM = 35;
%% Dynamic flow
% Constant Plant flow - > testing average data -> using gal/min going into
% NT converted to m3/day
%Qplant = 57622.44668;
Qplant = interp1(Var.Qt,Var.Qflow,t); % Interpolate data set of volumetric flow at specified time
%% Flow solver
opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
Qint = [20000 19000 100 60000 60000 60000 25000 17000 10000 10000 500 10000 600 0 600];
[Q,~] = fsolve(@(Q) myfunc(Q,Qplant,Var),Qint,opts);
function Qsolve = myfunc(Q,Qplant,Var)
%% Solves flow balance for given plant flow, recycle/waste fractions, and clarifier splits
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
Qsolve = [Qplant - Q(1);
Q(1) - Q(3) - Q(2);
Var.Rir*Q(2) - Q(8);
Var.Rr*Q(2) - Q(12);
Var.fpc*Q(1) - Q(3);
Var.fsc*Q(10) - Q(12);
Q(2) + Q(8) + Q(12) - Q(4);
Q(5) - Q(4);
Q(6) - Q(5);
Q(7) + Q(8) - Q(6);
Q(9) + Q(10) - Q(7);
Q(12) + Q(11) - Q(10);
Q(3) + Q(9) + Q(11) - Q(1);
Q(3) + Q(11) - Q(13);
Q(13) - Q(15);
];
end
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
Vol6 = Var.param(27); % AD volume


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

%% ODE for MLE system
% Concentration of component i in stream 1-12
% K(1,i)*theta1(1) + K(2,i)*theta1(2) + K(3,i)*theta1(3) + K(4,i)*theta1(4) + K(5,i)*theta1(5) + K(6,i)*theta1(6) + K(7,i)*theta1(7) + K(8,i)*theta1(8) 
i = 1;
Conc = zeros(length(dCdt),numel(Q));
while i < (CompASM + 1)
    % Model won't be used for my plant simulation
    %% Modeling primary clarifier - Otterpohl and Freund 1992
%     hrt = Vol1/Q(1); % Hydraulic residence time
%     n_COD = 2.7*log(hrt*hrt + 9)/100; % Removal efficiency
    XCOD1 = dCdt(3,1) + dCdt(4,1) + dCdt(5,1) + dCdt(6,1) + dCdt(7,1); % Particulate COD in influent
%     CODin = dCdt(1,1) + dCdt(2,1) + XCOD1; % Total COD in influent
%     n_x = (n_COD*CODin)/XCOD1;
%     if n_x > 0.95
%        n_x = 0.95;
%     elseif n_x < 0.05
%        n_x = 0.05;
%     else
%        n_x = n_x;
%     end
%     n_x = (1-n_x);
    %% TSS removal
    n_x = 0.533; % Fraction of TSS left in effluent, taken as average from ST and NT from Appendix B GPS-X files from B&V
    % Determine which components are separated
    % Comment out the two lines below to run constant influent data
    Dyn_conc = interp1(Var.Ct,Var.C,t); % Interpolate data set of concentration at specified time
    dCdt(i,1) = Dyn_conc(:,i);
    
    if i < 3
        dCdt(i,2) = dCdt(i,1);
    elseif (2<i) && (i<8)
        dCdt(i,2) = n_x*dCdt(i,1)*Q(1)/Q(2);
    elseif (7<i) && (i<12)
        dCdt(i,2) = dCdt(i,1); 
    elseif (11<i) && (i<13)
        dCdt(i,2) = n_x*dCdt(i,1)*Q(1)/Q(2);
    else
        dCdt(i,2) = dCdt(i,1);
    end
    if i < 3
        dCdt(i,3) = dCdt(i,1);
    elseif (2<i) && (i<8)
        dCdt(i,3) = (1-n_x)*dCdt(i,1)*Q(1)/Q(3);
    elseif (7<i) && (i<12)
        dCdt(i,3) = dCdt(i,1);
    elseif (11<i) && (i<13)
        dCdt(i,3) = (1-n_x)*dCdt(i,1)*Q(1)/Q(3);
    else
        dCdt(i,3) = dCdt(i,1);
    end
    dCdt(i,2) = (dCdt(i,1)*Q(1) - dCdt(i,3)*Q(3))/Q(2); % Mass balance for flow into/out of Primary Clarifier

    %% Recycle split, same concentration
    mass = zeros(13,12);
    mass(i,2) = Q(2)*dCdt(i,2);
    mass(i,3) = Q(3)*dCdt(i,3);
    mass(i,6) = Q(6)*dCdt(i,6);
    %disp(mass)
    m_rat_rec = Q(7)/Q(6); % Mass ratio
    dCdt(i,7) = m_rat_rec*mass(i,6)/Q(7); % Recycle Split
    dCdt(i,8) = (1-m_rat_rec)*mass(i,6)/Q(8); % Recycle Split

    % Model not worth implementing, values taken from B&V will be used from
    % App. B from GPS-X pdf file
    %% Secondary clarifier model
    c_x = 0.0015; % Fraction of TSS left in effluent, taken as average from ST and NT of B&V App. B file
    XCOD7 = dCdt(3,7) + dCdt(4,7) + dCdt(5,7) + dCdt(6,7) + dCdt(7,7); % Particulate COD in influent
    if i < 3
        dCdt(i,9) = dCdt(i,7);
    elseif (2<i) && (i<8)
        dCdt(i,9) = c_x*dCdt(i,7)*Q(7)/Q(9);
    elseif (7<i) && (i<12)
        dCdt(i,9) = dCdt(i,7);
    elseif (11<i) && (i<13)
        dCdt(i,9) = c_x*dCdt(i,7)*Q(7)/Q(9);
    else
        dCdt(i,9) = dCdt(i,7);
    end
    if i < 3
        dCdt(i,10) = dCdt(i,7);
    elseif (2<i) && (i<8)
        dCdt(i,10) = (1-c_x)*dCdt(i,7)*Q(7)/Q(10);
    elseif (7<i) && (i<12)
        dCdt(i,10) = dCdt(i,7);
    elseif (11<i) && (i<13)
        dCdt(i,10) = (1-c_x)*dCdt(i,7)*Q(7)/Q(10);
    else
        dCdt(i,10) = dCdt(i,7);
    end
    
    dCdt(i,7) = (dCdt(i,9)*Q(9) + dCdt(i,10)*Q(10))/Q(7); % Mass balance for flow into/out of Primary Clarifier
    XCOD9 = dCdt(3,9) + dCdt(4,9) + dCdt(5,9) + dCdt(6,9) + dCdt(7,9); % Particulate COD in effluent
    XCOD10 = dCdt(3,10) + dCdt(4,10) + dCdt(5,10) + dCdt(6,10) + dCdt(7,10); % Particulate COD in waste stream
    
    %% Waste split
    mass(i,10) = Q(10)*dCdt(i,10);
    m_rat_SS = Q(11)/Q(10); % Mass ratio
    dCdt(i,11) = m_rat_SS*mass(i,10)/Q(11); % WAS split
    dCdt(i,12) = (1-m_rat_SS)*mass(i,10)/Q(12); % RAS split
    
    %% Anox/Aer
    dCdt(i,4) = (Q(2)*dCdt(i,2) + Q(8)*dCdt(i,8) + Q(12)*dCdt(i,12))/Q(4); % mixing point
    
    Conc(i,5) = 1/Vol2*(Q(4)*dCdt(i,4) - Q(5)*dCdt(i,5)) + K(1,i)*theta1(1) + K(2,i)*theta1(2) + K(3,i)*theta1(3) + K(4,i)*theta1(4) + K(5,i)*theta1(5) + K(6,i)*theta1(6) + K(7,i)*theta1(7) + K(8,i)*theta1(8); % Anoxic balance
    Conc(i,6) = 1/Vol3*(Q(5)*dCdt(i,5) - Q(6)*dCdt(i,6)) + K(1,i)*theta2(1) + K(2,i)*theta2(2) + K(3,i)*theta2(3) + K(4,i)*theta2(4) + K(5,i)*theta2(5) + K(6,i)*theta2(6) + K(7,i)*theta2(7) + K(8,i)*theta2(8); % Aeration general balance
    if i == 8
        % Rough control of DO in aeration zone
        tadj = round(t*100)/100;
        t_find = find(tadj == Var.tsample);
        dCdt(10,6);
        dCdt(8,6);
        if t == 1
            KLa = Var.param(20); % Oxygen transfer coefficient
        else
        end
%         if numel(t_find) > 0
%             if dCdt(10,6) > 5
%                 KLa = Var.param(20);
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
%     XCOD6 = dCdt(3,6) + dCdt(4,6) + dCdt(5,6) + dCdt(6,6) + dCdt(7,6); % Particulate COD in MLSS
%     XCOD11 = dCdt(3,11) + dCdt(4,11) + dCdt(5,11) + dCdt(6,11) + dCdt(7,11); % Particulate COD in WAS
%     MCRT = (XCOD6*Vol3)/(XCOD11*Q(11))
%     assignin('base','MCRT_step',MCRT);
%     evalin('base','MCRT_out(end+1) = MCRT_step;');
%% Conversion from ASM1 to ADM1
    % Using paper: TOWARDS AN ASM1 - ADM1 STATE VARIABLE INTERFACE FOR
    % PLANT-WIDE WASTEWATER TREATMENT MODELING
    %% Waste sludge mixing
    dCdt(i,13) = (dCdt(i,11)*Q(11) + dCdt(i,3)*Q(3))/Q(13);
    
    %% Reducing total incoming COD for Ss,Xs,Xbh,Xba in that specific order
    % Maybe optimize for loop
    for lenComp = 1:13
        xtemp(lenComp) = dCdt(lenComp,13);
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
 CODdemand = dCdt(8,13) + CODequiv*dCdt(9,13);
     if CODdemand > (dCdt(2,13) + dCdt(4,13) + dCdt(5,13) + dCdt(6,13))
        disp('Warning! Influent characterization may need to be evaluated, not enough COD available')
     else
     end
 xtemp(8) = 0;
 xtemp(9) = 0;
    if CODdemand > dCdt(2,13)
        remaina = CODdemand - dCdt(2,13);
        xtemp(2) = 0;
        if remaina > dCdt(4,13)
            remainb = remaina - dCdt(4,13);
            xtemp(4) = 0;
                if remainb > dCdt(5,13)
                    remainc = remainb - dCdt(5,13);
                    xtemp(10) = xtemp(10) + dCdt(5,13)*fnbac;
                    xtemp(5) = 0;
                        if remainc > dCdt(6,13)
                            remaind = remainc - dCdt(6,13);
                            xtemp(10) = xtemp(10) + dCdt(6,13)*fnbac;
                            xtemp(6) = 0;
                            xtemp(8) = remaind;
                            disp('ERROR: COD shortage when removing inital oxygen and nitrate')
                        else
                            xtemp(6) = dCdt(6,13) - remainc;
                            xtemp(10) = xtemp(10) + remainc*fnbac;
                        end
                else
                    xtemp(5) = dCdt(5,13) - remainb;
                    xtemp(10) = xtemp(10) + remainb*fnbac;
                end
        else
            xtemp(4) = dCdt(4,13) - remaina;
        end
    else
        xtemp(2) = dCdt(2,13) - CODdemand;
    end
    
    sorgn = dCdt(11,13)/fnaa;
    if sorgn >= xtemp(2)
        adm(2) = xtemp(2);
        xtemp(11) = xtemp(11) - xtemp(2)*fnaa;
        xtemp(2) = 0;
    else
        adm(2) = sorgn;
        xtemp(2) = xtemp(2) - sorgn;
        xtemp(11) = 0;
    end
    
    xorgn = dCdt(12,13)/fnaa;
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
                    inertS = inertS + temp(10)/fsni_adm;
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
    chargeasm1 = (dCdt(13,13)*alfachalk + dCdt(10,13)*alfachnh +...
        dCdt(9,13)*alfachno)/1000;
    chargeadm1 = adm(7)*alfachac + adm(6)*alfachpro +...
        adm(5)*alfachbu + adm(4)*alfachva + adm(11)*alfachin;
    adm(10) = (chargeasm1 - chargeadm1)/alfachic;
    pkk = 14 - log10(exp(55900/bigr/100*tfac));
    ancat = chargeadm1 + adm(10)*alfachic - 10^(-ph) + 10^(-pkk + ph);
    
    % Check mass balances
    totCODin = 0;
    for vec = 1:7
        totCODin = totCODin + dCdt(vec,13);
    end
    
    totNin = dCdt(9,13) + dCdt(10,13) + dCdt(11,13) + dCdt(12,13)...
        + (nbac/1000)*(dCdt(5,13) + dCdt(6,13)) +  (sni/1000)*dCdt(1,13)...
        + (xni/1000)*(dCdt(3,13) + dCdt(7,13));
    
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
    MBCOD = totCODin - totCODout;
    MBN = totNin - totNout;
    
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
    % Test variables to see if AD works against original paper
%     S_su_in = 0; %l.S_su_in;
%     S_aa_in = 0.043902; %l.S_aa_in;
%     S_fa_in = 0; %l.S_fa_in;
%     S_va_in = 0; %l.S_va_in;
%     S_bu_in = 0; %l.S_bu_in;
%     S_pro_in = 0; %l.S_pro_in;
%     S_ac_in = 0; %l.S_ac_in;
%     S_h2_in = 0; %l.S_h2_in;
%     S_ch4_in = 0; %l.S_ch4_in;
%     S_IC_in = 0.007918; %l.S_IC_in;
%     S_IN_in = 0.001972; %l.S_IN_in;
%     S_I_in = 0.028067; %l.S_I_in;
%     X_c_in = 0; %l.X_xc_in;
%     X_ch_in = 3.7236; %l.X_ch_in;
%     X_li_in = 8.04686; %l.X_li_in;
%     X_I_in = 17.011; %l.X_I_in;
%     S_cat_in = 0; %l.S_cat_in;
%     S_an_in = 0.0052; %l.S_an_in;
%     X_pr_in = 15.9243; %l.X_pr_in;
%     X_su_in = 0; %l.X_su_in;
%     X_aa_in = 0; %l.X_aa_in;
%     X_fa_in = 0; %l.X_fa_in;
%     X_c4_in = 0; %l.X_c4_in;
%     X_pro_in = 0; %l.X_pro_in;
%     X_ac_in = 0; %l.X_ac_in;
%     X_h2_in = 0; %l.X_h2_in;
    % Save to concentration
    dCdt(14,13) = S_su_in;
    dCdt(15,13) = S_aa_in;
    dCdt(16,13) = S_fa_in;
    dCdt(17,13) = S_va_in;
    dCdt(18,13) = S_bu_in;
    dCdt(19,13) = S_pro_in;
    dCdt(20,13) = S_ac_in;
    dCdt(21,13) = S_h2_in;
    dCdt(22,13) = S_ch4_in;
    dCdt(23,13) = S_IC_in;
    dCdt(24,13) = S_IN_in;
    dCdt(25,13) = S_I_in;
    dCdt(26,13) = X_c_in;
    dCdt(27,13) = X_ch_in;
    dCdt(28,13) = X_pr_in;
    dCdt(29,13) = X_li_in;
    dCdt(30,13) = X_su_in;
    dCdt(31,13) = X_aa_in;
    dCdt(32,13) = X_fa_in;
    dCdt(33,13) = X_c4_in;
    dCdt(34,13) = X_pro_in;
    dCdt(35,13) = X_ac_in;
    dCdt(36,13) = X_h2_in;
    dCdt(37,13) = X_I_in;
    dCdt(38,13) = S_cat_in;
    dCdt(39,13) = S_an_in;
    % Variable change to check against AD paper
    %l.V_liq = 0.9189*Var.param(27);
    %l.V_gas = 0.0811*Var.param(27);
    %l.q_in = Q(13);
    
%% AD differential equations
% Data/Equations pulled from Aspects on ADM1 implementation within the BSM2 framework
% CHECK PARAMETERS -> CHANGE VOLUMES (headspace/liquid)
    % Intial conditions inside reactor
    S_su = dCdt(14,15);
    S_aa = dCdt(15,15);
    S_fa = dCdt(16,15);
    S_va = dCdt(17,15);
    S_bu = dCdt(18,15);
    S_pro = dCdt(19,15);
    S_ac = dCdt(20,15);
    S_h2 = dCdt(21,15);
    S_ch4 = dCdt(22,15);
    S_IC = dCdt(23,15);
    S_IN = dCdt(24,15);
    S_I = dCdt(25,15);
    X_c = dCdt(26,15);
    X_ch = dCdt(27,15);
    X_pr = dCdt(28,15);
    X_li = dCdt(29,15);
    X_su = dCdt(30,15);
    X_aa = dCdt(31,15);
    X_fa = dCdt(32,15);
    X_c4 = dCdt(33,15);
    X_pro = dCdt(34,15);
    X_ac = dCdt(35,15);
    X_h2 = dCdt(36,15);
    X_I = dCdt(37,15);
    S_cat = dCdt(38,15);
    S_an = dCdt(39,15);
    S_vam = dCdt(40,15);
    S_bum = dCdt(41,15);
    S_prom = dCdt(42,15);
    S_acm = dCdt(43,15);
    S_hco3m = dCdt(44,15);
    S_nh3 = dCdt(45,15);
    S_gas_h2 = dCdt(46,14);
    S_gas_ch4 = dCdt(47,14);
    S_gas_co2 = dCdt(48,14);
    
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
    Conc(14,15) = Q(13)/l.V_liq*(S_su_in - S_su) + rho_2 + (1 - l.f_fa_li)*rho_4 - rho_5;       % 1
    Conc(15,15) = Q(13)/l.V_liq*(S_aa_in - S_aa) + rho_3 - rho_6;                               % 2
    Conc(16,15) = Q(13)/l.V_liq*(S_fa_in - S_fa) + l.f_fa_li*rho_4 - rho_7;                     % 3
    Conc(17,15) = Q(13)/l.V_liq*(S_va_in - S_va) + (1 - l.Y_aa)*l.f_va_aa*rho_6 - rho_8;        % 4
    
    % Differential equations 5-8, soluble matter
    Conc(18,15) = Q(13)/l.V_liq*(S_bu_in - S_bu) + (1 - l.Y_su)*l.f_bu_su*rho_5 + ...           % 5
    (1 - l.Y_aa)*l.f_bu_aa*rho_6 - rho_9; 
    Conc(19,15) = Q(13)/l.V_liq*(S_pro_in - S_pro) + (1 - l.Y_su)*l.f_pro_su*rho_5 + ...          % 6 
    (1 - l.Y_aa)*l.f_pro_aa*rho_6 + (1-l.Y_c4)*0.54*rho_8 - rho_10;
    Conc(20,15) = Q(13)/l.V_liq*(S_ac_in - S_ac) + (1 - l.Y_su)*l.f_ac_su*rho_5 + ...             % 7
    (1 - l.Y_aa)*l.f_ac_aa*rho_6 + (1 - l.Y_fa)*0.7*rho_7 + ... 
    (1 - l.Y_c4)*0.31*rho_8 + (1 - l.Y_c4)*0.8*rho_9 + ... 
    (1 - l.Y_pro)*0.57*rho_10 - rho_11;
    Conc(21,15) = Q(13)/l.V_liq*(S_h2_in - S_h2) + (1 - l.Y_su)*l.f_h2_su*rho_5 + ...           % 8
    (1 - l.Y_aa)*l.f_h2_aa*rho_6 + (1 - l.Y_fa)*0.3*rho_7 + ... 
    (1 - l.Y_c4)*0.15*rho_8 + (1 - l.Y_c4)*0.2*rho_9 + ... 
    (1 - l.Y_pro)*0.43*rho_10 - rho_12 - rho_T_8;

    % Differential equations 9-12, soluble matter
    Conc(22,15) = Q(13)/l.V_liq*(S_ch4_in - S_ch4) + (1 - l.Y_ac)*rho_11 + ...                  % 9
    (1 - l.Y_h2)*rho_12 - rho_T_9; 
    Conc(23,15) = Q(13)/l.V_liq*(S_IC_in - S_IC) - (s_1*rho_1 + s_2*rho_2 + ...                 % 10
    s_3*rho_3 + s_4*rho_4 + s_5*rho_5 + s_6*rho_6 + s_7*rho_7 + ... 
    s_8*rho_8 + s_9*rho_9 + s_10*rho_10 + s_11*rho_11 + s_12*rho_12 + ... 
    s_13*(rho_13 + rho_14 + rho_15 + rho_16... 
    + rho_17 + rho_18 + rho_19)) - rho_T_10; 
    Conc(24,15) = Q(13)/l.V_liq*(S_IN_in - S_IN) - l.Y_su*l.N_bac*rho_5 + ...                   % 11
    (l.N_aa - l.Y_aa*l.N_bac)*rho_6 - l.Y_fa*l.N_bac*rho_7 - ... 
    l.Y_c4*l.N_bac*rho_8 - l.Y_c4*l.N_bac*rho_9 - ... 
    l.Y_pro*l.N_bac*rho_10 - l.Y_ac*l.N_bac*rho_11 - ... 
    l.Y_h2*l.N_bac*rho_12 + (l.N_bac - l.N_xc)*... 
    (rho_13 + rho_14 + rho_15 + rho_16 + rho_17 + rho_18 +... 
    rho_19) + (l.N_xc - l.f_xI_xc*l.N_I - l.f_sI_xc*l.N_I - l.f_pr_xc*l.N_aa)*rho_1; 
    Conc(25,15) = Q(13)/l.V_liq*(S_I_in - S_I) + l.f_sI_xc*rho_1;                               % 12
    
    % Differential equations 13-16, particulate matter
    Conc(26,15) = Q(13)/l.V_liq*(X_c_in - X_c) - rho_1 + ...                                   % 13
    rho_13 + rho_14 + rho_15 + rho_16 + rho_17 + rho_18 + rho_19; 
    Conc(27,15) = Q(13)/l.V_liq*(X_ch_in - X_ch) + l.f_ch_xc*rho_1 - rho_2;                     % 14
    Conc(28,15) = Q(13)/l.V_liq*(X_pr_in - X_pr) + l.f_pr_xc*rho_1 - rho_3;                     % 15
    Conc(29,15) = Q(13)/l.V_liq*(X_li_in - X_li) + l.f_li_xc*rho_1 - rho_4;                     % 16
    
    % Differential equations 17-20, particulate matter
    Conc(30,15) = Q(13)/l.V_liq*(X_su_in - X_su) + l.Y_su*rho_5 - rho_13;                       % 17
    Conc(31,15) = Q(13)/l.V_liq*(X_aa_in - X_aa) + l.Y_aa*rho_6 - rho_14;                       % 18
    Conc(32,15) = Q(13)/l.V_liq*(X_fa_in - X_fa) + l.Y_fa*rho_7 - rho_15;                       % 19
    Conc(33,15)  = Q(13)/l.V_liq*(X_c4_in - X_c4) + l.Y_c4*rho_8 + l.Y_c4*rho_9 - rho_16;       % 20
    
    % Differential equations 21-24, particulate matter
    Conc(34,15) = Q(13)/l.V_liq*(X_pro_in - X_pro) + l.Y_pro*rho_10 - rho_17;                   % 21
    Conc(35,15) = Q(13)/l.V_liq*(X_ac_in - X_ac) + l.Y_ac*rho_11 - rho_18;                      % 22
    Conc(36,15) = Q(13)/l.V_liq*(X_h2_in - X_h2) + l.Y_h2*rho_12 - rho_19;                      % 23
    Conc(37,15) = Q(13)/l.V_liq*(X_I_in - X_I) + l.f_xI_xc*rho_1;                               % 24
    
    % Differential equations 25-26, cations and anions
    Conc(38,15) = Q(13)/l.V_liq*(S_cat_in - S_cat);                                             % 25
    Conc(39,15) = Q(13)/l.V_liq*(S_an_in - S_an);                                               % 26
    
    % Differential equations 27-32, ion states
    Conc(40,15) = -rho_A_4;                                                                     % 27
    Conc(41,15) = -rho_A_5;                                                                     % 28
    Conc(42,15) = -rho_A_6;                                                                     % 29
    Conc(43,15) = -rho_A_7;                                                                     % 30
    Conc(44,15) = -rho_A_10;                                                                    % 31
    Conc(45,15) = -rho_A_11;                                                                    % 32
    
    % Differential equations 33-35, gas phase equations
    Conc(46,14) = -S_gas_h2*q_gas/l.V_gas + rho_T_8*l.V_liq/l.V_gas;                            % 33
    Conc(47,14) = -S_gas_ch4*q_gas/l.V_gas + rho_T_9*l.V_liq/l.V_gas;                           % 34
    Conc(48,14) = -S_gas_co2*q_gas/l.V_gas + rho_T_10*l.V_liq/l.V_gas;                          % 35
    
%% Conversion from ADM1 to ASM1
% Use fortran code
    for vec3 = 1:13
        asmm(vec3) = 0;
    end
    for vec4 = 1:24
        xtemp(vec4) = dCdt(13 + vec4,15);
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
        xptemp = biomass*fnbacfxni;
        biomass_bioN = 0;
    else
        xptemp = biomass_nobio;
    end
    if ((biomass_bioN/fnxc) <= (biomass - biomass_nobio))
        xstemp = biomass - biomass_nobio - xstemp;
        if ((xtemp(11)*14000/fnxc) > remainCOD)
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
    if (fsni_adm <fsni)
        inertS = xtemp(12)*fsni_adm.fsni;
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
        adm(11)*14000 + sni_adm*adm(12) + xni*adm(15);
    totCODout = 0;
    for vec11 = 1:7
        totCODout = totCODout + asmm(vec11);
    end
    % SI_N not included here below
    totTKNout = asmm(10) + asmm(11) + asmm(12) + fsni*asmm(1) +...
        fnbac*(asmm(5) + asmm(6)) + fxni*(asmm(3) + asmm(7));
    MBCOD2 = totCODin - totCODout;
    MBTNK = totTKNin - totTKNout;

    % CONVERT UNITS CORRECTLY BACK TO ASM1
    % CODconserved = CODt_anaerobic - Sh2 - Sch4;
    % COD Conversions
    dCdt(1,15) = asmm(1);
    dCdt(2,15) = asmm(2);
    dCdt(3,15) = asmm(3);
    dCdt(4,15) = asmm(4);
    dCdt(5,15) = asmm(5);
    dCdt(6,15) = asmm(6);
    dCdt(7,15) = asmm(7);
    dCdt(8,15) = asmm(8);
    dCdt(9,15) = asmm(9);
    dCdt(10,15) = asmm(10);
    dCdt(11,15) = asmm(11);
    dCdt(12,15) = asmm(12);
    dCdt(13,15) = asmm(13);
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
elseif t >= 1.05
    % Avoid large dataset and quicken solver time by only aquiring data at the specified time span
    adjT = round(t*100)/100;
    findT = find(adjT == Var.timespan);
    if numel(findT) > 0
        Array.mleArray = [Array.mleArray;dCdt];
        Array.tArray = [Array.tArray;t];
        Array.Qarray = [Array.Qarray;Q];
        assignin('base','Array',Array);
    else
    end
else
end
fprintf('Current simulation time is %6.6f days\n',t)
end