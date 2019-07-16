tic
%% TO-DO
% Adjust initial values for each stream throughout the entire plant
% Rename streams to actual processes
% Convert ASM values to ADM
    % UNITS NEED TO BE ADJUSTED
    % Essentially done
    % Need to verify and find missing variables
    
% Implement thickener(before AD)/dewatering(after AD)/dryer(after AD)
% Implement denitrification filter

% Create PID for oxygen transfer in aeration zone
% Create GUI for kinetics/ recycle ratios/ simulation time/ flow variables

%% Issues
% AD isn't working properly
% ***** CHECK UNITS *****

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
l = SE_IndataADM1_v2;

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
% % Units g/m3 = mg/L
% TSS = 82; % mg/L
% VSS = 61.5; % mg/L
% If no VSS data, assume TSS/VSS ratio of 0.69 or 0.75
% tv = 0.69;
% VSS = TSS/tv;
% FSS = TSS - VSS; % fixed suspended solids
% BOD5 = 155; % mg/L
% CODtotal = 325; % mg/L
% TKN = 43.5; % mg-N/L
% NH3 = 25; % mg-N/L
% NO3 = 0; % mg-N/L
% ALK = 200; % mg-N/L
% % Conversion to ASM1 Variables
% CODt = 2.1*BOD5; % mgCOD/L, Converts BOD5 to total COD, if not available
% CODbo = 1.71*BOD5; % mgCOD/L, biodegradable COD
% CODio = CODt - CODbo; % mgCOD/L, inert COD
% Xio = 0.375*1.5*VSS; % mgCOD/L, particulate inert COD
% Sio = CODio - Xio; % mgCOD/L, soluble inert COD
% f_readily = 0.43; % fraction of biodegradable COD that is readily biodegradable
% Sso = CODbo*f_readily; % mgCOD/L, readily biodegradable substrate
% Xso = CODbo - Sso; % mgCOD/L, slowly biodegradable substrate
% ONtotal = TKN - NH3; % mg-N/L, total organic nitrogen
% Snio = 1.5; % mg-N/L, soluble inert organic nitrogen
% in_xd = 0.06; % mass of nitrogen per mass of COD in biomass
% Xnio = in_xd*Xio; % mg-N/L, 
% Snso_Xnso = ONtotal - Snio - Xnio; % mg-N/L, biodegradable organic nitrogen
% Sndo = Snso_Xnso*(Sso/(Sso+Xso)); % mg-N/L, soluble biodegradable nitrogen
% Xndo = Snso_Xnso - Sndo; % mg-N/L, particulate biodegradable nitrogen
% Xbho = 0.000000001; %0; % mgCOD/L, heterotrophic active biomass -> cant be zero, but very close to it
% Xbao = 0.000000001; %0; % mgCOD/L, autrophic active biomass -> cant be zero, but very close to it
% Soo = 0; % mgO2/L, oxygen concentration
% Xpo = 0; % mgCOD/L, biomass debris 
% Salko = ALK/100; % mM/L, alkalinity
% Snho = NH3; % Initial ammonia
% Snoo = 0; % Initial nitrite/nitrate
% alt method (poopyLab github)
%nb_TKN = TKN*0.03;
%sol_bio_orgN_ratio = Sso/(Sso+ Xso);
%inf_S_NS = (TKN - Snho - nb_TKN)*(sol_bio_orgN_ratio);
%inf_X_NS = (TKN - Snho - nb_TKN)*(1 - sol_bio_orgN_ratio);

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
dat_dry=importdata('datos/Inf_dry_2006.txt','\t',1);
dat_rain=importdata('datos/Inf_rain_2006.txt','\t',1);
dat_strm=importdata('datos/Inf_strm_2006.txt','\t',1);

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

%% Carbon content fractions
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

%% Carbon (grams) in system
% Create arrays for carbon
% Each carbon component relates to a component in the ASM1 Concentration
% data set, which will later be multiplied by the carbon content fraction
CSI = 1; Carbon.CSI = [];
CSS = 2; Carbon.CSS = [];
CIP = 3; Carbon.CIP = [];
CSBS = 4; Carbon.CSBS = [];
CHB = 5; Carbon.CHB = [];
CAB = 6; Carbon.CAB = [];
CUB = 7; Carbon.CUB = [];
CCC = 13; Carbon.CCC = [];
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
    CSI = CCC + 1;
    CSS = CCC + 2;
    CIP = CCC + 3;
    CSBS = CCC + 4;
    CHB = CCC + 5;
    CAB = CCC + 6;
    CUB = CCC + 7;
    CCC = CCC + 13;
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

% Determine total carbon mass flow rate [gC/day] in stream m
% Multiply carbon fraction by the corresponding carbon component
% Methane, methanol, and CO2 still need to be implemented
m = 1;
C = ones(length(time),col); % Preallocate array
while m < (col + 1)
    C(:,m) = Carbon.CSIgrams(:,m).*CF.C_S_I + Carbon.CIPgrams(:,m).*CF.C_I_P + Carbon.CSBSgrams(:,m).*CF.C_SB_S +... 
    Carbon.CSSgrams(:,m).*CF.C_S_S + Carbon.CCCgrams(:,m).*CF.C_CC+ Carbon.CHBgrams(:,m).*CF.C_HB + Carbon.CABgrams(:,m).*CF.C_AB +...
    Carbon.CUBgrams(:,m).*CF.C_UB;
    m = m + 1;
end
% Mass balance check
error(:,1) = C(:,1) - (C(:,2) + C(:,3));
error(:,2) = C(:,4) - (C(:,12) + C(:,8)+ C(:,2));
error(:,3) = C(:,4) - C(:,5);
error(:,4) = C(:,5) - C(:,6);
error(:,5) = C(:,6) - (C(:,8) + C(:,7));
error(:,6) = C(:,7) - (C(:,9) + C(:,10));
error(:,7) = C(:,10) - (C(:,12) + C(:,11));
error(:,8) = C(:,1) - (C(:,11) + C(:,9) + C(:,3));
error(:,9) = C(:,13) - (C(:,11) + C(:,3));
error(:,10) = C(:,15) - C(:,15);
% Plot error and Carbon flow
figure(1)
plot(time,error)
title('Mass Balance Error')
ylabel('Error')
xlabel('Time, days')
figure(2)
plot(time,C)
title('Carbon Flow Throughout System')
ylabel('Carbon Mass Flowrate, g/day')
xlabel('Time, days')

%% Plotting separate results
j = 1;
stream = [];
figure(3)
% Plot results
while j < (17+1)
    if j < (col + 1)
    rev = num2str(j);
    stream.rev = reshape(Concentration(:,j),[compASM,length(time)])';
    subplot(4,5,j);
    plot(time,stream.rev);
    title(['Stream ' rev])
    ylabel('Concentration, mg/L')
    xlabel('Time, days')
    else
    end
    if j == 6
        subplot(4,5,16)
        plot(time,stream.rev(:,8))
        title('DO')
        ylabel('Concentration, mg/L')
        xlabel('Time, days')
    elseif j == 9
        subplot(4,5,17)
        plot(time,stream.rev(:,10))
        title('Effluent Ammonia')
        ylabel('Concentration, mg/L')
        xlabel('Time, days')
        subplot(4,5,18)
        plot(time,stream.rev(:,9))
        title('Effluent Nitrate/Nitrite')
        ylabel('Concentration, mg/L')
        xlabel('Time, days')
    elseif j == 19
        subplot(4,5,j)
        plot(time,Array.Qarray)
        title('Plant flow')
        ylabel('Volumetric Flow, m3/day')
        xlabel('Time, days')
    else
    end
    j = j + 1;
end

AD_14 = reshape(Conc_AD(:,2),[compADM,length(time)])';
gasAD = AD_14(:,33:35);
S_gas_h2 = gasAD(:,1);
S_gas_ch4 = gasAD(:,2);
S_gas_co2 = gasAD(:,3);

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

k = 1;
figure(4)
while k < 4
    if k == 1
    revAD = num2str(k);
    ADstream.revAD = reshape(Conc_AD(:,k),[compADM,length(time)])';
    subplot(2,2,k);
    plot(time,ADstream.revAD);
    title(['Stream ' revAD])
    ylabel('Concentration, kg/m3')
    xlabel('Time, days')
    elseif k == 2
    subplot(2,2,k);
    plot(time,gasAD);
    title(['Stream ' 2])
    ylabel('Concentration, kg/m3')
    xlabel('Time, days')
    else
    revAD = num2str(k);
    ADstream.revAD = reshape(Conc_AD(:,k),[compADM,length(time)])';
    subplot(2,2,k);
    plot(time,ADstream.revAD);
    title(['Stream ' revAD])
    ylabel('Concentration, kg/m3')
    xlabel('Time, days')
    end
    k = k + 1;
end
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
Vol1 = Var.param(22); % primary clarifier m3
Vol2 = Var.param(23); % anoxic m3
Vol3 = Var.param(24); % aeration m3
Vol4 = Var.param(25); % secondary clarifier m3
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
    CODc = 1/1000; % Convert COD ASM1 variables to ADM1 -> g to kg
    Nc = 1/14/1000; % Convert g N to kmol of N
    % Assuming alkalinity is calcium carbonate for C conversion
    Alkc = (0.001/1000)*(12/(12 + 16*3 + 40.078)); % Convert mol/L of Alkalinity to kmol-C/m3
    CODdemand = dCdt(8,13) + 2.86*dCdt(9,13);
%     if CODdemand < (dCdt(2,13) + dCdt(4,13) + dCdt(5,13) + dCdt(6,13))
%         disp('Warning! Influent characterization may need to be evaluated, not enough COD available')
%     else
%     end
    if CODdemand > dCdt(2,13)
    CODdemand = CODdemand - dCdt(2,13);
    dCdt(2,13) = - CODdemand;
        if dCdt(2,13) < 0
            dCdt(2,13) = 0;
        else
        end
    elseif CODdemand < dCdt(2,13)
        dCdt(2,13) = dCdt(2,13) - CODdemand;
        CODdemand = - dCdt(2,13);
        if CODdemand < 0
            CODdemand = 0;
        else
        end
    elseif CODdemand == dCdt(2,13)
        dCdt(2,13) = dCdt(2,13) - CODdemand;
        CODdemand = - dCdt(2,13);
    end
    if CODdemand > dCdt(4,13)
        CODdemand = CODdemand - dCdt(4,13);
        dCdt(4,13) = - CODdemand;
            if dCdt(4,13) < 0
                dCdt(4,13) = 0;
            else
            end
    elseif CODdemand < dCdt(4,13)
        dCdt(4,13) = dCdt(4,13) - CODdemand;
        CODdemand = - dCdt(4,13);
        if CODdemand < 0
            CODdemand = 0;
        else
        end
    elseif CODdemand == dCdt(4,13)
        dCdt(4,13) = dCdt(4,13) - CODdemand;
        CODdemand = - dCdt(4,13);
    end
    if CODdemand > dCdt(5,13)
    CODdemand = CODdemand - dCdt(5,13);
    dCdt(5,13) = - CODdemand;
        if dCdt(5,13) < 0
            dCdt(5,13) = 0;
        else
        end
    elseif CODdemand < dCdt(5,13)
        dCdt(5,13) = dCdt(5,13) - CODdemand;
        CODdemand = - dCdt(5,13);
        if CODdemand < 0
            CODdemand = 0;
        else
        end
    elseif CODdemand == dCdt(5,13)
        dCdt(5,13) = dCdt(5,13) - CODdemand;
        CODdemand = - dCdt(5,13);
    end
    if CODdemand > dCdt(6,13)
    CODdemand = CODdemand - dCdt(6,13);
    dCdt(6,13) = - CODdemand;
        if dCdt(6,13) < 0
            dCdt(6,13) = 0;
        else
        end
    elseif CODdemand < dCdt(6,13)
        dCdt(6,13) = dCdt(6,13) - CODdemand;
        CODdemand = - dCdt(6,13);
        if CODdemand < 0
            CODdemand = 0;
        else
        end
    elseif CODdemand == dCdt(6,13)
        dCdt(6,13) = dCdt(6,13) - CODdemand;
        CODdemand = - dCdt(6,13);
    end
    %% Soluble Organic Nitrogen
    ReqCODs = (dCdt(11,13)*Nc)/l.N_aa;
    if dCdt(2,13)*CODc > ReqCODs
        S_aa_in = ReqCODs;
        S_su_in = dCdt(2,13)*CODc - ReqCODs;
    else
        S_aa_in = dCdt(2,13)*CODc;
        S_su_in = 0.00000000001;
    end
    
    %% Soluble Inert Organic Material
    CODin = dCdt(1,13) + dCdt(2,13) + dCdt(3,13) + dCdt(4,13) + dCdt(5,13) + dCdt(6,13) + dCdt(7,13);
    CODremain = (CODin - dCdt(2,13))*CODc;
    OrgN = dCdt(11,13) + dCdt(12,13)- dCdt(10,13);
    OrgNremain = OrgN*Nc - S_aa_in*l.N_aa;
    ReqOrgNs = l.N_I*dCdt(1,13)*CODc;
    if OrgNremain > ReqOrgNs
        S_I_in = dCdt(1,13)*CODc;
        S_su_in = S_su_in;
    else
        S_I_in = OrgNremain/l.N_I;
        S_su_in = S_su_in + dCdt(1,13)*CODc - S_I_in;
    end
    CODremain = CODremain - dCdt(1,13)*CODc;
    OrgNremain = OrgNremain - S_I_in*l.N_I;
    
    %% Particulate Inert  COD mapping
    ReqOrgNx = l.f_xI_xc*(dCdt(3,13) + dCdt(7,13))*l.N_I*CODc;
    if OrgNremain > ReqOrgNx
        X_I_in = l.f_xI_xc*(dCdt(3,13) + dCdt(7,13))*CODc;
    else
        X_I_in = OrgNremain/l.N_I;
    end
    CODremain = CODremain - X_I_in;
    OrgNremain = OrgNremain - X_I_in*l.N_I;
    
    %% Partitioning of Remaining COD and TKN
    ReqCODXc = OrgNremain/l.N_xc;
    if CODremain > ReqCODXc
        X_c_in = ReqCODXc;
        X_ch_in = (l.f_ch_xc/(l.f_ch_xc + l.f_li_xc))*(CODremain - X_c_in);
        X_li_in = (l.f_li_xc/(l.f_ch_xc + l.f_li_xc))*(CODremain - X_c_in);
        S_IN_in = dCdt(10,13)*Nc;
    else
        X_c_in = CODremain;
        X_ch_in = 0.0000000001; % Cant set to exactly zero
        X_li_in = 0.0000000001; % Cant set to exactly zero
        S_IN_in = dCdt(10,13)*Nc + OrgNremain - X_c_in*l.N_xc;
    end
    % CODremain and TNKremain should be zero.
    %disp(CODremain)
    %disp(OrgNremain)
    S_IC_in = dCdt(13,13)*Alkc;
    S_cat_in = S_IC_in + 0.035;
    S_an_in = S_IN_in;
    % Remaing influent components unaccounted for in the ASM/ADM conversion
    % Cant set to exactly zero
    S_fa_in = 0.0000000001;
    S_va_in = 0.0000000001;
    S_bu_in = 0.0000000001;
    S_pro_in = 0.0000000001;
    S_ac_in = 0.0000000001;
    S_h2_in = 0.0000000001;
    S_ch4_in = 0.0000000001;
    X_pr_in = 0.0000000001;
    X_su_in = 0.0000000001;
    X_aa_in = 0.0000000001;
    X_fa_in = 0.0000000001;
    X_c4_in = 0.0000000001;
    X_pro_in = 0.0000000001;
    X_ac_in = 0.0000000001;
    X_h2_in = 0.0000000001;
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
    % Non-incoming variables
%     S_vam = dCdt(40,13);
%     S_bum = dCdt(41,13);
%     S_prom = dCdt(42,13);
%     S_acm = dCdt(43,13);
%     S_hco3m = dCdt(44,13);
%     S_nh3_in = dCdt(45,13);
%     S_gas_h2 = dCdt(46,13);
%     S_gas_ch4 = dCdt(47,13);
%     S_gas_co2 = dCdt(48,13);
    
%% AD differential equations
% Data/Equations pulled from Aspects on ADM1 implementation within the BSM2 framework
% CHECK PARAMETERS -> CHANGE VOLUMES (headspace/liquid)
    %Separate the components under ENTIRE while function, this part under
    %while (i > compASM) && (i < (compADM + 1))
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
    % CONVERT UNITS CORRECTLY BACK TO ASM1
    % CODconserved = CODt_anaerobic - Sh2 - Sch4;
    % COD Conversions
    dCdt(1,15) = S_I/CODc;
    dCdt(3,15) = X_I/CODc;
    dCdt(2,15) = (S_su + S_aa + S_fa + S_va + S_bu + S_pro + S_ac)/CODc;
    dCdt(4,15) = (X_c + X_ch + X_pr + X_li + X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2)/CODc;
    % TKN Conversions
    dCdt(10,15) = S_IN/Nc;
    dCdt(11,15) = (S_I*l.N_I + S_aa*l.N_aa)/Nc;
    % Ixe is unknown, make it same as ixb
    ixe = ixb/1000; % gN/gCOD -> kgN/kgCOD
    dCdt(12,15) = (l.N_bac*(X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2) + l.N_I*X_I + l.N_xc*X_c + l.N_aa*X_pr - ixe*X_I)/Nc;
    dCdt(5,15) = 0;
    dCdt(6,15) = 0;
    dCdt(7,15) = 0;
    dCdt(8,15) = 0;
    dCdt(13,15) = S_IC/Alkc;
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