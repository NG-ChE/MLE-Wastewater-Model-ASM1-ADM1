tic
%% TO-DO

% Look into changing intial conditions for each stream to possibly quicken
% ODE solving time
% Dynamic influent slowing down ODE by A LOT -> 80X
% Convert ASM values to ADM
% Create GUI for kinetics/ recycle ratios/ simulation time/ flow variables

%% Issues
% ***** CHECK UNITS *****

clc
%clear all

%% Simple simulation of MLE system with no digester.

%% Solver for liquid flow balance of MLE system
Rir = 2.5; % Specify internal recycle
Rr = 0.95; % Specify Waste recycle (has to be less than 1)
if Rr > 1
    Rr = 1;
elseif Rr < 0
    Rr = 0;
else 
end
fpc = 0.01; % Primary Clarifier flow separation
if fpc > 1
    fpc = 1;
elseif fpc < 0
    fpc = 0;
else
end
fsc = 0.389; % Secondary Clarifier flow separation
if fsc > 1
    fsc = 1;
elseif fsc < 0
    fsc = 0;
else
end

%% biological parameters and volumes
% Check parameter values with other reported values
% Make them a function of temperature
% Correct for units
param = [0.67 0.24 0.08 0.08 ...
0.06 4 10 0.2 ...
0.5 0.8 0.8 1 ...
0.4 0.3 0.05 0.05 ...
3 0.1 0.8 40 ...
10 416.83 5522.15 36339.95 ....
4803.84]';

%% simulation time span (days)
t = 1:29;

%% Intial conditions coming into plant
% Units g/m3 = mg/L
% Intial values
TSS = 82; % mg/L
VSS = 61.5; % mg/L
FSS = TSS - VSS; % fixed suspended solids
BOD5 = 155; % mg/L
CODtotal = 325; % mg/L
TKN = 43.5; % mg-N/L
NH3 = 25; % mg-N/L
NO3 = 0; % mg-N/L
ALK = 200; % mg-N/L
% Conversion to ASM1 Variables
CODt = 2.1*BOD5; % mgCOD/L, Converts BOD5 to total COD, if not available
CODbo = 1.71*BOD5; % mgCOD/L, biodegradable COD
CODio = CODt - CODbo; % mgCOD/L, inert COD
Xio = 0.375*1.5*VSS; % mgCOD/L, particulate inert COD
Sio = CODio - Xio; % mgCOD/L, soluble inert COD
f_readily = 0.43; % fraction of biodegradable COD that is readily biodegradable
Sso = CODbo*f_readily; % mgCOD/L, readily biodegradable substrate
Xso = CODbo - Sso; % mgCOD/L, slowly biodegradable substrate
ONtotal = TKN - NH3; % mg-N/L, total organic nitrogen
Snio = 1.5; % mg-N/L, soluble inert organic nitrogen
in_xd = 0.06; % mass of nitrogen per mass of COD in biomass
Xnio = in_xd*Xio; % mg-N/L, 
Snso_Xnso = ONtotal - Snio - Xnio; % mg-N/L, biodegradable organic nitrogen
Sndo = Snso_Xnso*(Sso/(Sso+Xso)); % mg-N/L, soluble biodegradable nitrogen
Xndo = Snso_Xnso - Sndo; % mg-N/L, particulate biodegradable nitrogen
Xbho = 0.000000001; %0; % mgCOD/L, heterotrophic active biomass -> cant be zero, but very close to it
Xbao = 0.000000001; %0; % mgCOD/L, autrophic active biomass -> cant be zero, but very close to it
Soo = 0; % mgO2/L, oxygen concentration
Xpo = 0; % mgCOD/L, biomass debris 
Salko = ALK/100; % mM/L, alkalinity
Snho = NH3; % Initial ammonia
Snoo = 0; % Initial nitrite/nitrate

%% alt method (poopyLab github)
%nb_TKN = TKN*0.03;
%sol_bio_orgN_ratio = Sso/(Sso+ Xso);
%inf_S_NS = (TKN - Snho - nb_TKN)*(sol_bio_orgN_ratio);
%inf_X_NS = (TKN - Snho - nb_TKN)*(1 - sol_bio_orgN_ratio);

%% Intial conditions for function
% assume the initial conditions into plant are those of the reactors

Sio  = 30; %Soluble inert organic matter
Sso  = 69.5; %Readily biodegradable substrate
Xio  = 51.2; %Particulate inert organic matter
Xso  = 202.32; %Slowly biodegradable substrate
Xbho = 28.17; %Active heterotrophic biomass
Xbao = 25; %Active autotrophic biomass
Xpo  = 0; %Particulate products arising from biomass decay
Soo  = 0; %Oxygen
Snoo = 0; %Nitrate and nitrite nitrogen
Snho = 5; %NH 4+ + NH 3 nitrogen
Sndo = 6.95; %Soluble biodegradable organic nitrogen
Xndo = 5; %Particulate biodegradable organic nitrogen
Salko = 7; %Alkalinity
x_int = [Sio Sso Xio Xso Xbho Xbao Xpo Soo Snoo Snho Sndo Xndo Salko]; % Create vector of initial values for the MLE system
x = x_int(:)*ones(1,12); % Format to an array of [components,streams]

%% Import data for varying influent flow

dat_dry=importdata('datos/Inf_dry_2006.txt','\t',1);
dat_rain=importdata('datos/Inf_rain_2006.txt','\t',1);
dat_strm=importdata('datos/Inf_strm_2006.txt','\t',1);
% Append data
aaa=dat_dry.data;
aaa_=dat_strm.data;
aaa_(:,1)=aaa_(:,1)+aaa(end,1);
Qi=[aaa;aaa_];

%% solve diff equation

VarConc_out = zeros(13,12); % Intialize value for assignin
RowTime_out = [];
ColTime_out = [];
StepTime_out = [];
ismemberTime_out = [];
PlantInflow_out = [];
T_out = 0; % Intialize value for assignin
ODE_sol = ode15s(@(t,x) MLE(t,x,param,Qi,Rir,Rr,fpc,fsc),t,x);

%% Import/Manipulate data
%conc_time = importdata('dCdt_var.txt');% Pull the time values out%time = conc_time(14:14:end,1);% Delete every row with time in it%conc_time(14:14:end,:) = [];%Concentration = conc_time;

% Get the data for each stream from the column of Concentration, then
% manipulate each array to be used for graphing
VarConc_out(1:13,:) = []; % Delete dummy rows
Concentration = VarConc_out; 
T_out(1) = []; % Delete dummy value
time = T_out;
% To comment out block of code -> select the code and type "Ctrl" + "R". 
% To uncomment the selected text, click the "Uncomment" button or type "Ctrl" + "T"

streamOne = reshape(Concentration(:,1),[13,length(time)])';
streamTwo = reshape(Concentration(:,2),[13,length(time)])';
streamThree = reshape(Concentration(:,3),[13,length(time)])';
streamFour = reshape(Concentration(:,4),[13,length(time)])';
streamFive = reshape(Concentration(:,5),[13,length(time)])';
streamSix = reshape(Concentration(:,6),[13,length(time)])';
streamSeven = reshape(Concentration(:,7),[13,length(time)])';
streamEight = reshape(Concentration(:,8),[13,length(time)])';
streamNine = reshape(Concentration(:,9),[13,length(time)])';
streamTen = reshape(Concentration(:,10),[13,length(time)])';
streamEleven = reshape(Concentration(:,11),[13,length(time)])';
streamTwelve = reshape(Concentration(:,12),[13,length(time)])';
subplot(4,4,1)
plot(time,streamOne)
title('Stream One')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,2)
plot(time,streamTwo)
title('Stream Two')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,3)
plot(time,streamThree)
title('Stream Three')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,4)
plot(time,streamFour)
title('Stream Four')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,5)
plot(time,streamFive)
title('Stream Five')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,6)
plot(time,streamSix)
title('Stream Six')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,7)
plot(time,streamSeven)
title('Stream Seven')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,8)
plot(time,streamEight)
title('Stream Eight')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,9)
plot(time,streamNine)
title('Stream Nine')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,10)
plot(time,streamTen)
title('Stream Ten')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,11)
plot(time,streamEleven)
title('Stream Eleven')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,12)
plot(time,streamTwelve)
title('Stream Twelve')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,13)
plot(time,streamNine(:,10))
title('Effluent Ammonia')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,14)
plot(time,streamNine(:,9))
title('Effluent Nitrate/Nitrite')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,4,15)
plot(time,PlantInflow_out)
title('Plant flow')
ylabel('Volumetric Flow, m3/day')
xlabel('Time, days')

toc
function [Conc] = MLE(t,dCdt,param,Qi,Rir,Rr,fpc,fsc)
% Fixes structure of intial values from vector to an array
dCdt = reshape(dCdt,[13,12]);
%% Dynamic flow
format long g
t_val = Qi(:,1);
t_val = round(t_val*1000)/1000;
t_step = t;
t_step = round(t_step*1000)/1000;
%assignin('base','DynamicTime',t_val);
%[row,~] = find(t_val == (t_step-1));
tol = 0.0052;
[row,~] = find(abs(t_val - (t_step-1))<tol);

% ** Rather than choose one randomly, implement interpolation between
% the points -> this will also help with data points that are missing
% due to floating point errors, so you can possibly increase the
% tolerance and use interpolation between all values.**

% If the number of elements in row is greater than one, just use first
% element
if numel(row) > 1
    row = row(1);
else
end
%assignin('base','RowTime',row);
%evalin('base','RowTime_out = [RowTime_out;RowTime];');
%assignin('base','ColTime',col);
%evalin('base','ColTime_out = [ColTime_out;ColTime];');
%assignin('base','StepTime',(t_step)-1);
%evalin('base','StepTime_out = [StepTime_out;StepTime];');
%Same_val = ismember(t_val,(1-t_step),'rows');
%Same_val = ismembertol(t_val,(1-t_step),0.001);
%Display = find(Same_val == 1);
%assignin('base','ismemberTime',Same_val);
%evalin('base','ismemberTime_out = [ismemberTime_out;ismemberTime];');
Qi(:,1) = [];
Q_val = Qi(:,14);
Qi(:,14) = [];
Dyn_conc = Qi;
%Influent_conc = Dyn_conc(row,:);
Qplant = Q_val(row);
%Interp_Dyn_conc = interp1(Dyn_conc,
assignin('base','Dynamic_Conc',Dyn_conc);
%% Flow solver
opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
Qint = [100000 99000 1000 300000 300000 300000 200000 250000 100000 50000 50000 300000];
[Q,~] = fsolve(@(Q) myfunc(Q,Qplant,Rir,Rr,fpc,fsc),Qint,opts);
assignin('base','PlantInflow',Q);
evalin('base','PlantInflow_out = [PlantInflow_out;PlantInflow];');
%flows = {'Qpi','Qpo','Qpu','Qi','Qm','Qio','Qo','Qir','Qe','Qu','Qw','Qr'};
%table(Q','RowNames',flows)
function Qsolve = myfunc(Q,Qplant,Rir,Rr,fpc,fsc)
%% Solves flow balance for given plant flow, recycle/waste fractions, and clarifier splits
Qsolve = [Qplant - Q(1);
Q(1) - Q(3) - Q(2);
Rir*Q(2) - Q(8);
Rr*Q(10) - Q(12);
fpc*Q(1) - Q(3);
fsc*Q(7) - Q(10);
Q(2) + Q(8) + Q(12) - Q(4);
Q(5) - Q(4);
Q(6) - Q(5);
Q(7) + Q(8) - Q(6);
Q(9) + Q(10) - Q(7);
Q(12) + Q(11) - Q(10);
Q(3) + Q(9) + Q(11) - Q(1)
];
end
%% Assigning Parameters
% Stoichiometric
Yh = param(1); % gCOD(biomassformed)/gCOD(substrate removed)
Ya = param(2); % gCOD(biomassformed)/gN(N oxidized)
ixb = param(3); % gN/gCOD
fp = param(4);
ixp = param(5); % gN/gCOD

% Kinetic parameters
muh = param(6); % 1/day
Ks = param(7); % gCOD/m3
Koh = param(8); % gN/m3
Kno = param(9); % 0.5 gO2/m3
ng = param(10); 
mua = param(11); % 1/day
Knh = param(12); % gN/m3
Koa = param(13); % gO2/m3
bh = param(14); % 1/day
ba = param(15); % 1/day
ka = param(16); % m3/gCOD/day
kh = param(17); % 1/day
Kx = param(18); 
nh = param(19);

%O2 transference
KLa = param(20); % Oxygen transfer coefficient
So_sat = param(21); % gO2/m3
Vol1 = param(22); % primary clarifier m3
Vol2 = param(23); % anoxic m3
Vol3 = param(24); % aeration m3
Vol4 = param(25); % secondary clarifier m3
%Vol5 = param(26); % denit filters m3

%% State vector
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
%Q(1) - Q(3) - Q(2); % primary clarifier
%Q(2) + Q(8) + Q(12) - Q(4); % inflow to anox 
%Q(5) - Q(4); % anox reactor
%Q(6) - Q(5); % aeration reactor
%Q(7) + Q(8) - Q(6); % recycle split
%Q(9) + Q(10) - Q(7); % secondary clarifier
%Q(12) + Q(11) - Q(10); % RAS/WAS split

%% ODE for MLE system
% Concentration of component i in stream 1-12
% K(1,i)*theta1(1) + K(2,i)*theta1(2) + K(3,i)*theta1(3) + K(4,i)*theta1(4) + K(5,i)*theta1(5) + K(6,i)*theta1(6) + K(7,i)*theta1(7) + K(8,i)*theta1(8) 
i = 1;
Conc = zeros(13,12);
while i < (length(dCdt)+1)
    % Model is under performing compared to simulation results from B&V
    %% Modeling primary clarifier - Otterpohl and Freund 1992
    %hrt = Vol1/Q(1); % Hydraulic residence time
    %n_COD = 2.7*log(hrt*hrt + 9)/100; % Removal efficiency
    XCOD1 = dCdt(3,1) + dCdt(4,1) + dCdt(5,1) + dCdt(6,1) + dCdt(7,1); % Particulate COD in influent
    %CODin = dCdt(1,1) + dCdt(2,1) + XCODin; % Total COD in influent
    %n_x = (n_COD*CODin)/XCODin;
    %if n_x > 0.95
    %    n_x = 0.95;
    %elseif n_x < 0.05
    %    n_x = 0.05;
    %else
    %    n_x = n_x;
    %end
    %% TSS removal
    n_x = 0.533; % Fraction of TSS left in waste effluent, taken as average from ST and NT from Appendix B GPS-X files from B&V
    % Mass of TSS in influent
    massTSS_inf = XCOD1*Q(1);
    % Determine which components are separated

    dCdt(i,1) = Dyn_conc(row,i);
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
    dCdt(i,2) = (dCdt(i,1)*Q(1) - dCdt(i,3)*Q(3))/Q(2); % Mass balance for flow Into/Out of Primary Clarifier

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
    c_x = 0.0015; % Mass fraction of TSS left in waste effluent, taken as average from ST and NT of B&V App. B file
    XCOD7 = dCdt(3,7) + dCdt(4,7) + dCdt(5,7) + dCdt(6,7) + dCdt(7,7); % Particulate COD in influent
    massTSS_7 = XCOD7*Q(7); % Mass of TSS in influent
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
    
    dCdt(i,7) = (dCdt(i,9)*Q(9) + dCdt(i,10)*Q(10))/Q(7); % Mass balance for flow Into/Out of Primary Clarifier
    XCOD9 = dCdt(3,9) + dCdt(4,9) + dCdt(5,9) + dCdt(6,9) + dCdt(7,9); % Particulate COD in influent
    XCOD10 = dCdt(3,10) + dCdt(4,10) + dCdt(5,10) + dCdt(6,10) + dCdt(7,10); % Particulate COD in influent
    
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
        Conc(8,6) = Conc(8,6) + KLa*(So_sat - dCdt(8,6)); % Effect of aireation on the Oxygen concentration
    else 
        Conc(i,6) = Conc(i,6);
    end
    COD = [XCOD1,XCOD7,XCOD9,XCOD10];
    %disp(COD)
    i = i + 1;
end
vec_len = numel(Conc);
Conc = reshape(Conc,[vec_len,1]);
% Save dCdt to file
%if t == 1
%    save('dCdt_var.txt','dCdt','t','-ascii')
%else 
%    save('dCdt_var.txt','dCdt','t','-ascii','-append')
%end
% Due to Conc only changing in the reactors, dCdt variable needs to be
% saved and sent to the workspace.
fprintf('Current simulation time is %6.2f days\n',t)
assignin('base','VarConc',dCdt);
evalin('base','VarConc_out = [VarConc_out;VarConc];');
assignin('base','T_step',t);
evalin('base','T_out(end+1) = T_step;');
end