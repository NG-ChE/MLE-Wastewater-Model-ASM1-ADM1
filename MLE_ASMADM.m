tic
%% TO-DO
% Convert ASM values to ADM
    % UNITS NEED TO BE ADJUSTED
    % Essentially done
    % Need to verify and find missing variables
    
% Implement ODE/system equations for ADM1
    % d_1_code,IndataADM1_v2,ADM1_fun_v2_ODE
% Implement thickener(before AD)/dewatering(after AD)/dryer(after AD)
% Implement denitrification filter

% Create PID for oxygen transfer in aeration zone
% Create GUI for kinetics/ recycle ratios/ simulation time/ flow variables

%% Issues
% ***** CHECK UNITS *****

clc
%clear all

%% Simple simulation of MLE system with no digester.

%% Recycle ratios for MLE system
Rir = 1.88; % Specify internal recycle ratio MLR/PC_effluent 
Rr = 0.71; % Specify return activated sludge recycle ratio RAS/PC_effluent 
fpc = 0.017554; % Primary Clarifier flow separation (has to be less than or equal to 1)
if fpc > 1
    fpc = 1;
elseif fpc < 0
    fpc = 0;
else
end
fsc = 0.974; % Secondary Clarifier underflow separation (has to be less than or equal to 1)
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
t = 1:0.1:29;
% Sample rate
    % Decrease sample rate more better DO control, but ODE takes longer
sp = 1/5;
tsample = t + sp*t;
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
x_int = [Sio Sso Xio Xso Xbho Xbao Xpo Soo Snoo Snho Sndo Xndo Salko]; % Create vector of initial values for the MLE system
x = x_int(:)*ones(1,12); % Format to an array of [components,streams] -> maybe send to ODE function to allow the initial conditions to be adjusted for added streams

%% Import data for varying influent flow
dat_dry=importdata('datos/Inf_dry_2006.txt','\t',1);
dat_rain=importdata('datos/Inf_rain_2006.txt','\t',1);
dat_strm=importdata('datos/Inf_strm_2006.txt','\t',1);
% Append data
drydata = dat_dry.data;
stormdata = dat_strm.data;
stormdata(:,1) = stormdata(:,1) + drydata(end,1);
InfluentData = [drydata;stormdata];
DataTime = InfluentData(:,1); % Time
% Manipulate data for ODE input
Qt = InfluentData(:,1); % Time, days
Ct = Qt; % Time, days
InfluentData(:,1) = []; % Remove column of time data
Qflow = InfluentData(:,14); % Flow at time, Qt
InfluentData(:,14) = []; % Remove column of influent flow data
C = InfluentData; % Conc at time, Ct
% Correct for duplicates in data by adding a small increment to each element
k = 1;
while k < (length(Qflow)+1)
    if k == 1
        Qt(k) = 1; % Set time to exactly 1 to avoid interpolation errors on initial points for Qflow and C in ODE
        Ct(k) = 1;
    else
    Qt(k) = Qt(k) + (rand*rand)*1E-5 + 1; 
    Ct(k) = Qt(k);
    end
Qflow(k) = Qflow(k) + (rand*rand)*1E-3;
C(k,:) = C(k,:) + (rand*rand)*1E-7;
k = k + 1;
end
  
%% ODE
% Intialize variables for assignin
VarConc_out = zeros(13,12); 
PlantInflow_out = [];
T_out = [];
MCRT_out = [];
% Optimize ODE solver
figure(1)
opts = odeset('MStateDependence','JPattern','Stats','on','OutputFcn',@odeplot);
ODE_sol = ode15s(@(t,x) MLE(t,x,param,Rir,Rr,fpc,fsc,Qflow,Qt,Ct,C,tsample),t,x,opts);

%% Results
% Get the data for each stream from the column of Concentration, then
% manipulate each array to be used for graphing
VarConc_out(1:13,:) = []; % Delete dummy rows
Concentration = VarConc_out; 
time = T_out;
j = 1;
stream = [];
[row,col] = size(Concentration);
figure(2)
% Plot results
while j < (16+1)
    if j < (col+1)
    rev = num2str(j);
    stream.rev = reshape(Concentration(:,j),[13,length(time)])';
    subplot(4,4,j);
    plot(time,stream.rev);
    title(['Stream ' rev])
    ylabel('Concentration, mg/L')
    xlabel('Time, days')
    else
    end
    if j == 6
        subplot(4,4,16)
        plot(time,stream.rev(:,8))
        title('DO')
        ylabel('Concentration, mg/L')
        xlabel('Time, days')
    elseif j == 9
        subplot(4,4,13)
        plot(time,stream.rev(:,10))
        title('Effluent Ammonia')
        ylabel('Concentration, mg/L')
        xlabel('Time, days')
        subplot(4,4,14)
        plot(time,stream.rev(:,9))
        title('Effluent Nitrate/Nitrite')
        ylabel('Concentration, mg/L')
        xlabel('Time, days')
    elseif j == 15
        subplot(4,4,j)
        plot(time,PlantInflow_out)
        title('Plant flow')
        ylabel('Volumetric Flow, m3/day')
        xlabel('Time, days')

    else
    end
    j = j + 1;
end

toc
function [Conc] = MLE(t,dCdt,param,Rir,Rr,fpc,fsc,Qflow,Qt,Ct,C,tsample)
% Fixes structure of intial values from vector to an array
dCdt = reshape(dCdt,[13,12]);
%% Dynamic flow
Qplant = interp1(Qt,Qflow,t); % Interpolate data set of volumetric flow at specified time
%% Flow solver
opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
Qint = [20000 19000 100 60000 60000 60000 25000 17000 10000 10000 500 10000];
[Q,~] = fsolve(@(Q) myfunc(Q,Qplant,Rir,Rr,fpc,fsc),Qint,opts);
assignin('base','PlantInflow',Q);
evalin('base','PlantInflow_out = [PlantInflow_out;PlantInflow];');
function Qsolve = myfunc(Q,Qplant,Rir,Rr,fpc,fsc)
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
Qsolve = [Qplant - Q(1);
Q(1) - Q(3) - Q(2);
Rir*Q(2) - Q(8);
Rr*Q(2) - Q(12);
fpc*Q(1) - Q(3);
fsc*Q(10) - Q(12);
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

% O2 transference
KLa = param(20); % Oxygen transfer coefficient
So_sat = param(21); % gO2/m3

% Equipment volumes
Vol1 = param(22); % primary clarifier m3
Vol2 = param(23); % anoxic m3
Vol3 = param(24); % aeration m3
Vol4 = param(25); % secondary clarifier m3
%Vol5 = param(26); % denit filters m3

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
while i < (length(dCdt)+1)
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
    Dyn_conc = interp1(Ct,C,t); % Interpolate data set of concentration at specified time
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
    
    dCdt(i,7) = (dCdt(i,9)*Q(9) + dCdt(i,10)*Q(10))/Q(7); % Mass balance for flow Into/Out of Primary Clarifier
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
%         % Rough control of DO in aeration zone
%         tadj = round(t*100)/100;
%         t_find = find(tadj == tsample);
%         dCdt(10,6);
%         dCdt(8,6);
%         if numel(t_find) > 0
%             if dCdt(10,6) > 5
%                 KLa = param(20);
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
    % Not executable yet
    % Add new streams into Q and initial conditions
    % COD demand
%     dCdt(i,13) = (dCdt(i,12)*Q(12) + dCdt(i,3)*Q(3))/Q(13);
%     CODdemand = dCdt(8,13) + 2.86*dCdt(9,13);
%     % Reducing total incoming COD for Ss,Xs,Xbh,Xba in that specific order
%     if CODdemand < (dCdt(2,13) + dCdt(4,13) + dCdt(5,13) + dCdt(6,13))
%         disp('Warning! Influent characterization may need to be evaluated, not enough COD available')
%     else
%     end
%     if CODdemand > dCdt(2,13)
%     CODdemand = CODdemand - dCdt(2,13);
%     dCdt(2,13) = -CODdemand;
%         if dCdt(2,13) < 0
%             dCdt(2,13) = 0;
%         else
%         end
%     elseif CODdemand < dCdt(2,13)
%         dCdt(2,13) = dCdt(2,13) - CODdemand;
%         CODdemand = - dCdt(2,13);
%         if CODdemand < 0
%             CODdemand = 0;
%         else
%         end
%     elseif CODdemand == dCdt(2,13)
%         dCdt(2,13) = dCdt(2,13) - CODdemand;
%         CODdemand = - dCdt(2,13);
%     end
%     if CODdemand > dCdt(4,13)
%         CODdemand = CODdemand - dCdt(4,13);
%         dCdt(4,13) = -CODdemand;
%             if dCdt(4,13) < 0
%                 dCdt(4,13) = 0;
%             else
%             end
%     elseif CODdemand < dCdt(4,13)
%         dCdt(4,13) = dCdt(4,13) - CODdemand;
%         CODdemand = - dCdt(4,13);
%         if CODdemand < 0
%             CODdemand = 0;
%         else
%         end
%     elseif CODdemand == dCdt(4,13)
%         dCdt(4,13) = dCdt(4,13) - CODdemand;
%         CODdemand = - dCdt(4,13);
%     end
%     if CODdemand > dCdt(5,13)
%     CODdemand = CODdemand - dCdt(5,13);
%     dCdt(5,13) = -CODdemand;
%         if dCdt(5,13) < 0
%             dCdt(5,13) = 0;
%         else
%         end
%     elseif CODdemand < dCdt(5,13)
%         dCdt(5,13) = dCdt(5,13) - CODdemand;
%         CODdemand = - dCdt(5,13);
%         if CODdemand < 0
%             CODdemand = 0;
%         else
%         end
%     elseif CODdemand == dCdt(5,13)
%         dCdt(5,13) = dCdt(5,13) - CODdemand;
%         CODdemand = - dCdt(5,13);
%     end
%     if CODdemand > dCdt(6,13)
%     CODdemand = CODdemand - dCdt(6,13);
%     dCdt(6,13) = -CODdemand;
%         if dCdt(6,13) < 0
%             dCdt(6,13) = 0;
%         else
%         end
%     elseif CODdemand < dCdt(6,13)
%         dCdt(6,13) = dCdt(6,13) - CODdemand;
%         CODdemand = - dCdt(6,13);
%         if CODdemand < 0
%             CODdemand = 0;
%         else
%         end
%     elseif CODdemand == dCdt(6,13)
%         dCdt(6,13) = dCdt(6,13) - CODdemand;
%         CODdemand = - dCdt(6,13);
%     end
%     % Soluble Organic Nitrogen
%     ReqCODs = dCdt(11,13)/Naa; % Naa is nitrogen fraction of the amino acid state variable, Saa
%     if dCdt(2,13) > ReqCODs
%         Saa = ReqCODs;
%         Ssu_A = dCdt(2,13) - ReqCODs;
%     else
%         Saa = dCdt(2,13);
%         Ssu_A = 0;
%     end
%     % Soluble Inert Organic Material
%     CODin = dCdt(1,13) + dCdt(2,13) + dCdt(3,13) + dCdt(4,13) + dCdt(5,13) + dCdt(6,13) + dCdt(7,13);
%     CODremain = CODin - dCdt(2,13);
%     OrgN = dCdt(11,13) + dCdt(12,13)- dCdt(10,13);
%     OrgNremain = OrgN - Saa*Naa;
%     ReqOrgNs = Ni*dCdt(1,13);
%     if OrgNremain > ReqOrgNs
%         Si_ADM = dCdt(1,13);
%         Ssu = Ssu_A;
%     else
%         Si_ADM = OrgNremain/Ni;
%         Ssu = Ssu_A + dCdt(1,13) - Si_ADM;
%     end
%     CODremain = CODremain - dCdt(1,13);
%     OrgNremain = OrgNremain - Si_ADM*Ni;
%     % Particulate Inert  COD mapping
%     ReqOrgNx = fxi*(dCdt(3,13) + dCdt(7,13))*Ni;
%     if OrgNremain > ReqOrgNx
%         Xi_ADM = fxi*(dCdt(3,13) + dCdt(7,13));
%     else
%         Xi_ADM = OrgNremain/Ni;
%     end
%     CODremain = CODremain - Xi_ADM;
%     OrgNremain = OrgNremain - Xi_ADM*Ni;
%     % Partitioning of Remaining COD and TKN
%     ReqCODXc = OrgNremain/Nxc;
%     if CODremain > ReqCODXc
%         Xc = ReqCODXc;
%         Xch = (fch_xc/(fch_xc + fli_xc))*(CODremain - Xc);
%         Xli = (fli_xc/(fch_xc + fli_xc))*(CODremain - Xc);
%         Sin = dCdt(10,13);
%     else
%         Xc = CODremain;
%         Xch = 0;
%         Xli = 0;
%         Sin = dCdt(10,13) + OrgNremain - Xc*Nxc;
%     end
%     % CODremain and TNKremain should be zero.
%     Scat = dCdt(13,13) + 0.035;
%     San = Sin;
%          
%     %% Conversion from ADM1 to ASM1
%     CODconserved = CODt_anaerobic - Sh2 - Sch4;
%     % COD Conversions
%     dCdt(1,effluent) = Si_ADM;
%     dCdt(3,eff) = Xi_ADM;
%     dCdt(2,eff) = Ssu + Saa + Sfa + Sva + Sbu + Spro + Sac;
%     dCdt(4,eff) = Xc + Xch + Xpr + Xli + Xsu + Xaa + Xfa + Xc4 + Xpro + Xac + Xh2;
%     % TKN Conversions
%     dCdt(10,eff) = Sin;
%     dCdt(11,eff) = Si_ADM*Ni + Saa*Naa;
%     dCdt(12,eff) = Nbac*(Xsu + Xaa + Xfa + Xc4 + Xpro + Xac + Xh2) + Ni*Xi + Nxc*Xc + Naa*Xpr - ixe*Xi;
%     dCdt(5,eff) = 0;
%     dCdt(6,eff) = 0;
%     dCdt(7,eff) = 0;
%     dCdt(8,eff) = 0;
%     dCdt(13,eff) = Sic; 
    i = i + 1;
end
vec_len = numel(Conc);
Conc = reshape(Conc,[vec_len,1]);
% Due to Conc only changing in the reactors, dCdt variable needs to be
% saved and sent to the workspace.
fprintf('Current simulation time is %6.2f days\n',t)
assignin('base','VarConc',dCdt);
evalin('base','VarConc_out = [VarConc_out;VarConc];');
assignin('base','T_step',t);
evalin('base','T_out(end+1) = T_step;');
end