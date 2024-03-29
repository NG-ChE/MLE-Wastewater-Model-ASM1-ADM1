tic
%% TO-DO

% Implement FOG for AD
    % No chemical characeristics given, not sure what to do or if even
    % implement
% Implement primary bypass
% Implement dryer(after AD)
    % Thickener/dewatering implemented
% Implement proper secondary clarifier model
    % Primary clarifier is done but is not removing enough TSS than
    % whats actually happening at the plant

% Create PID for oxygen transfer in aeration zone
    % Close to being done
% Create GUI for kinetics/ recycle ratios/ simulation time/ flow variables

%% Issues
% High nitrate levels in effluent
    % Alkalinity had to be adjusted (twice that of the original)
    % otherwise the ODE solver would stop as the alkalinity would reach 0
    % Change value in InflChar.m
% Denitrification filter not doing anything
    % High O2 concentration undoubtably the cause
    % O2 controller should resolve this issue


clc

%% Pull influent parameters from file
[Var,x] = influentParam(); 

%% Pull AD parameters from file
l = IndataADM1_v2;

%% Simulation time span (days)
% Decrease time_span for more data points in result section, can
% signicantly increase time if interval is very small
time_span = 0.1;
% Uncomment to run simulation for all days from excel file
%t = 1:time_span:fixData(end,1);

t = 1:time_span:75;

% Sample rate for certain variables in ODE
sp = 1/5;
Var.tsample = t + sp*t;
Var.timespan = t;

%% ODE
% Calls the MLE function file
odefunc = @(t,x) MLE(t,x,Var,l);
opts = odeset('MStateDependence','JPattern');
ODE_sol = ode15s(odefunc,t,x,opts);
toc
ODEToc = toc;

%% Results
[ASMstream,ADMstream,Cflow] = MLEresults(Array,ODEToc);

function [Conc] = MLE(t,dCdt,Var,l)
persistent Array
persistent ft
persistent preverrDO_NT
persistent prevt_NT
persistent preverrDO_ST
persistent prevt_ST
CompASM = 13;

%% Dynamic flow
% Interpolate data set of volumetric flow at specified time
Q_NT = interp1(Var.Qt,Var.QflowNT,t); % North train flow, m3/day
Q_ST = interp1(Var.Qt,Var.QflowST,t); % South train flow, m3/day

%% Average flows for steady state simulation
% Comment out previous variables, and uncomment the two lines below, for this case
%Q_NT = 57622.45; % m3/day
%Q_ST = 35982.00; % m3/day
%% Solve flow balance
if t == 1
ft = Var.ft;
else
end
Q = zeros(1,34);
% Stream identification
% Flow rates are in m3/day

% North train
% Q(1) = north train/primary clarifier influent 
% Q(2) = primary clarifier effluent
% Q(3) = primary clarifier waste sludge
% Q(4) = anoxic influent
% Q(5) = anoxic effluent
% Q(6) = aeration Effluent
% Q(7) = secondary clarifier influent
% Q(8) = internal recycle (MLR)
% Q(9) = secondary clarifier effluent
% Q(10) = secondary clarifier underflow
% Q(11) = WAS
% Q(12) = RAS

% South train
% Q(35) = south train influent
% Q(21) = primary clarifier influent
% Q(22) = primary clarifier effluent
% Q(23) = primary clarifier waste sludge
% Q(24) = anoxic influent
% Q(25) = anoxic effluent
% Q(26) = aeration Effluent
% Q(27) = secondary clarifier influent
% Q(28) = internal recycle (MLR)
% Q(29) = secondary clarifier effluent
% Q(30) = secondary clarifier underflow
% Q(31) = WAS
% Q(32) = RAS

% Other streams
% Q(13) = mixing of WAS and PC waste sludge, sent to thickener
% Q(14) = centrate stream
% Q(15) = thickened sludge
% Q(16) = gas flow -> this is generated -> not in mass balance
% Q(17) = liquid flow, sent to dewatering
% Q(18) = denitrification filter stream
% Q(19) = filtered plant effluent stream
% Q(20) = denit filter waste stream
% Q(33) = dewatering cake stream
% Q(34) = dewatering centrate stream

% Flow balance

% North Train
% Q(1) - Q(3) - Q(2); % primary clarifier
% Q(2) + Q(8) + Q(12) - Q(4); % inflow to anox 
% Q(5) - Q(4); % anox basin
% Q(6) - Q(5); % aeration basin
% Q(7) + Q(8) - Q(6); % recycle split
% Q(9) + Q(10) - Q(7); % secondary clarifier
% Q(12) + Q(11) - Q(10); % RAS/WAS split

% South Train
% Q(35) + Q(14) + Q(34) - Q(21); Mixing of Plant flow/ centrate
% Q(21) - Q(22) - Q(23); % PC
% Q(24) - Q(28) - Q(22) - Q(32); % Inflow to anox
% Q(25) - Q(24); % Anox
% Q(26) - Q(25); % Aer
% Q(27) + Q(28) - Q(26); % Recycle split
% Q(29) + Q(30) - Q(27); % SC
% Q(32) + Q(31) - Q(30); % RAS/WAS split

% Other streams
% Q(3) + Q(11) + Q(23) + Q(31) - Q(13) % Mixing of WAS and PC sludge from
% both process trains
% Q(14) + Q(15) - Q(13); % Thickener balance
% Q(15) - Q(17); % Digester balance
% Q(18) - Q(9) - Q(29) - QExtC % Denite filter
% Q(20) + Q(19) - Q(18); % Denite filter TSS removal

% North Train
Q(1) = Q_NT;
Q(3) = 189.27; % Set waste sludge flow, taken from B&V as a control variable
Q(2) = Q(1) - Q(3);
Q(8) = Var.RirNT*Q(2);
Q(12) = Var.RrNT*Q(2);
Q(10) = Q(12)/Var.fscNT;
Q(4) = Q(2) + Q(8) + Q(12);
Q(5) = Q(4);
Q(6) = Q(5);
Q(7) = Q(6) - Q(8);
Q(9) = Q(7) - Q(10);
Q(11) = Q(10) - Q(12);

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
Vol1 = Var.param(22); % NT primary clarifier m3, taken as (total surface area)*(basin depth), from B&V
Vol2 = Var.param(23); % NT anoxic m3, taken as total volume of NT anoxic tanks, from B&V
Vol3 = Var.param(24); % NT aeration m3, taken as total volume of NT aeration tanks, from B&V
Vol4 = Var.param(25); % NT secondary clarifier m3, taken as (total surface area)*(basin depth), from B&V
Vol5 = Var.param(26); % denit filters m3
Vol6 = Var.param(27); % Total AD volume, taken as two digesters in parallel as per B&V study, with a 5200 gal headspace for gas
Vol7 = Var.param(28); % ST PC, m3
Vol8 = Var.param(29); % ST anoxic, m3
Vol9 = Var.param(30); % ST aeration, m3
Vol10 = Var.param(31); % ST SC, m3


%% Component Identification
% Si  = dCdt(1,i); % Soluble inert organic matter
% Ss  = dCdt(2,i); % Readily biodegradable substrate
% Xi  = dCdt(3,i); % Particulate inert organic matter
% Xs  = dCdt(4,i); % Slowly biodegradable substrate
% Xbh = dCdt(5,i); % Active heterotrophic biomass
% Xba = dCdt(6,i); % Active autotrophic biomass
% Xp  = dCdt(7,i); % Particulate products arising from biomass decay
% So  = dCdt(8,i); % Oxygen
% Sno = dCdt(9,i); % Nitrate and nitrite nitrogen
% Snh = dCdt(10,i); % NH 4+ + NH 3 nitrogen
% Snd = dCdt(11,i); % Soluble biodegradable organic nitrogen
% Xnd = dCdt(12,i); % Particulate biodegradable organic nitrogen
% Salk = dCdt(13,i); % Alkalinity

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
% Can be found in the paper: ACTIVATED SLUDGE MODELS ASM1, ASM2, ASM2d AND ASM3

% NT Anoxic tank dCdt(i,5)
theta1 = [muh*(dCdt(2,5)/(Ks+dCdt(2,5)))*(dCdt(8,5)/(Koh+dCdt(8,5)))*dCdt(5,5);...
    muh*(dCdt(2,5)/(Ks+dCdt(2,5)))*(Koh/(Koh+dCdt(8,5)))*(dCdt(9,5)/(Kno+dCdt(9,5)))*ng*dCdt(5,5);...
    mua*(dCdt(10,5)/(Knh+dCdt(10,5)))*(dCdt(8,5)/(Koa+dCdt(8,5)))*dCdt(6,5);...
    bh*dCdt(5,5);...
    ba*dCdt(6,5);...
    ka*dCdt(11,5)*dCdt(5,5);...
    kh*dCdt(4,5)/dCdt(5,5)/(Kx+dCdt(4,5)/dCdt(5,5))*((dCdt(8,5)/(Koh+dCdt(8,5)))+nh*Koh/(Koh+dCdt(8,5))*dCdt(9,5)/(Kno+dCdt(9,5)))*dCdt(5,5);...
    kh*dCdt(4,5)/dCdt(5,5)/(Kx+dCdt(4,5)/dCdt(5,5))*((dCdt(8,5)/(Koh+dCdt(8,5)))+nh*Koh/(Koh+dCdt(8,5))*dCdt(9,5)/(Kno+dCdt(9,5)))*dCdt(5,5)*dCdt(12,5)/dCdt(4,5)];

% NT Aeration tank dCdt(i,6)
theta2 = [muh*(dCdt(2,6)/(Ks+dCdt(2,6)))*(dCdt(8,6)/(Koh+dCdt(8,6)))*dCdt(5,6);...
    muh*(dCdt(2,6)/(Ks+dCdt(2,6)))*(Koh/(Koh+dCdt(8,6)))*(dCdt(9,6)/(Kno+dCdt(9,6)))*ng*dCdt(5,6);...
    mua*(dCdt(10,6)/(Knh+dCdt(10,6)))*(dCdt(8,6)/(Koa+dCdt(8,6)))*dCdt(6,6);...
    bh*dCdt(5,6);...
    ba*dCdt(6,6);...
    ka*dCdt(11,6)*dCdt(5,6);...
    kh*dCdt(4,6)/dCdt(5,6)/(Kx+dCdt(4,6)/dCdt(5,6))*((dCdt(8,6)/(Koh+dCdt(8,6)))+nh*Koh/(Koh+dCdt(8,6))*dCdt(9,6)/(Kno+dCdt(9,6)))*dCdt(5,6);...
    kh*dCdt(4,6)/dCdt(5,6)/(Kx+dCdt(4,6)/dCdt(5,6))*((dCdt(8,6)/(Koh+dCdt(8,6)))+nh*Koh/(Koh+dCdt(8,6))*dCdt(9,6)/(Kno+dCdt(9,6)))*dCdt(5,6)*dCdt(12,6)/dCdt(4,6)];

% Denit Filter dCdt(i,18)
theta3 = [muh*(dCdt(2,18)/(Ks+dCdt(2,18)))*(dCdt(8,18)/(Koh+dCdt(8,18)))*dCdt(5,18);...
    muh*(dCdt(2,18)/(Ks+dCdt(2,18)))*(Koh/(Koh+dCdt(8,18)))*(dCdt(9,18)/(Kno+dCdt(9,18)))*ng*dCdt(5,18);...
    mua*(dCdt(10,18)/(Knh+dCdt(10,18)))*(dCdt(8,18)/(Koa+dCdt(8,18)))*dCdt(6,18);...
    bh*dCdt(5,18);...
    ba*dCdt(6,18);...
    ka*dCdt(11,18)*dCdt(5,18);...
    kh*dCdt(4,18)/dCdt(5,18)/(Kx+dCdt(4,18)/dCdt(5,18))*((dCdt(8,18)/(Koh+dCdt(8,18)))+nh*Koh/(Koh+dCdt(8,18))*dCdt(9,18)/(Kno+dCdt(9,18)))*dCdt(5,18);...
    kh*dCdt(4,18)/dCdt(5,18)/(Kx+dCdt(4,18)/dCdt(5,18))*((dCdt(8,18)/(Koh+dCdt(8,18)))+nh*Koh/(Koh+dCdt(8,18))*dCdt(9,18)/(Kno+dCdt(9,18)))*dCdt(5,18)*dCdt(12,18)/dCdt(4,18)];

% ST Anoxic tank dCdt(i,25)
theta4 = [muh*(dCdt(2,25)/(Ks+dCdt(2,25)))*(dCdt(8,25)/(Koh+dCdt(8,25)))*dCdt(5,25);...
    muh*(dCdt(2,25)/(Ks+dCdt(2,25)))*(Koh/(Koh+dCdt(8,25)))*(dCdt(9,25)/(Kno+dCdt(9,25)))*ng*dCdt(5,25);...
    mua*(dCdt(10,25)/(Knh+dCdt(10,25)))*(dCdt(8,25)/(Koa+dCdt(8,25)))*dCdt(6,25);...
    bh*dCdt(5,25);...
    ba*dCdt(6,25);...
    ka*dCdt(11,25)*dCdt(5,25);...
    kh*dCdt(4,25)/dCdt(5,25)/(Kx+dCdt(4,25)/dCdt(5,25))*((dCdt(8,25)/(Koh+dCdt(8,25)))+nh*Koh/(Koh+dCdt(8,25))*dCdt(9,25)/(Kno+dCdt(9,25)))*dCdt(5,25);...
    kh*dCdt(4,25)/dCdt(5,25)/(Kx+dCdt(4,25)/dCdt(5,25))*((dCdt(8,25)/(Koh+dCdt(8,25)))+nh*Koh/(Koh+dCdt(8,25))*dCdt(9,25)/(Kno+dCdt(9,25)))*dCdt(5,25)*dCdt(12,25)/dCdt(4,25)];

% ST Aeration tank dCdt(i,26)
theta5 = [muh*(dCdt(2,26)/(Ks+dCdt(2,26)))*(dCdt(8,26)/(Koh+dCdt(8,26)))*dCdt(5,26);...
    muh*(dCdt(2,26)/(Ks+dCdt(2,26)))*(Koh/(Koh+dCdt(8,26)))*(dCdt(9,26)/(Kno+dCdt(9,26)))*ng*dCdt(5,26);...
    mua*(dCdt(10,26)/(Knh+dCdt(10,26)))*(dCdt(8,26)/(Koa+dCdt(8,26)))*dCdt(6,26);...
    bh*dCdt(5,26);...
    ba*dCdt(6,26);...
    ka*dCdt(11,26)*dCdt(5,26);...
    kh*dCdt(4,26)/dCdt(5,26)/(Kx+dCdt(4,26)/dCdt(5,26))*((dCdt(8,26)/(Koh+dCdt(8,26)))+nh*Koh/(Koh+dCdt(8,26))*dCdt(9,26)/(Kno+dCdt(9,26)))*dCdt(5,26);...
    kh*dCdt(4,26)/dCdt(5,26)/(Kx+dCdt(4,26)/dCdt(5,26))*((dCdt(8,26)/(Koh+dCdt(8,26)))+nh*Koh/(Koh+dCdt(8,26))*dCdt(9,26)/(Kno+dCdt(9,26)))*dCdt(5,26)*dCdt(12,26)/dCdt(4,26)];


%% ODE for MLE system
% Concentration of component i in any given stream k
% dCdt(i,k)
i = 1; 
Conc = zeros(length(dCdt),numel(Q)); % Initialize Conc array
while i < (CompASM + 1)
    %% North train
    
    %% Modeling primary clarifier - Otterpohl and Freund 1992
    % Model won't be used for the plant simulation
%     hrt = Vol1/Q(1) % Hydraulic residence time
%     n_COD = 2.7*log(hrt*hrt + 9)/100; % Removal efficiency
%     % Particulate COD in influent
%     XCOD1 = dCdt(3,1) + dCdt(4,1) + dCdt(5,1) + dCdt(6,1) + dCdt(7,1); 
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
    n_x = 0.52506; % Fraction of TSS left in effluent, 
    % taken as average of TSS in lbs removed from NT 
    % from Appendix B GPS-X files from B&V
    
    % Determine which components are separated
    % Comment the following 3 lines out if you want to run steady state
    % Interpolate data set of concentration at specified time
    Dyn_conc = interp1(Var.Ct,Var.C,t); 
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
    % Mass balance for flow into/out of Primary Clarifier
    dCdt(i,2) = (dCdt(i,1)*Q(1) - dCdt(i,3)*Q(3))/Q(2); 

    %% Initialize convergence
    for intRec = 1:50 % Number of iterations for convergence solver
    if intRec == 1
    GC8(i,1) = dCdt(i,8); % For first iteration, set concenration guess to that of the initial value for that stream
    GC12(i,1) = dCdt(i,12);
    else
    end
    
    %% Mixing point before anoxic basin
    % mixing point mass balance
    dCdt(i,4) = (Q(2)*dCdt(i,2) + Q(8)*GC8(i,intRec) + ...
        Q(12)*GC12(i,intRec))/Q(4); 
    
    %% Anoxic and aeration zone reactions
    % Reaction rate for anox and aer
    r_anox = K(1,i)*theta1(1) + K(2,i)*theta1(2) + K(3,i)*theta1(3) + ...
        K(4,i)*theta1(4) + K(5,i)*theta1(5) + K(6,i)*theta1(6) + ...
        K(7,i)*theta1(7) + K(8,i)*theta1(8); % anoxic rxn
    r_aer =  K(1,i)*theta2(1) + K(2,i)*theta2(2) + K(3,i)*theta2(3) + ...
        K(4,i)*theta2(4) + K(5,i)*theta2(5) + K(6,i)*theta2(6) + ...
        K(7,i)*theta2(7) + K(8,i)*theta2(8); % aer rxn
    
    Conc(i,5) = 1/Vol2*(Q(4)*dCdt(i,4) - Q(5)*dCdt(i,5)) + r_anox; % Anoxic balance
    Conc(i,6) = 1/Vol3*(Q(5)*dCdt(i,5) - Q(6)*dCdt(i,6)) + r_aer; % Aeration general balance
    if i == 8
        % DO control, needs to be adjusted
        % The set value doesnt correspond to the final result as there is
        % some error between set point and final value
        set_DO = 1; % g/m3, 1 is set -> 1.164 is result under steady state
        % Error between set point and variable result
        errDO_NT = set_DO - dCdt(8,6);
        mbias = 50; % percent of total Kla
        % Controller gains
        Kc = 1.1; Ki = 2.5;
        if t == 1
          prevt_NT = 0.9999999999;
          preverrDO_NT = 0;
        else
        end
        Xvec_NT = [preverrDO_NT,errDO_NT];
        Yvec_NT = [prevt_NT,t];
        controller_NT = mbias + Kc*errDO_NT+ Ki*trapz(Xvec_NT,Yvec_NT,2); 
        preverrDO_NT = errDO_NT; %[preverrDO_NT,errDO_NT];
        prevt_NT = t; %[prevt_NT,t];
        % Prevent large data storage
        if numel(preverrDO_NT) > 100
            preverrDO_NT = preverrDO_NT(end);
        else
        end
        if numel(prevt_NT) > 100
            prevt_NT = prevt_NT(end);
        else
        end
        max_Kla = 80; % Maximum mass transfer of oxygen, [1/day]
        if controller_NT > 100
            controller_NT = 100;
        elseif controller_NT < 0
            controller_NT = 0;
        else
        end
        kla_set = controller_NT*max_Kla/100; % Set KLa value from controller
        % Effect of aeration on the Oxygen concentration
        Conc(8,6) = Conc(8,6) + kla_set*(So_sat - dCdt(8,6)); 
    else 
        Conc(i,6) = Conc(i,6);
    end 
    
    %% Recycle split, same concentration

    dCdt(i,8) = dCdt(i,6);
    dCdt(i,7) = dCdt(i,6);

    %% Secondary clarifier model
    % Model not worth implementing, values taken from B&V will be used from
    % App. B from GPS-X pdf file
    c_x = 0.0017; % Fraction of TSS left in effluent, 
    % taken as lbs of TSS removed from NT of B&V App. B file

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
    % Mass balance for flow into/out of Primary Clarifier
    dCdt(i,7) = (dCdt(i,9)*Q(9) + dCdt(i,10)*Q(10))/Q(7); 

    %% Waste split
    dCdt(i,11) = dCdt(i,10);
    dCdt(i,12) = dCdt(i,10);
    
    %% Convergence solver
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
    
    errT = abs(err1) + abs(err2);
    tol = 0.01;
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
    if intRec > 50
        disp('Convergence failure')
        return
    end
    
    %% Calculate MCRT north train
    % Particulate COD weight in aeration tank, grams
    TSSas = (dCdt(3,6) + dCdt(4,6) + dCdt(5,6) + dCdt(6,6) + dCdt(7,6))*Vol3;
    % Particulate COD concentration in WAS, g/m3
    TSSsc = (dCdt(3,11) + dCdt(4,11) + dCdt(5,11) + dCdt(6,11) + dCdt(7,11)); 
    SLw = TSSsc*Q(11); % Mass flow rate of WAS g/day
    MCRT_NT = (TSSas)/(SLw);
    
    %% South train
    
    Q(35) = Q_ST; % Flow rate after splitter box for south train
    dCdt(i,35) = dCdt(i,1); % Same concentration as north train
    %% Initialize convergence for outer loop (centrate/dewatering stream being recycled back into south train)
    for outRec = 1:50
    if outRec == 1
    GC14(i,1) = dCdt(i,14); % Guess concentration, for centrate stream
    GC34(i,1) = dCdt(i,34); % Guess concentration, for dewatering stream
    GQ14(outRec) = 2430; % Guess flow rate, m3/day, for centrate stream
    GQ34(outRec) = 10; % Guess flow rate, m3/day, for dewatering stream
    else
    end
    % Mixing of south train influent and centrate/dewatering stream
    Q(21) = Q(35) + GQ14(outRec) + GQ34(outRec);  % ST primary clarifier influent
    
    % Mass balance for mixing of centrate streams and ST flow
    dCdt(i,21) = (dCdt(i,35)*Q(35) + GC14(i,outRec)*GQ14(outRec) + GC34(i,outRec)*GQ34(outRec))/Q(21); 

    %% South train flow balances. Refer to beginning of ODE for identification/clarification on streams/variables
    Q(23) = 189.27; % Set waste sludge flow m3/day, taken from B&V as a control variable
    Q(22) = Q(21) - Q(23);
    Q(28) = Var.RirST*Q(22);
    Q(32) = Var.RrST*Q(22);
    Q(30) = Q(32)/Var.fscST;
    Q(24) = Q(22) + Q(28) + Q(32);
    Q(25) = Q(24);
    Q(26) = Q(25);
    Q(27) = Q(26) - Q(28);
    Q(29) = Q(27) - Q(30);
    Q(31) = Q(30) - Q(32);
    
    %% Primary clarifier
    s_x = 0.540723; % Fraction of TSS left in effluent, taken as average 
    % of TSS in lbs removed from ST from Appendix B GPS-X files from B&V
    if i < 3
        dCdt(i,22) = dCdt(i,21);
    elseif (2 < i) && (i < 8)
        dCdt(i,22) = s_x*dCdt(i,21)*Q(21)/Q(22);
    elseif (7 < i) && (i < 12)
        dCdt(i,22) = dCdt(i,21); 
    elseif (11 < i) && (i < 13)
        dCdt(i,22) = s_x*dCdt(i,21)*Q(21)/Q(22);
    else
        dCdt(i,22) = dCdt(i,21);
    end
    if i < 3
        dCdt(i,23) = dCdt(i,21);
    elseif (2 < i) && (i < 8)
        dCdt(i,23) = (1 - s_x)*dCdt(i,21)*Q(21)/Q(23);
    elseif (7 < i) && (i < 12)
        dCdt(i,23) = dCdt(i,21);
    elseif (11 < i) && (i < 13)
        dCdt(i,23) = (1 - s_x)*dCdt(i,21)*Q(21)/Q(23);
    else
        dCdt(i,23) = dCdt(i,21);
    end
    % Mass balance for flow into/out of Primary Clarifier
    dCdt(i,22) = (dCdt(i,21)*Q(21) - dCdt(i,23)*Q(23))/Q(22); 

    %% Initialize convergence for inner loop (MLE south train)
    for innRec = 1:50
    %% Anox/Aer
    if innRec == 1
    GC28(i,1) = dCdt(i,28); % Concentraction guess for streams 28 and 32
    GC32(i,1) = dCdt(i,32); % For first iteration, use the initial value of those streams
    else
    end
    % mixing point
    dCdt(i,24) = (Q(22)*dCdt(i,22) + Q(28)*GC28(i,innRec) + ...
        Q(32)*GC32(i,innRec))/Q(24); 
    
    rAnox =  K(1,i)*theta4(1) + K(2,i)*theta4(2) + K(3,i)*theta4(3) + ...
        K(4,i)*theta4(4) + K(5,i)*theta4(5) + K(6,i)*theta4(6) + ...
        K(7,i)*theta4(7) + K(8,i)*theta4(8); % Anox rxn
    rAer =   K(1,i)*theta5(1) + K(2,i)*theta5(2) + K(3,i)*theta5(3) + ...
        K(4,i)*theta5(4) + K(5,i)*theta5(5) + K(6,i)*theta5(6) + ...
        K(7,i)*theta5(7) + K(8,i)*theta5(8); % Aer rxn
    
    Conc(i,25) = 1/Vol8*(Q(24)*dCdt(i,24) - Q(25)*dCdt(i,25)) + rAnox; % Anoxic balance
    Conc(i,26) = 1/Vol9*(Q(25)*dCdt(i,25) - Q(26)*dCdt(i,26)) + rAer; % Aeration general balance
    
    if i == 8
        % DO control, needs to be adjusted
        % The set value doesnt correspond to the final result as there is
        % some error between set point and final value
        setDO = 0.8; % g/m3 0.8 is set -> 1.847 is result under steady state
        % Error between set point and variable result
        errDO_ST = setDO - dCdt(8,26);
        mbias = 50; % percent of total Kla
        % Controller gains
        Kc = 1.1; Ki = 2.5;
        if t == 1
          prevt_ST = 0.9999999999;
          preverrDO_ST = 0;
        else
        
        end
        Xvec_ST = [preverrDO_ST,errDO_ST];
        Yvec_ST = [prevt_ST,t];
        controller_ST = mbias + Kc*errDO_ST + Ki*trapz(Xvec_ST,Yvec_ST,2);
        preverrDO_ST = errDO_ST; %[preverrDO_ST,errDO_ST];
        prevt_ST = t; %[prevt_ST,t];
        % Prevent large data storage
        if numel(preverrDO_ST) > 100
            preverrDO_ST = preverrDO_ST(end);
        else
        end
        if numel(prevt_ST) > 100
            prevt_ST = prevt_ST(end);
        else
        end
        maxKla = 75; % Maximum mass transfer of oxygen, [1/day]
        if controller_ST > 100
            controller_ST = 100;
        elseif controller_ST < 0
            controller_ST = 0;
        else
        end
        klaSet = controller_ST*maxKla/100; % Set KLa value from controller
        Conc(8,26) = Conc(8,26) + klaSet*(So_sat - dCdt(8,26)); 
    else 
        Conc(i,26) = Conc(i,26);
    end 

    %% Recycle split, same concentration

    dCdt(i,28) = dCdt(i,26);
    dCdt(i,27) = dCdt(i,26);
 
    % Secondary clarifier
    b_x = 0.0012354; % Fraction of TSS left in effluent, taken as average 
    % of TSS in lbs removed from ST from Appendix B GPS-X files from B&V
    
    if i < 3
        dCdt(i,29) = dCdt(i,27);
    elseif (2 < i) && (i < 8)
        dCdt(i,29) = b_x*dCdt(i,27)*Q(27)/Q(29);
    elseif (7 < i) && (i < 12)
        dCdt(i,29) = dCdt(i,27);
    elseif (11 < i) && (i < 13)
        dCdt(i,29) = b_x*dCdt(i,27)*Q(27)/Q(29);
    else
        dCdt(i,29) = dCdt(i,27);
    end
    if i < 3
        dCdt(i,30) = dCdt(i,27);
    elseif (2 < i) && (i < 8)
        dCdt(i,30) = (1 - b_x)*dCdt(i,27)*Q(27)/Q(30);
    elseif (7 < i) && (i < 12)
        dCdt(i,30) = dCdt(i,27);
    elseif (11 < i) && (i < 13)
        dCdt(i,30) = (1 - b_x)*dCdt(i,27)*Q(27)/Q(30);
    else
        dCdt(i,30) = dCdt(i,27);
    end
    % Mass balance for flow into/out of Primary Clarifier
    dCdt(i,27) = (dCdt(i,29)*Q(29) + dCdt(i,30)*Q(30))/Q(27); 

    %% Waste split
    dCdt(i,31) = dCdt(i,30);
    dCdt(i,32) = dCdt(i,30);

    %% Convergence solver
    if dCdt(i,28) <= 0 
        GC28(i,innRec) = 0;
        err1ST = 0;
    else
        err1ST = 100.*((GC28(i,innRec) - dCdt(i,28))./dCdt(i,28));
    end
    if dCdt(i,32) <= 0
        GC32(i,innRec) = 0;
        err2ST = 0;
    else
        err2ST = 100.*((GC32(i,innRec) - dCdt(i,32))./dCdt(i,32));
    end
    
    errTST = abs(err1ST) + abs(err2ST);
    tol = 0.01;
    if errTST > tol
        if err1ST > 0
            q1 = 0.5;
            GC28(i,innRec + 1) = q1.*GC28(i,innRec) + (1 - q1)*dCdt(i,28);
        else
            q1 = -0.25;
            GC28(i,innRec + 1) = q1.*GC28(i,innRec) + (1 - q1)*dCdt(i,28);
        end
        if err2ST > 0
            q2 = 0.5;
            GC32(i,innRec + 1) = q2.*GC32(i,innRec) + (1 - q2)*dCdt(i,32);
        else
            q2 = -0.25;
            GC32(i,innRec + 1) = q2.*GC32(i,innRec) + (1 - q2)*dCdt(i,32);
        end
    else
        break
    end
    end
    if innRec > 50
        disp('Convergence failure')
        return
    end

    %% Calculate MCRT south train
    % Particulate COD in aeration tank, grams
    TSSas = (dCdt(3,26) + dCdt(4,26) + dCdt(5,26) + dCdt(6,26) + dCdt(7,26))*Vol9;
    % Particulate COD in WAS, g/m3
    TSSsc = (dCdt(3,31) + dCdt(4,31) + dCdt(5,31) + dCdt(6,31) + dCdt(7,31)); 
    SLw = TSSsc*Q(31);
    MCRT_ST = (TSSas)/(SLw);

    
    %% Waste sludge mixing
    Q(13) = Q(3) + Q(11) + Q(23) + Q(31); % Combines both WAS and PC sludge from both trains
    % Mass balance on mixing of streams
    dCdt(i,13) = (dCdt(i,11)*Q(11) + dCdt(i,3)*Q(3) + dCdt(i,23)*Q(23) + dCdt(i,31)*Q(31))/Q(13); 
    
    %% Denit. filter reaction
    % External carbon source
    MethC = 1185000; % COD content of methanol g/m3
    ExtC = [0 MethC 0 0 0 0 0 0 0 0 0 0 0];
    QextC = 165; % Flow rate of methanol m3/day, needs to be adjusted
    Q(18) = Q(9) + QextC + Q(29); % Adding the two process train flow rates plus the external carbon source
    dCdt(i,18) = (Q(9)*dCdt(i,9) + QextC*ExtC(i) + dCdt(i,29)*Q(29))/Q(18); 
    
    % Reaction rate of denitrification filter
    r_den =  K(1,i)*theta3(1) + K(2,i)*theta3(2) + K(3,i)*theta3(3) + ...
        K(4,i)*theta3(4) + K(5,i)*theta3(5) + K(6,i)*theta3(6) + ...
        K(7,i)*theta3(7) + K(8,i)*theta3(8);
    
    Conc(i,18) = 1/Vol5*((Q(9)*dCdt(i,9)+ QextC*ExtC(i) + ...
        dCdt(i,29)*Q(29)) - Q(18)*dCdt(i,18)) + r_den;
    
    %% Denit. filter TSS removal
    % Control variable unknown,waste stream set to 0.05 percent of inflow
    Q(20) = 0.0005*Q(18); % Waste stream
    Q(19) = Q(18) - Q(20); % Plant effluent
    DFx = 1 - 0.9; % 90 percent of solids removed, taken from B&V
    if i < 3
        dCdt(i,19) = dCdt(i,18);
    elseif (2 < i) && (i < 8)
        dCdt(i,19) = DFx*dCdt(i,18)*Q(18)/Q(19);
    elseif (7 < i) && (i < 12)
        dCdt(i,19) = dCdt(i,18);
    elseif (11 < i) && (i < 13)
        dCdt(i,19) = DFx*dCdt(i,18)*Q(18)/Q(19);
    else
        dCdt(i,19) = dCdt(i,18);
    end
    if i < 3
        dCdt(i,20) = dCdt(i,18);
    elseif (2 < i) && (i < 8)
        dCdt(i,20) = (1 - DFx)*dCdt(i,18)*Q(18)/Q(20);
    elseif (7 < i) && (i < 12)
        dCdt(i,20) = dCdt(i,18);
    elseif (11 < i) && (i < 13)
        dCdt(i,20) = (1 - DFx)*dCdt(i,18)*Q(18)/Q(20);
    else
        dCdt(i,20) = dCdt(i,18);
    end
    % Mass balance across denit. filter removal
    dCdt(i,18) = (dCdt(i,19)*Q(19) + dCdt(i,20)*Q(20))/Q(18); 

    %% Thickener
    % TSS concentration is control variable
    TSS15 = dCdt(3,15) + dCdt(4,15) + dCdt(5,15) + dCdt(6,15) + dCdt(7,15);
    % The set TSS is not exact for some reason, but is nonetheless fairly constant
    % The set TSS will need to be adjusted manually for the actual target
    % TSS desired
    maxTSS = 46000/1.55;
    minTSS = 44000/1.55;
    setTSS = 45000/1.55;
    if TSS15 > maxTSS || TSS15 < minTSS
        ft = (dCdt(3,13) + dCdt(4,13) + dCdt(5,13) + dCdt(6,13) ...
            + dCdt(7,13))/(setTSS);
            if ft > 1
                ft = 1;
            elseif ft < 0
                ft = 0.001;
            else
            end
    else
    end
%     If you want to check TSS concentration vs flow rate during
%     simulation, pause and insert the uncommented lines below in the command
%     window
%     yyaxis left
%     plot(Array.tArray,Array.Qarray(:,15));
%     yyaxis right
%     plot(Array.tArray,Array.TSS15)
    Q(15) = ft*Q(13); % Digester influent
    Q(14) = Q(13) - Q(15); % Centrate stream
    Q(17) = Q(15); % Digester effluent
    Q(16) = 0; % Digester gas stream not in flow balance
   
    % Removal efficiency of 0.8 of TSS, taken from B&V
    TSSx = 1 - 0.8;
    if i < 3
        dCdt(i,14) = dCdt(i,13);
    elseif (2 < i) && (i < 8)
        dCdt(i,14) = TSSx*dCdt(i,13)*Q(13)/Q(14);
    elseif (7 < i) && (i < 12)
        dCdt(i,14) = dCdt(i,13); 
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
    % Mass balance of thickener
    dCdt(i,13) = (dCdt(i,14)*Q(14) + dCdt(i,15)*Q(15))/Q(13);
    
    %% Conversion from ASM1 to ADM1
    % Using paper: Benchmark Simulation Model No.2 (BSM2), source: http://iwa-mia.org/wp-content/uploads/2019/04/BSM_TG_Tech_Report_no_3_BSM2_General_Description.pdf
    
    % Initialize adm and xtemp array
    % xtemp consistents of asm variables
    lenComp = 1:13;
    lenADM = 1:24;
    xtemp = zeros(1,numel(lenComp));
    adm = zeros(1,numel(lenADM));
    for lenComp = 1:13
        xtemp(lenComp) = dCdt(lenComp,15);
    end
    for lenADM = 1:24
        adm(lenADM) = 0;
    end
    
    %% Parameters
    tfac = 1/298.15 - 1/(273.15 + 35); % Operating temperature 35C
    bigr = 0.08314; % R gas constant [bar M^-1 K^-1]
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
    
    % If extreme detail was used then some extra NH4 would be transformed
    % into N bound in biomass and some biomass would be formed when
    % removing the CODdemand (based on the yield). But on a total COD balance
    % approach the below is correct (neglecting the N need for biomass growth)
    %% Reducing total incoming COD for Ss,Xs,Xbh,Xba in that specific order
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
    
    % SS becomes part of amino acids when transformed into ADM
    % and any remaining SS is mapped to monosacharides (no N contents)
    % Enough SND must be available for mapping to amino acids
    % Saa COD equivalent to SND
    
    sorgn = dCdt(11,15)/fnaa;
    if sorgn >= xtemp(2)
        % Not all SND-N fits into amino acids
        % Map all SS COD into Saa 
        adm(2) = xtemp(2);
        % The remaining part of SND 
        xtemp(11) = xtemp(11) - xtemp(2)*fnaa;
        % All SS used
        xtemp(2) = 0;
    else
        % All SND-N fits into amino acids
        % Map all SND related COD into Saa
        adm(2) = sorgn;
        % The rest of the SS COD, later mapped into sugar
        xtemp(2) = xtemp(2) - sorgn;
        % All SND used
        xtemp(11) = 0;
    end
    
    % XS becomes part of Xpr (proteins) when transformed into ADM and any
    % remaining XS is mapped to Xch and Xli (no N contents). Enough XND
    % must be available for mapping to XPr
    % Xpr COD equivalent to XND
    xorgn = dCdt(12,15)/fnaa;
    if xorgn >= xtemp(4)
        % Not all XND-N fits into Xpr
        % Map all XS COD into Xpr
        xprtemp = xtemp(4);
        % The remaining part of XND
        xtemp(12) = xtemp(12) - xtemp(4)*fnaa;
        % All XS used
        xtemp(4) = 0;
        xlitemp = 0;
        xchtemp = 0;
    else
        % All XND-N fits into Xpr
        % Map all XND related COD into Xpr
        xprtemp = xorgn;
        % Part of XS COD not associated with N
        xlitemp = frlixs*(xtemp(4) - xorgn);
        % Part of XS COD not associated with N
        xchtemp = (1 - frlixs)*(xtemp(4) - xorgn);
        % All XS used
        xtemp(4) = 0;
        % All XND used
        xtemp(12) = 0;
    end
    
    % Biomass becomes part of Xpr and XI when transformed into ADM
    % and any remaining XBH and XBA is mapped to Xch and Xli (no N contents)
    % Remaining XND-N can be used as nitrogen source to form Xpr
    biomass = xtemp(5) + xtemp(6);
    % Part which is mapped to XI
    biomass_nobio = biomass*(1 - frxs);
    biomass_bioN = biomass*fnbac - biomass_nobio*fxni;
    xchtemp2 = 0;
    xlitemp2 = 0;
    xprtemp2 = 0;
    if biomass_bioN < 0
        disp('ERROR: not enough biomass N to map the requested inert part')
    end
    if ((biomass_bioN/fnaa) <= (biomass - biomass_nobio))
        % All biomass N used
        xprtemp2 = biomass_bioN/fnaa;
        remainCOD = biomass - biomass_nobio - xprtemp2;
        % Use part of remaining XND-N to form proteins
        if ((xtemp(12)/fnaa) > remainCOD)
            xprtemp2 = xprtemp2 + remainCOD;
            xtemp(12) = xtemp(12) - remainCOD*fnaa;
            remainCOD = 0;
            xtemp(5) = 0;
            xtemp(6) = 0;
        % Use all remaining XND-N to form proteins    
        else
            xprtemp2 = xprtemp2 + xtemp(12)/fnaa;
            remainCOD = remainCOD - xtemp(12)/fnaa;
            xtemp(12) = 0;
        end
        % Part of the COD not associated with N 
        xlitemp2 = frlixb*remainCOD;
        % Part of the COD not associated with N 
        xchtemp2 = (1 - frlixb)*remainCOD;
    else
        % All biomass COD used
        xprtemp2 = biomass - biomass_nobio;
        % Any remaining N in XND
        xtemp(12) = xtemp(12) + biomass*fnbac - biomass_nobio*fxni - xprtemp2*fnaa;
    end
    
    xtemp(5) = 0;
    xtemp(6) = 0;
    
    % Direct mapping of XI and XP to ADM1 XI
    % Assumption: same N content in both ASM1 and ADM1 particulate inerts
    inertX = (1 - fdegrade)*(xtemp(3) + xtemp(7));
    
    % Special case: if part of XI and XP in ASM can be degraded in AD
    % we have no knowledge about the contents so we put it in as composits (xc)
    % we need to keep track of the associated nitrogen
    % N content may be different, take first from XI&XP-N, then XND-N, then
    % SND-N, then SNH. Similar principle could be used for other states.
    xc = 0;
    xlitemp3 = 0;
    xchtemp3 = 0;
    if fdegrade > 0
        noninertX = fdegrade*(xtemp(3) + xtemp(7));
        % N in XI & XP (ASM) not enough
        if ((noninertX*fxni) < (noninertX*fnxc))
            xc = noninertX*fxni/fnxc;
            noninertX = noninertX - noninertX*fxni/fnxc;
            % N in XND not enough
            if xtemp(12) < noninertX*fnxc
                xc = xc + xtemp(12)/fnxc;
                noninertX = noninertX - xtemp(12)/fnxc;
                xtemp(12) = 0;
                % N in SND not enough
                if (xtemp(11) < noninertX*fnxc)
                    xc = xc + temp(11)/fnxc;
                    noninertX = nonintertX - xtemp(11)/fnxc;
                    xtemp(11) = 0;
                    % N in SNH not enough
                    if (xtemp(10) < noninertX*fnxc)
                        xc = xc + xtemp(10)/fnxc;
                        noninertX = noninertX - xtemp(10)/fnxc;
                        xtemp(10) = 0;
                        disp('ERROR: Nitrogen shortage when converting biodegradable XI & XP')
                        disp('Putting remaining XI & XP as lipids (50%) and carbohydrates (50%)') 
                        xlitemp3 = 0.5*noninertX;
                        xchtemp3 = 0.5*noninertX;
                        noninertX = 0;
                    % N in SNH enough for mapping
                    else
                        xc = xc + noninertX;
                        xtemp(10) = xtemp(10) - noninertX*fnxc;
                        noninertX = 0;
                    end
                % N in SNH enough for mapping
                else
                    xc = xc + noninertX;
                    xtemp(11) = xtemp(11) - noninertX*fnxc;
                    noninertX = 0;
                end
            % N in SND enough for mapping    
            else
                xc = xc + noninertX;
                xtemp(12) = xtemp(12) - noninertX*fnxc;
                noninertX = 0;
            end
        % N in XI & XP (ASM) enough for mapping
        else
            xc = xc + noninertX;
            % Put remaining N as XND
            xtemp(12) = xtemp(12) + noninertX*(fxni - fnxc);
            noninertX = 0;
        end
    end
    
    % Mapping of ASM SI to ADM1 SI
    % N content may be different, take first from SI-N, then SND-N, then XND-N,
    % then SNH. Similar principle could be used for other states.
    inertS = 0;
    % N in SI (ASM) not enough
    if ((xtemp(1)*fsni) < (xtemp(1)*fsni_adm))
        inertS = xtemp(1)*fsni/fsni_adm;
        xtemp(1) = xtemp(1) - xtemp(1)*fsni/fsni_adm;
        % N in SND not enough
        if (xtemp(11) < (xtemp(1)*fsni_adm))
            inertS = inertS + xtemp(11)/fsni_adm;
            xtemp(1) = xtemp(1) - xtemp(11)/fsni_adm;
            xtemp(11) = 0;
            % N in XND not enough
            if (xtemp(12) < (xtemp(1)*fsni_adm))
                inertS = inertS + xtemp(12)/fsni_adm;
                xtemp(1) = xtemp(1) - xtemp(12)/fsni_adm;
                xtemp(12) = 0;
                % N in SNH not enough
                if (xtemp(10) < (xtemp(1)*fsni_adm))
                    inertS = inertS + xtemp(10)/fsni_adm;
                    xtemp(1) = xtemp(1) - xtemp(10)/fsni_adm;
                    xtemp(10) = 0;
                    disp('ERROR: Nitrogen shortage when converting SI')
                    disp('Putting remaining SI as monosacharides') 
                    xtemp(2) = xtemp(2) + xtemp(1);
                    xtemp(1) = 0;
                % N in SNH enough for mapping
                else
                    inertS = inertS + xtemp(1);
                    xtemp(10) = xtemp(10) - xtemp(1)*fsni_adm;
                    xtemp(1) = 0;
                end
            % N in XND enough for mapping
            else
                inertS = inertS + xtemp(1);
                xtemp(12) = xtemp(12) - xtemp(1)*fsni_adm;
                xtemp(1) = 0;
            end
        % N in SND enough for mapping
        else
            inertS = inertS + xtemp(1);
            xtemp(11) = xtemp(11) - xtemp(1)*fsni_adm;
            xtemp(1) = 0;
        end
    % N in SI (ASM) enough for mapping
    else
        inertS = inertS + xtemp(1);
        % Put remaining N as SND
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

%% AD differential equations
% Full list of equations can be found in the paper: Aspects on ADM1 Implementation within the BSM2 Framework
% Paper pdf source: https://www.iea.lth.se/publications/Reports/LTH-IEA-7224.pdf
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

    q_gasAlt = ((l.R*l.T_op)/(l.P_atm - l.p_gas_h2o))*l.V_liq*((rho_T_8/16) + (rho_T_9/64) + rho_T_10);
    if q_gasAlt < 0 
        q_gasAlt = 0; % Running initial conditions as "start up", 
                      % under influentParam.m, can cause weird side
                      % effects.
                      % This accounts for the odd case of a single instance
                      % of a negative value.
    else
    end
    
%% Conversion from ADM1 to ASM1
    % Initialize asmm array and xtemp array
    % xtemp is full of the adm components to be used for conversion to asm
    vec3 = 1:13;
    vec4 = 1:24;
    asmm = zeros(1,numel(vec3)); 
    xtemp = zeros(1,numel(vec4));
    for vec3 = 1:13
        asmm(vec3) = 0;
    end
    for vec4 = 1:24
        xtemp(vec4) = dCdt(13 + vec4,17);
    end
% Set parameter values
    % Refer to ASM to ADM conversion above for specifics about each
    % parameter. Some parameters have been multipled by 1000 to account for
    % units. kg -> g
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
    
    % Biomass becomes part of XS and XP when transformed into ASM
    % Assume Npart of formed XS to be fnxc and Npart of XP to be fxni
    % Remaining N goes into ammonia pool
    biomass = 0;
    for vec5 = 17:23
        biomass = biomass + xtemp(vec5)*1000;
    end
    % Part of biomass mapping into XP
    biomass_nobio = biomass*(1 - frxs_AS);
    biomass_bioN = biomass*fnbac - biomass_nobio*fxni;
    remainCOD = 0;
    if biomass_bioN < 0
        disp('Warning 1')
        xptemp = biomass*fnbac/fxni;
        biomass_nobio = xptemp; % added from matlab code of BSM2
        biomass_bioN = 0;
    else
        xptemp = biomass_nobio;
    end
    if ((biomass_bioN/fnxc) <= (biomass - biomass_nobio))
        % All biomass used
        xstemp = biomass_bioN/fnxc;
        remainCOD = biomass - biomass_nobio - xstemp;
        % Use part of remaining S_IN to form XS
        if ((xtemp(11)*14000/fnxc) >= remainCOD)% added >= instead of > from matlab code of BSM2 
            xstemp = xstemp + remainCOD;
        else
            disp('Warning 2')
            disp('System failure')
        end
    else
        % All biomass COD used
        xstemp = biomass - biomass_nobio;
    end
    % Any remaining N in S_IN
    xtemp(11) = xtemp(11) + biomass*fnbac/14000 - xptemp*fxni/14000 ...
        - xstemp*fnxc/14000;
    
    % XS = all X from ADM except Xi and biomass
    asmm(4) = xstemp;
    for vec6 = 13:16
        asmm(4) = asmm(4) + xtemp(vec6)*1000;
    end
    
    % Inert part of biomass
    asmm(7) = xptemp;
    
    % Mapping of intert XI in AD into XI and possibly XS in AS
    % Assumption: same N content in both ASM1 and ADM1 particulate inerts
    % Special case: if part of XI in AD can be degraded in AS
    % We have no knowledge about the contents so we put it in as part substrate
    % (XS)
    % We need to keep track of the associated nitrogen
    % N content might be different, take first from XI-N then S_IN
    inertX = (1 - fdegrade_AS)*xtemp(24)*1000;
    xstemp2 = 0;
    noninertX = 0;
    if (fdegrade_AS > 0)
        noninertX = fdegrade_AS*xtemp(24)*1000;
        % N in XI (AD) not enough
        if fxni < fnxc
            xstemp2 = noninertX*fxni/fnxc;
            noninertX = noninertX - noninertX*fnxi/fnxc;
            if ((xtemp(11)*14000) < (noninertX*fnxc))
                % N in SNH not enough
                xstemp2 = xstemp2 + xtemp(11)*14000/fnxc;
                noninertX = noninertX - xtemp(11)*14000/fnxc;
                xtemp(11) = 0;
                disp('Warning 3')
                inertX = inertX + noninertX;
                % N in S_In enough for mapping
            else
                xstemp2 = xstemp2 + noninertX;
                xtemp(11) = xtemp(11) - noninertX*fnxc/14000;
                noninertX = 0;
            end
        % N in XI (AD) enough for mapping
        else
            xstemp2 = xstemp2 + noninertX;
            % Put remaining N as S_IN
            xtemp(11) = xtemp(11) + noninertX*(fnxi - fnxc)/14000;
            noninertX = 0;
        end
    end
    
    % Xi = Xi*fdegrade_AS + possibly nonmappable XS
    asmm(3) = inertX;
    % Extra added XS (biodegradable XI)
    asmm(4) = asmm(4) + xstemp2;
    
    % Mapping of ADM SI into ASM1 SI
    % N content may be different, take first from SI_N then S_IN
    inertS = 0;
    % N in SI (AD) not enough
    if (fsni_adm < fsni)
        inertS = xtemp(12)*fsni_adm/fsni;
        xtemp(12) = xtemp(12) - xtemp(12)*fsni_adm/fsni; % added from matlab code of BSM2
        if ((xtemp(11)*14) < (xtemp(12)*fsni))
            % N in S_IN not enough
            inertS = inertS + xtemp(11)*14/fsni;
            xtemp(12) = xtemp(12) - xtemp(11)*14/fsni;
            xtemp(11) = 0;
            disp('Warning 5')
        % N in S_IN enough for mapping
        else
            inertS = inertS + xtemp(12);
            xtemp(11) = xtemp(11) - xtemp(12)*fsni/14;
            xtemp(12) = 0;
        end
    % N in SI (AD) enough for mapping
    else
        inertS = inertS + xtemp(12);
        % Put remaining N as S_IN
        xtemp(11) = xtemp(11) + xtemp(12)*(fsni_adm - fsni)/14;
        xtemp(12) = 0;
    end
    
    % Si as SI
    asmm(1) = inertS*1000;
    
    % Nitrogen in biomass, composites, proteins
    % Xnd is the nitrogen part of XS in ASM1. Should be based on
    % the same variables as constitutes XS, xc and xpo (no nitrogen in lipids
    % and carbohydrates)
    asmm(12) = fnxc*(xstemp + xstemp2) + nxc*xtemp(13) + naa*xtemp(15);
    
    % Sh2 = x(8) and Sch4 = x(9) assumed to be stripped upon reentry to ASM side
    % Ss = sum of all solubles except Sh2, Sch4, Si, Sic, Sin
    for vec7 = 1:7
        asmm(2) = asmm(2) + xtemp(vec7)*1000;
    end
    
    % Snd is the nitrogen part of Ss in ASM1. Therefore Snd should be based on
    % the same variables as constitutes Ss, and we assume
    % there is only nitrogen in the amino acids. The N content of Si is
    % not included in ASM1 
    asmm(11) = naa*xtemp(2);
    
    % Snh = Sin including adjustments above
    asmm(10) = xtemp(11)*14000;
    
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

    % CODconserved = CODt_anaerobic - Sh2 - Sch4;
    % Set the asmm values to their respective components in the dCdt vector
    % for stream 17
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
    
    %% Dewatering unit
    d_x = 1 - 0.85; % TSS left in centrate, taken from B&V
    % Control variable is cake lb/day -> set to 26,000 lb/day or 11,793,402
    % grams
    set_cakeWT = 11793402; % grams
    concT = dCdt(1,33) + dCdt(2,33) + dCdt(3,33) + dCdt(4,33) + ...
      dCdt(5,33) + dCdt(6,33) + dCdt(7,33) + dCdt(8,33) + ...
      dCdt(9,33) + dCdt(10,33) + dCdt(11,33) + ...
      dCdt(12,33) + dCdt(13,33); % Sum of all concentrations
    CakeWT = concT*Q(33); % Weight of cake in g/day
    cakeFr = (set_cakeWT/concT)/Q(17); % New flow split fraction for specified cake weight
    if cakeFr > 1
      %disp('Error: cake weight requirement too high')
      cakeFr = 0.99;
    elseif cakeFr < 0
      %disp('Error: cake weight requirement too low')
      cakeFr = 0.01;
    else
    end
    
    Q(33) = cakeFr*Q(17); % Cake flow
    Q(34) = Q(17) - Q(33); % Centrate balance
    
    % Separate the TSS components
    if i < 3
        dCdt(i,34) = dCdt(i,17);
    elseif (2 < i) && (i < 8)
        dCdt(i,34) = d_x*dCdt(i,17)*Q(17)/Q(34);
    elseif (7 < i) && (i < 12)
        dCdt(i,34) = dCdt(i,17); 
    elseif (11 < i) && (i < 13)
        dCdt(i,34) = d_x*dCdt(i,17)*Q(17)/Q(34);
    else
        dCdt(i,34) = dCdt(i,17);
    end
    if i < 3
        dCdt(i,33) = dCdt(i,17);
    elseif (2 < i) && (i < 8)
        dCdt(i,33) = (1 - d_x)*dCdt(i,13)*Q(17)/Q(33);
    elseif (7 < i) && (i < 12)
        dCdt(i,33) = dCdt(i,17);
    elseif (11 < i) && (i < 17)
        dCdt(i,33) = (1 - d_x)*dCdt(i,17)*Q(17)/Q(33);
    else
        dCdt(i,33) = dCdt(i,13);
    end
    
    % Mass balance on dewatering stream
    dCdt(i,17) = (dCdt(i,34)*Q(34) + dCdt(i,33)*Q(33))/Q(17);
    
    %% Convergence solver for outer loop
    if dCdt(i,14) <= 0 
        GC14(i,outRec) = 0;
        err14 = 0;
    else
        err14 = 100.*((GC14(i,outRec) - dCdt(i,14))./dCdt(i,14));
    end
    if dCdt(i,34) <= 0
        GC34(i,outRec) = 0;
        err34 = 0;
    else
        err34 = 100.*((GC34(i,outRec) - dCdt(i,34))./dCdt(i,34));
    end
    if Q(14) <= 0 
        GQ14(outRec) = 0;
        errQ14 = 0;
    else
        errQ14 = 100.*((GQ14(outRec) - Q(14))./Q(14));
    end
    if Q(34) <= 0
        GQ34(outRec) = 0;
        errQ34 = 0;
    else
        errQ34 = 100.*((GQ34(outRec) - Q(34))./Q(34));
    end
    
    errTout = abs(err14) + abs(err34) + abs(errQ14) +abs(errQ34);
    tol = 0.01;
    if errTout > tol
        % Concentration
        if err14 > 0
            q1 = 0.5;
            GC14(i,outRec + 1) = q1.*GC14(i,outRec) + (1 - q1)*dCdt(i,14);
        else
            q1 = -0.25;
            GC14(i,outRec + 1) = q1.*GC14(i,outRec) + (1 - q1)*dCdt(i,14);
        end
        if err34 > 0
            q2 = 0.5;
            GC34(i,outRec + 1) = q2.*GC34(i,outRec) + (1 - q2)*dCdt(i,34);
        else
            q2 = -0.25;
            GC34(i,outRec + 1) = q2.*GC34(i,outRec) + (1 - q2)*dCdt(i,34);
        end
        % Flow
        if errQ14 > 0
            q1 = 0.5;
            GQ14(outRec + 1) = q1.*GQ14(outRec) + (1 - q1)*Q(14);
        else
            q1 = -0.25;
            GQ14(outRec + 1) = q1.*GQ14(outRec) + (1 - q1)*Q(14);
        end
        if errQ34 > 0
            q2 = 0.5;
            GQ34(outRec + 1) = q2.*GQ34(outRec) + (1 - q2)*Q(34);
        else
            q2 = -0.25;
            GQ34(outRec + 1) = q2.*GQ34(outRec) + (1 - q2)*Q(34);
        end
    else
        break
    end
    end
    if outRec > 50
        disp('Convergence failure')
        return
    end
    
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
    Array.QgasAlt = q_gasAlt;
    Array.SRT_NT = MCRT_NT;
    Array.SRT_ST = MCRT_ST;
    Array.CakeWT = CakeWT;
elseif t >= 1.01
    % Avoid large dataset and quicken solver time by only aquiring data at the specified time span
    AdjT = round(t*100)/100;
    FindT = find(AdjT == Var.timespan);
    if numel(FindT) > 0
        Array.mleArray = [Array.mleArray;dCdt];
        Array.tArray = [Array.tArray;t];
        Array.Qarray = [Array.Qarray;Q];
        Array.TSS15 = [Array.TSS15;TSS15];
        Array.ft = [Array.ft;ft];
        Array.pH = [Array.pH;pH];
        Array.QgasAlt = [Array.QgasAlt;q_gasAlt];
        Array.SRT_NT = [Array.SRT_NT;MCRT_NT];
        Array.SRT_ST = [Array.SRT_ST;MCRT_ST];
        Array.CakeWT = [Array.CakeWT;CakeWT];
        assignin('base','Array',Array);
    else
    end
else
end
fprintf('Current simulation time is %6.6f days\n',t)
end