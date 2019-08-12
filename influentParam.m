%% This file pulls data from InflChar.m and initialValuesAll.m
% The Var structure is used for biological parameters,
% recycle ratios, the influent concentration and flow,
% as well as the initial values of the ODE.
% The initial values can be setup as "start-up",
% where the initial values of all the ASM streams are the same,
% or can be setup as "steady-state", where intial values were taken from a 
% previous start-up run. The steady state values can be checked by running a
% start-up simulation and checking/changing the values of each stream with
% the initialvaluesAll.m file.
% Start-up should be chosen as a first choice, then adjusted.


function [Var,x] =  influentParam()
% Parameters correspond to those of ASM1, at 20 degree C
% Values taken from GPS-X, (B&V)
% Except for ng and nh -> Taken from book: ACTIVATED SLUDGE MODELS
% ASM1, ASM2, ASM2d AND ASM3
Var.param = [0.666 ... % Yh [gCOD/gCOD]
    0.24 ... % Ya [gCOD/gCOD]
    0.086 ... % ixb [gN/gCOD]
    0.08 ... % fp [gCOD/gCOD]
    0.06 ... % ixp [gN/gCOD]
    6 ... % muh [1/day]
    20 ... % Ks [gCOD/m3]
    0.2 ... % Koh [gO2/m3]
    0.5 ... % Kno [gN/m3]
    0.8 ... % ng [-]
    0.8 ... % mua [1/day]
    1 ... % Knh [gN/m3]
    0.4 ... % Koa [gO2/m3]
    0.62 ... % bh [1/day]
    0.04 ... % ba [1/day]
    0.08 ... % ka [m3/gCOD/day]
    3 ... % kh [1/day]
    0.03 ... % Kx [gCOD/gCOD]
    0.4 ... % nh [-]
    100 ... % KLa [1/day]
    9.42 ... % O2 saturation concentration mg/L
    1782.92895 ... % NT Primary clarifier volume m3 
    5522.15 ... % NT Anoxic tank volume m3 
    36339.95 .... % NT Aeration tank volume m3 
    19229.891863 ... % NT Secondary clarifier volume m3
    916.024*7.13232 ... % Denite filter volume, taken as surface area x depth, values from B&V
    6056.65888 ... % Digester total volume in m3 (taken as 1.6 MG)
    1688.1371 ... % ST PC, m3
    4904.42114688 ... % ST anoxic, m3
    30283.29427 ... % ST aeration, m3
    19215.3873]'; % ST SC, m3

%% Recycle ratios for MLE system
% It should be noted that although the WAS is a control variable, with a
% value of 160 gpm, it can't be specified if the RAS recycle ratio is
% already specified, which is the case. Even though this is the case, the
% WAS should be very close to that value.
% North Train
Var.RirNT = 1.85; % Specify internal recycle ratio MLR/PC_influent, calculated from flow data average [Flow_June2018_2019.xlsx]
Var.RrNT = 0.70; % Specify return activated sludge recycle ratio RAS/PC_influent , calculated from flow data average [Flow_June2018_2019.xlsx]
Var.fscNT = 0.974026; % Secondary Clarifier underflow separation (has to be less than or equal to 1) (taken from B&V results)
% South Train
Var.RirST = 1.97; % Specify internal recycle ratio MLR/PC_influent, calculated from flow data average  [Flow_June2018_2019.xlsx]
Var.RrST = 0.69; % Specify return activated sludge recycle ratio RAS/PC_influent, calculated from flow data average  [Flow_June2018_2019.xlsx]
Var.fscST = 0.947867; % Secondary Clarifier underflow separation (has to be less than or equal to 1) (taken from B&V results)

Var.ft = 0.1976; % Flow fraction of concentrated TSS stream with respect to inflow. It is just a set fraction, the value is dynamically changed inside ODE function

% "start up", initial conditions for each stream are the same, except for
% digester streams
[sys_int,Var1] = InflChar; % Pull influent characteristics from file
Var = catstruct(Var,Var1); % Combine two structures using castruct: https://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct
x = sys_int(:)*ones(1,34); % Format to an array of [components,streams]

% "steady state"
% Uncomment line below, and comment out the two above for steady state
% values as initial conditions for each stream
%x = initialValuesAll(sys_int);
end