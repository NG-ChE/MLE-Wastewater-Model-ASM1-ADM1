%% This file pulls data from InflChar.m
% The Var structure is used for biological parameters,
% recycle ratios, the influent concentration and flow,
% as well as the initial values of the ODE.
% The initial values can be setup as "start-up",
% where the initial values of all the ASM streams are the same,
% or can be setup as "steady-state", where intial values were taken from a 
% previous start-up run. The steady state values can be checked by running a
% start-up simulation and checking/changing the values of each stream with
% the initialvaluesAll.m file.
% Start-up should be chosen as a first choice.


function [Var,x] =  influentParam()
% Check parameter values with other reported values
% Make them a function of temperature
% Correct for units
Var.param = [0.67 0.24 0.08 0.08 ...
0.06 4 10 0.2 ...
0.5 0.8 0.8 1 ...
0.4 0.3 0.05 0.05 ...
3 0.1 0.8 100 ...
10 ...
1782.92895 ... % NT Primary clarifier volume m3 (470,000 gal) taken from B&V
5522.15 ... % NT Anoxic tank volume m3 
36339.95 .... % NT Aeration tank volume m3 
19229.891863 ... % NT Secondary clarifier volume m3 (5,080,000 gal) taken from B&V
916.024*7.13232 ... % Denite filter volume, taken as surface area x depth, values from B&V
6056.65888 ... % Digester total volume in m3 (taken as 1.6 MG)
1688.1371 ... % ST PC, m3
4904.42114688 ... % ST anoxic, m3
30283.29427 ... % ST aeration, m3
19215.3873]'; % ST SC, m3

%% Recycle ratios for MLE system
% North Train
Var.RirNT = 1.85; % Specify internal recycle ratio MLR/PC_influent 
Var.RrNT = 0.70; % Specify return activated sludge recycle ratio RAS/PC_influent 
Var.fscNT = 0.974026; % Secondary Clarifier underflow separation (has to be less than or equal to 1)
% South Train
Var.RirST = 1.97; % Specify internal recycle ratio MLR/PC_influent 
Var.RrST = 0.69; % Specify return activated sludge recycle ratio RAS/PC_influent 
Var.fscST = 0.947867; % Secondary Clarifier underflow separation (has to be less than or equal to 1)

Var.ft = 0.1976; % Flow fraction of concentrated TSS stream with respect to inflow

% "start up"
[sys_int,Var1] = InflChar; % Pull influent characteristics from file
Var = catstruct(Var,Var1); % Combine two structures using castruct:https://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct
x = sys_int(:)*ones(1,20); % Format to an array of [components,streams]
% Uncomment line below, and comment out the two above for "steady state"
% variables at t = 1;
%x = initialValuesAll(sys_int);
end