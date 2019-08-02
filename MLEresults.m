%% This file pulls data form IndataADM1_v2.m
% The purpose of this file is to take the ODE results and plot certain
% streams, as well as do a carbon mass balance.


function [ASMstream,ADMstream,Cflow] = MLEresults(Array,ODEToc)
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
% North Train
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
% South Train
ASMstream.TwentyOne = reshape(Concentration(:,21),[compASM,length(time)])';
ASMstream.TwentyTwo = reshape(Concentration(:,22),[compASM,length(time)])';
ASMstream.TwentyThree = reshape(Concentration(:,23),[compASM,length(time)])';
ASMstream.TwentyFour = reshape(Concentration(:,24),[compASM,length(time)])';
ASMstream.TwentyFive = reshape(Concentration(:,25),[compASM,length(time)])';
ASMstream.TwentySix = reshape(Concentration(:,26),[compASM,length(time)])';
ASMstream.TwentySeven = reshape(Concentration(:,27),[compASM,length(time)])';
ASMstream.TwentyEight = reshape(Concentration(:,28),[compASM,length(time)])';
ASMstream.TwentyNine = reshape(Concentration(:,29),[compASM,length(time)])';
ASMstream.Thirty = reshape(Concentration(:,30),[compASM,length(time)])';
ASMstream.ThirtyOne = reshape(Concentration(:,31),[compASM,length(time)])';
ASMstream.ThirtyTwo = reshape(Concentration(:,32),[compASM,length(time)])';

% Other Streams
ASMstream.Thirteen = reshape(Concentration(:,13),[compASM,length(time)])';
ASMstream.Fourteen = reshape(Concentration(:,14),[compASM,length(time)])';
ASMstream.Fifteen = reshape(Concentration(:,15),[compASM,length(time)])';
% Sixteen is missing due to it only being the digestion gas stream
ASMstream.Seventeen = reshape(Concentration(:,17),[compASM,length(time)])'; % ASM1 variables that have been converted from ADM1
ASMstream.Eightteen = reshape(Concentration(:,18),[compASM,length(time)])';
ASMstream.Nineteen = reshape(Concentration(:,19),[compASM,length(time)])';
ASMstream.Twenty = reshape(Concentration(:,20),[compASM,length(time)])';

figure(1) % North train
subplot(4,3,1)
plot(time,ASMstream.One)
title('Plant Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,2)
plot(time,ASMstream.Two)
title('Primary Clarifier Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,3)
plot(time,ASMstream.Three)
title('Primary Clarifier WAS')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,4)
plot(time,ASMstream.Four)
title('Mixing Point Before Anoxic Tank')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,5)
plot(time,ASMstream.Five)
title('Anoxic Tank Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,6)
plot(time,ASMstream.Six)
title('Aeration Tank Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,7)
plot(time,ASMstream.Seven)
title('Aeration Tank Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,8)
plot(time,ASMstream.Eight)
title('Internal Recycle')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,9)
plot(time,ASMstream.Nine)
title('Secondary Clarifier Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,10)
plot(time,ASMstream.Ten)
title('Secondary Clarifier Underflow')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,11)
plot(time,ASMstream.Eleven)
title('WAS Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,12)
plot(time,ASMstream.Twelve)
title('RAS Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')

figure(2) % South train 
subplot(4,3,1)
plot(time,ASMstream.TwentyOne)
title('Plant Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,2)
plot(time,ASMstream.TwentyTwo)
title('Primary Clarifier Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,3)
plot(time,ASMstream.TwentyThree)
title('Primary Clarifier WAS')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,4)
plot(time,ASMstream.TwentyFour)
title('Mixing Point Before Anoxic Tank')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,5)
plot(time,ASMstream.TwentyFive)
title('Anoxic Tank Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,6)
plot(time,ASMstream.TwentySix)
title('Aeration Tank Influent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,7)
plot(time,ASMstream.TwentySeven)
title('Aeration Tank Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,8)
plot(time,ASMstream.TwentyEight)
title('Internal Recycle')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,9)
plot(time,ASMstream.TwentyNine)
title('Secondary Clarifier Effluent')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,10)
plot(time,ASMstream.Thirty)
title('Secondary Clarifier Underflow')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,11)
plot(time,ASMstream.ThirtyOne)
title('WAS Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,12)
plot(time,ASMstream.ThirtyTwo)
title('RAS Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')

%% Anaerobic Digester Results
l = IndataADM1_v2; % Read data file
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

% Gas flow
ADMstream.P_gas = ADMstream.P_gas_h2 + ADMstream.P_gas_ch4 + ...
    ADMstream.P_gas_co2 + l.p_gas_h2o; % Total gas pressure
ADMstream.q_gas = Array.QgasAlt; % Total gas flow

% Calculate gas component gas flow
ADMstream.q_gas_h2 = ADMstream.P_gas_h2./ADMstream.P_gas.*ADMstream.q_gas;
ADMstream.q_gas_ch4 = ADMstream.P_gas_ch4./ADMstream.P_gas.*ADMstream.q_gas;
ADMstream.q_gas_co2 = ADMstream.P_gas_co2./ADMstream.P_gas.*ADMstream.q_gas;
ADMstream.q_gas_h2o = l.p_gas_h2o./ADMstream.P_gas.*ADMstream.q_gas;

% Mole fraction of gas components
ADMstream.PyH2 = ADMstream.P_gas_h2./ADMstream.P_gas;
ADMstream.PyCO2 = ADMstream.P_gas_co2./ADMstream.P_gas;
ADMstream.PyCH4 = ADMstream.P_gas_ch4./ADMstream.P_gas;

figure(3)
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
    % Adjusted to only col 2, was all col before
    Carbon.CSM = [Carbon.CSM;Conc_AD(loop.CSM,2)];
    Carbon.CCO2 = [Carbon.CCO2;Conc_AD(loop.CCO2,2)];
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
% Removed extra CSI array, and changed everything besides CSI/CSM/CCO2 to
% their correct structure name (was all CSI)
Carbon.CSIgrams = Carbon.CSI.*Array.Qarray;
Carbon.CSSgrams = Carbon.CSS.*Array.Qarray;
Carbon.CIPgrams = Carbon.CIP.*Array.Qarray;
Carbon.CSBSgrams = Carbon.CSBS.*Array.Qarray;
Carbon.CHBgrams = Carbon.CHB.*Array.Qarray;
Carbon.CABgrams = Carbon.CAB.*Array.Qarray;
Carbon.CUBgrams = Carbon.CUB.*Array.Qarray;
Carbon.CCCgrams = Carbon.CCC.*Array.Qarray;
Carbon.CSMgram = Carbon.CSM.*ADMstream.q_gas.*1000; % Convert to grams
Carbon.CCO2moleC = Carbon.CCO2.*ADMstream.q_gas.*1000; % Convert to mol

% Determine total carbon mass flow rate [gC/day] in stream m
% Multiply carbon fraction by the corresponding carbon component
% Methanol still needs to be implemented
m = 1;
Cflow = ones(length(time),col); % Preallocate array
while m < (col + 1)
    Cflow(:,m) = Carbon.CSIgrams(:,m).*CF.C_S_I + Carbon.CIPgrams(:,m).*CF.C_I_P + Carbon.CSBSgrams(:,m).*CF.C_SB_S +... 
    Carbon.CSSgrams(:,m).*CF.C_S_S + Carbon.CCCgrams(:,m).*CF.C_CC + Carbon.CHBgrams(:,m).*CF.C_HB + Carbon.CABgrams(:,m).*CF.C_AB +...
    Carbon.CUBgrams(:,m).*CF.C_UB;
    m = m + 1;
end
% Methanol stream, as soluble substrate, will need to pull flow from ODE
MethC = 1185000; % COD content of methanol g/m3
QextC = 165;%2; % Flow rate of methanol m3/day
Cflow(:,21) = MethC*QextC*CF.C_S_S;
% Adjust for gas stream
% Convert from COD to grams and mol to grams
Cflow(:,16) = Carbon.CSMgram.*CF.C_S_M + Carbon.CCO2moleC.*12.01; 
% Convert grams to lbs 
Cflow = 0.00220462.*Cflow;
% IGNORE MASS BALANCE CHECK
    % Dynamic system
% Mass balance check
PerError = [];
PerError(:,1) = 100.*(Cflow(:,1) - (Cflow(:,2) + Cflow(:,3)))./Cflow(:,1);
PerError(:,2) = 100.*(Cflow(:,4) - (Cflow(:,12) + Cflow(:,8)+ Cflow(:,2)))./Cflow(:,4);
PerError(:,3) = 100.*(Cflow(:,4) - Cflow(:,5))./Cflow(:,4);
PerError(:,4) = 100.*(Cflow(:,5) - Cflow(:,6))./Cflow(:,5);
PerError(:,5) = 100.*(Cflow(:,6) - (Cflow(:,8) + Cflow(:,7)))./Cflow(:,6);
PerError(:,6) = 100.*(Cflow(:,7) - (Cflow(:,9) + Cflow(:,10)))./Cflow(:,7);
PerError(:,7) = 100.*(Cflow(:,10) - (Cflow(:,12) + Cflow(:,11)))./Cflow(:,10);
PerError(:,8) = 100.*((Cflow(:,1) + Cflow(:,21)) - (Cflow(:,16) + Cflow(:,17) + Cflow(:,19) + Cflow(:,20) + Cflow(:,14)))./(Cflow(:,1) + Cflow(:,21)) ;
PerError(:,9) = 100.*(Cflow(:,13) - (Cflow(:,11) + Cflow(:,3) + Cflow(23) + Cflow(31)))./Cflow(:,13);
PerError(:,10) = 100.*(Cflow(:,13) - Cflow(:,15) - Cflow(:,14))./Cflow(:,13);
PerError(:,11) = 100.*(Cflow(:,15) - Cflow(:,17) - Cflow(:,16))./Cflow(:,15);
PerError(:,12) = 100.*((Cflow(:,9) + Cflow(:,21) + Cflow(29)) - Cflow(:,18))./((Cflow(:,9) + Cflow(:,21)));
PerError(:,13) = 100.*(Cflow(:,18) - Cflow(:,19) - Cflow(:,20))./Cflow(:,18);
% Plot percent error
figure(4)
plot(time,PerError)
title('Mass Balance Percent Error')
ylabel('Error, %')
xlabel('Time, days')

% Plot carbon flow for each stream separately
figure(5)
subplot(4,4,1)
plot(time,Cflow(:,1))
title('Plant Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,2)
plot(time,Cflow(:,2))
title('Primary Clarifier Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,3)
plot(time,Cflow(:,3))
title('Primary Clarifier WAS')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,4)
plot(time,Cflow(:,4))
title('Anoxic Tank Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,5)
plot(time,Cflow(:,5))
title('Anoxic Tank Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,6)
plot(time,Cflow(:,6))
title('Aeration Tank Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,7)
plot(time,Cflow(:,7))
title('SC Influent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,8)
plot(time,Cflow(:,8))
title('Internal Recycle')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,9)
plot(time,Cflow(:,9))
title('SC Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,10)
plot(time,Cflow(:,10))
title('SC Underflow')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,11)
plot(time,Cflow(:,11))
title('WAS Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,12)
plot(time,Cflow(:,12))
title('RAS Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,13)
plot(time,Cflow(:,13))
title('Mixing of WAS and PS')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,14)
plot(time,Cflow(:,16))
title('Anaerobic Digester Gas Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,15)
plot(time,Cflow(:,17))
title('Anaerobic Digester Sludge Stream')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')
subplot(4,4,16)
plot(time,Cflow(:,19))
title('Plant Effluent')
ylabel('Carbon Flow, lb/day')
xlabel('Time, days')

%% Plant Effluent
figure(6)
subplot(2,2,1)
plot(time,ASMstream.Nineteen(:,10))
title('Plant Effluent Ammonia')
ylabel('Concentration, g/m3')
xlabel('Time, days')
subplot(2,2,2)
plot(time,ASMstream.Nineteen(:,9))
title('Plant Effluent Nitrate/Nitrite')
ylabel('Concentration, g/m3')
xlabel('Time, days')
subplot(2,2,3)
plot(time,ASMstream.Nineteen(:,8))
title('Plant Effluent DO')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(2,2,4)
title('AD Gas Production & Inflow')
yyaxis left
plot(time,ADMstream.q_gas)
ylabel('Flowrate, m3/day')
xlabel('Time, days')
yyaxis right
plot(time,Array.Qarray(:,15));

figure(7) %% Random plots
subplot(4,3,1)
plot(time,Array.SRT_NT)
title('MCRT NT')
ylabel('MCRT, days')
xlabel('Time, days')
subplot(4,3,2)
plot(time,ASMstream.Thirteen)
title('Mixing of WAS and PS')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,3)
plot(time,ASMstream.Fourteen)
title('Thickener Centrate Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,4)
plot(time,ASMstream.Fifteen)
title('Thickener Sludge Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,5)
plot(time,ASMstream.Seventeen)
title('Digester Sludge Stream')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,6)
plot(time,Array.Qarray(:,1))
title('Plant Influent Flowrate')
ylabel('Volumetric Flow, m3/day')
xlabel('Time, days')
subplot(4,3,7)
plot(time,ASMstream.Nine(:,10))
title('NT Secondary Clarifier Ammonia')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,8)
plot(time,ASMstream.Nine(:,9))
title('NT Secondary Clarifier Nitrate/Nitrite')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,9)
plot(time,ASMstream.Nine(:,8))
title('NT Secondary Clarifier DO')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,10)
plot(time,ASMstream.TwentyNine(:,10))
title('ST Secondary Clarifier Ammonia')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,11)
plot(time,ASMstream.TwentyNine(:,9))
title('ST Secondary Clarifier Nitrate/Nitrite')
ylabel('Concentration, mg/L')
xlabel('Time, days')
subplot(4,3,12)
plot(time,ASMstream.TwentyNine(:,8))
title('ST Secondary Clarifier DO')
ylabel('Concentration, mg/L')
xlabel('Time, days')

toc
resultsToc = toc;
totalTime = ODEToc + resultsToc;
fprintf('The total elapsed time is %.2f minutes',(totalTime/60))
end