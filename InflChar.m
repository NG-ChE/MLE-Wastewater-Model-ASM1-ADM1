%% This file pulls data from simuPlantData.xlsx
% This files purpose is to pull all the influent data, with typical
% characteristics (cBOD,TKN,NH3,TSS,ALK), and convert it to ASM1 variables

function [sys_int,Var1] = InflChar()
%% Import data for varying influent flow
[~,sheets] = xlsfinfo('simuPlantData.xlsx');
sumData = [];
    for s = 1:numel(sheets)
        [data,~] = xlsread('simuPlantData.xlsx',s);
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
infChar.BOD = 1.5.*infChar.cBOD + 4.6.*infChar.TKN;
%infChar.BOD = infChar.cBOD;
% Alkalinity not given, will have Thursday from Albert
infChar.ALK = InfluentData(:,7); % mg/L
% % Conversion to ASM1 variables using Biological Wastewater Treatment by
% % Grady, Daigger, Love and Filipe, from Chapter 9.6
infChar.CODt = 2.1.*infChar.BOD; % mgCOD/L, Converts BOD5 to total COD, if not available
infChar.CODbo = 1.71.*infChar.BOD; % mgCOD/L, biodegradable COD
infChar.CODio = infChar.CODt - infChar.CODbo; % mgCOD/L, inert COD
% COD fractions taken from B&V study
% The collodial fraction of slowly biodegradable COD wasn't taken into
% account, although this can be implemented by simply multiplying that
% fraction frcol during the clarifier, or anything TSS related, in the main
% ODE function. frcol was equal to 0.312 from the study

frsi = 0.11; % Soluble inert fraction of total COD
frss = 0.1; % Readily biodegradable fraction of total COD
frxi = 0.09; % Particulate inert fraction of total COD
% It should be noted that frxi is very low when compared to other
% wastewater characterstic studies
frsnh = 0.9; % Ammonium fraction of soluble TKN 
insi = 0.05; % N content of soluble inert material, gN/gCOD
inxi = 0.05; % N content of inert particulate material, gN/gCOD

% ASM.Sio = frsi.*infChar.CODt; % mgCOD/L, soluble inert COD
% ASM.Sso = frss.*infChar.CODt;  % mgCOD/L, readily biodegradable substrate
% ASM.Xio = frxi.*infChar.CODt; % mgCOD/L, particulate inert COD
ASM.Xso = (1 - frsi - frss - frxi).*infChar.CODt; % mgCOD/L, slowly biodegradable substrate
ASM.Xbho = repelem(0.000000001,length(InfluentData))'; % mgCOD/L, heterotrophic active biomass -> cant be zero, but very close to it
ASM.Xbao = repelem(0.000000001,length(InfluentData))'; % mgCOD/L, autrophic active biomass -> cant be zero, but very close to it
ASM.Soo = repelem(0,length(InfluentData))'; % mgO2/L, oxygen concentration
ASM.Xpo = repelem(0,length(InfluentData))'; % mgCOD/L, biomass debris 
ASM.Salko = infChar.ALK./100; % mol/L, alkalinity
ASM.Snho = infChar.NH3; % Initial ammonia
ASM.Snoo = infChar.NO3; % Initial nitrite/nitrate
% stkn = (ASM.Snho)./frsnh; % Filtered TKN
% ASM.Sndo = stkn - ASM.Snho - insi.*ASM.Sio; % Soluble organic Nitrogen
% ASM.Xndo = infChar.TKN - stkn - inxi.*ASM.Xio;
% % If any are less than zero
% ASM.Sndo(ASM.Sndo<0) = 0.000000001;
% ASM.Xndo(ASM.Xndo<0) = 0.000000001;
 
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
sys_int = [MLE_int AD_int]; % Combine initial conditions

% Manipulate data for ODE input
Var1.Qt = InfluentData(:,1); % Time, days
Var1.Ct = Var1.Qt; % Time, days
Var1.QflowNT = InfluentData(:,11); % North Train flow m3/day
Var1.QflowST = InfluentData(:,10); % South train flow m3/day
Var1.C = MLE_influent; % Conc at time, Ct

% Correct for duplicates in data by adding a small increment to each element
k = 1;
    while k < (length(Var1.QflowNT)+1)
        if k == 1
           Var1.Qt(k) = 1; % Set time to exactly 1 to avoid interpolation errors on initial points for Qflow and C in ODE
           Var1.Ct(k) = 1;
        else
        Var1.Qt(k) = Var1.Qt(k) + (rand*rand)*1E-5 + 1; 
        Var1.Ct(k) = Var1.Qt(k);
        end
    Var1.QflowNT(k) = Var1.QflowNT(k) + (rand*rand)*1E-3;
    Var1.QflowST(k) = Var1.QflowST(k) + (rand*rand)*1E-3;
    Var1.C(k,:) = Var1.C(k,:) + (rand*rand)*1E-7;
    k = k + 1;
    end
end