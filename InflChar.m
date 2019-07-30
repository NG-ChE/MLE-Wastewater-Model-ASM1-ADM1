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
sys_int = [MLE_int AD_int]; % Combine initial conditions

% Manipulate data for ODE input
Var1.Qt = InfluentData(:,1); % Time, days
Var1.Ct = Var1.Qt; % Time, days
Var1.Qflow = InfluentData(:,11); % Flow at time, Qt
Var1.C = MLE_influent; % Conc at time, Ct

% Correct for duplicates in data by adding a small increment to each element
k = 1;
    while k < (length(Var1.Qflow)+1)
        if k == 1
           Var1.Qt(k) = 1; % Set time to exactly 1 to avoid interpolation errors on initial points for Qflow and C in ODE
           Var1.Ct(k) = 1;
        else
        Var1.Qt(k) = Var1.Qt(k) + (rand*rand)*1E-5 + 1; 
        Var1.Ct(k) = Var1.Qt(k);
        end
    Var1.Qflow(k) = Var1.Qflow(k) + (rand*rand)*1E-3;
    Var1.C(k,:) = Var1.C(k,:) + (rand*rand)*1E-7;
    k = k + 1;
    end
end