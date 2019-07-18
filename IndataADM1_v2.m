%Master Thesis
% 
%Applied Physics
%Chalmers University of Technolopgy
%Spring 2014
% 
%Oskar Danielsson
%oskard@student.chalmers.se
% 
%File for writing indata to ODE implementation of ADM1. The file saves all
%valuable information in the struct *l*, which later can be used through
%the rest of the implementation.
% 
%Names of variables and parameters are used from Rosen and Jeppsson
%implementation in report "Aspect of ADM1 Implementaion within the BSM2
%Framework", 2006.
% 
%Values of variables do allso comes from the report by Rosen and Jeposson,
%unless stated otherwise.
function [l]=IndataADM1_v2
%==========================================================================
%Stoichiometric parameter values
%==========================================================================
%Parameter Value Unit
%...........................................
l.f_sI_xc = 0.1; % ?
l.f_xI_xc = 0.2; % ?
l.f_ch_xc = 0.2; % ?
l.f_pr_xc = 0.2; % ?
l.f_li_xc = 0.3; % ?
l.N_xc = 0.0376/14; %kmole N (kg COD)^?1
l.N_I = 0.06/14; %kmole N (kg COD)^?1
l.N_aa = 0.007; %kmole N (kg COD)^?1
l.C_xc = 0.02786; %kmole C (kg COD)^?1
l.C_sI = 0.03; %kmole C (kg COD)^?1
l.C_ch = 0.0313; %kmole C (kg COD)^?1
l.C_pr = 0.03; %kmole C (kg COD)^?1
l.C_li = 0.022; %kmole C (kg COD)^?1
l.C_xI = 0.03; %kmole C (kg COD)^?1
%...........................................
l.C_su = 0.0313; %kmole C (kg COD)^?1
%...........................................
l.C_aa = 0.03; %kmole C (kg COD)^?1
%...........................................
l.f_fa_li = 0.95; % ?
l.C_fa = 0.0217; %kmole C (kg COD)^?1D.3. ADM1 indata file 86
%...........................................
l.f_h2_su = 0.19; % ?
l.f_bu_su = 0.13; % ?
l.f_pro_su = 0.27; % ?
l.f_ac_su = 0.41; % ?
l.N_bac = 0.08/14; %kmole N (kg COD)^?1
l.C_bu = 0.025; %kmole C (kg COD)^?1
l.C_pro = 0.0268; %kmole C (kg COD)^?1
l.C_ac = 0.0313; %kmole C (kg COD)^?1
l.C_bac = 0.0313; %kmole C (kg COD)^?1
l.Y_su = 0.1; % ?
%...........................................
l.f_h2_aa = 0.06; % ?
l.f_va_aa = 0.23; % ?
l.f_bu_aa = 0.26; % ?
l.f_pro_aa = 0.05; % ?
l.f_ac_aa = 0.40; % ?
l.C_va = 0.024; %kmole C (kg COD)^?1
l.Y_aa = 0.08; % ?
%...........................................
l.Y_fa = 0.06; % ?
%...........................................
l.Y_c4 = 0.06; % ?
%...........................................
l.Y_pro = 0.04; % ?
%...........................................
l.C_ch4 = 0.0156; %kmole C (kg COD)^?1
l.Y_ac = 0.05; % ?
%...........................................
l.Y_h2 = 0.06; % ?
%...........................................
%Extra input to ADM1 model (added more substartes)
%...........................................
l.C_lac = 0.0313; %kmole C (kg COD)^?1
l.Y_lac_f = 0.055; % ?
l.Y_lac_o = 0.055; % ?
%...........................................
%==========================================================================
%Physicochemical parameter values
%==========================================================================
%Parameter Value Unit
%..........................................................................
l.R = 0.083145; %bar M^?1 K^?1
l.T_base = 298.15; %K
l.T_op = 308.15; %K
%..........................................................................
l.K_w = 1e-14*exp(55900/(l.R*100)*(1/l.T_base - 1/l.T_op)); %M
l.K_a_va = 10^-4.86; %M
l.K_a_bu = 10^-4.82; %M
l.K_a_pro = 10^-4.88; %M
l.K_a_ac = 10^-4.76; %M
l.K_a_co2 = 10^-6.35*exp(7646/(l.R*100)*(1/l.T_base - 1/l.T_op)); %M
l.K_a_IN = 10^-9.25*exp(51965/(l.R*100)*(1/l.T_base - 1/l.T_op));%M
%..........................................................................
l.k_A_B_va = 1e10; %M^?1 d^?1
l.k_A_B_bu = 1e10; %M^?1 d^?1
l.k_A_B_pro = 1e10; %M^?1 d^?1
l.k_A_B_ac = 1e10; %M^?1 d^?1
l.k_A_B_co2 = 1e10; %M^?1 d^?1
l.k_A_B_IN = 1e10; %M^?1 d^?1
%..........................................................................
l.P_atm = 1.013; %bar
l.p_gas_h2o = 0.0313*exp(5290*(1/l.T_base - 1/l.T_op)); %bar
l.k_p = 5e4; %m^3 d^?1 bar^?1
%..........................................................................
l.k_L_a = 200; %d^?1
%..........................................................................
l.K_H_co2 = 0.035*exp(-19410/(l.R*100)*(1/l.T_base - 1/l.T_op)); %M_liq bar^?1
l.K_H_ch4 = 0.0014*exp(-14240/(l.R*100)*(1/l.T_base - 1/l.T_op)); %M_liq bar^?1
l.K_H_h2 = 7.8e-4*exp(-4180/(l.R*100)*(1/l.T_base - 1/l.T_op)); %M_liq bar^?1
%..........................................................................
%Extra input to ADM1 model (added more substartes)
%..........................................................................
l.K_a_lac = 10^-3.86;
%..........................................................................
%==========================================================================
%Physical parameter values
%==========================================================================
%Parameter Value Unit
%...........................................
l.V_liq = 3400; %m^3
l.V_gas = 300; %m^3
%...........................................
%==========================================================================
%Biochemical parameter values Comment
%==========================================================================
%Parameter Value Unit
%...........................................
l.k_dis = 0.5; %d^?1
l.k_hyd_ch = 10; %d^?1
l.k_hyd_pr = 10; %d^?1
l.k_hyd_li = 10; %d^?1
%...........................................
l.K_S_IN = 1e-4; %M
%...........................................
l.k_m_su = 30; %d^?1
l.K_S_su = 0.5; %kg COD m^?3
l.pH_UL_aa = 5.5; % ?
l.pH_LL_aa = 4; % ?
%...........................................
l.k_m_aa = 50; %d^?1
l.K_S_aa = 0.3; %kg COD m^?3
%...........................................
l.k_m_fa = 6; %d^?1
l.K_S_fa = 0.4; %kg COD m^?3
l.K_I_h2_fa = 5e-6; %kg COD m^?3
%...........................................
l.k_m_c4 = 20; %d^?1
l.K_S_c4 = 0.2; %kg COD m^?3
l.K_I_h2_c4 = 1e-5; %kg COD m^?3
%...........................................
l.k_m_pro = 13; %d^?1
l.K_S_pro = 0.1; %kg COD m^?3
l.K_I_h2_pro = 3.5e-6; %kg COD m^?3
%...........................................
l.k_m_ac = 8; %d^?1
l.K_S_ac = 0.15; %kg COD m^?3
l.K_I_nh3 = 0.0018; %M
l.pH_UL_ac = 7; % ?
l.pH_LL_ac = 6; % ?
%...........................................
l.k_m_h2 = 35; %d^?1
l.K_S_h2 = 7e-6; %kg COD m^?3
l.pH_UL_h2 = 6; % ?
l.pH_LL_h2 = 5; % ?
%...........................................
l.k_dec_X_su = 0.02; %d^?1
l.k_dec_X_aa = 0.02; %d^?1
l.k_dec_X_fa = 0.02; %d^?1
l.k_dec_X_c4 = 0.02; %d^?1
l.k_dec_X_pro= 0.02; %d^?1
l.k_dec_X_ac = 0.02; %d^?1
l.k_dec_X_h2 = 0.02; %d^?1
%...........................................
%Extra input to ADM1 model (added more substartes)
%...........................................
l.k_m_lac_f = 16; %d^?1
l.K_S_lac_f = 3.5169; %kg COD m^?3
l.k_m_lac_o = 16; %d^?1
l.K_S_lac_o = 0.6432; %kg COD m^?3
l.K_I_h2_lac_o = 1.4e-4; %kg COD m^?3
l.K_I_vfa = 3.5; %kg COD m^?3
%...........................................
l.k_dec_X_lac_f = 0.02; %d^?1
l.k_dec_X_lac_o = 0.02; %d^?1
%...........................................
l.K_S_p_caco3 = exp(-0.01183*l.T_op-8.03)/0.665^2; %kg CaCO3 m^?3 (equilibrium constan
l.K_r_caco3 = 1477.44; %d^?1 (specific rate precipitation)
%...........................................
%==========================================================================
%Steady?state input values
%==========================================================================
%Parameter Value Unit State No.
%.................................................
l.S_su_in = 0.01; %kg COD m^?3 1
l.S_aa_in = 0.001; %kg COD m^?3 2
l.S_fa_in = 0.001; %kg COD m^?3 3
l.S_va_in = 0.001; %kg COD m^?3 4
l.S_bu_in = 0.001; %kg COD m^?3 5
l.S_pro_in = 0.001; %kg COD m^?3 6
l.S_ac_in = 0.001; %kg COD m^?3 7
l.S_h2_in = 1.0e-8; %kg COD m^?3 8
l.S_ch4_in = 1.0e-5; %kg COD m^?3 9
l.S_IC_in = 0.04; %kmole C m^?3 10
l.S_IN_in = 0.01; %kmole N m^?3 11
l.S_I_in = 0.02; %kg COD m^?3 12
%.................................................
l.X_xc_in = 2.0; %kg COD m^?3 13
l.X_ch_in = 5.0; %kg COD m^?3 14
l.X_pr_in = 20.0; %kg COD m^?3 15
l.X_li_in = 5.0; %kg COD m^?3 16
l.X_su_in = 0.0; %kg COD m^?3 17
l.X_aa_in = 0.01; %kg COD m^?3 18
l.X_fa_in = 0.01; %kg COD m^?3 19
l.X_c4_in = 0.01; %kg COD m^?3 20
l.X_pro_in = 0.01; %kg COD m^?3 21
l.X_ac_in = 0.01; %kg COD m^?3 22
l.X_h2_in = 0.01; %kg COD m^?3 23
l.X_I_in = 25.0; %kg COD m^?3 24
%.................................................
l.S_cat_in = 0.04; %kmole m^?3 25
l.S_an_in = 0.02; %kmole m^?3 26
%.................................................
l.q_in = 170; %m^3 d^?1 ?
%l.T_op = 35.0; %degree C ? (defined above in Kelvin)
%.................................................
%Extra input to ADM1 model (added more substartes)
%.................................................
l.S_lac_in = 0;%0.001; %kg COD m^?3 36
%.................................................
l.X_lac_f_in= 0;%0.01; %kg COD m^?3 37
l.X_lac_o_in= 0;%0.01; %kg COD m^?3 38
%.................................................
l.S_ca_in = 0;%0.01; %kg CaCO3 m^?3 39