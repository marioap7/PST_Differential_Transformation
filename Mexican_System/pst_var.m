% file: pst_var.m
%
% Syntax: pst_var
%
% Purpose: Define global variables for power system 
%          simulation.

%
% Input:
%
% Output:
%
% Files:
%
% See Also:

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     January 1991

% system variables
global  basmva basrad syn_ref mach_ref sys_freq Mod
global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot mac_int
global  mac_ang mac_spd eqprime edprime psikd psikq
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq 
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime dpsikd dpsikq n_em

% excitation system variables
global  exc_con exc_pot Efd V_R V_A V_As R_f V_FB V_TR V_B
global  dEfd dV_R dV_As dR_f dV_TR
global  exc_sig

% non-conforming load variables
global  load_con load_pot

% svc variables
global  svc_con svc_pot B_cv dB_cv
global  svc_sig

% pss variables
global  pss_con pss_pot 
global  pss1 pss2 pss3 dpss1 dpss2 dpss3

% turbine-governor variables
global  tg_con tg_pot 
global  tg1 tg2 tg3 dtg1 dtg2 dtg3 
