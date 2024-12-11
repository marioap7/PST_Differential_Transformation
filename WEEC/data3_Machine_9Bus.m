% Two Area Test Case 

% bus data format
% bus: 
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu)
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%       bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu
%       1    2        3       4       5       6     7     8    9  10  11    12
bus = [ 1  1.040     0.0     0.716   0.27    0.0   0.0   0.0  0.0  1;
	    2  1.0250    9.4869  1.63    0.067   0.0   0.0   0.0  0.0  2;
	    3  1.0250    4.7742  0.85    0.109   0.0   0.0   0.0  0.0  2;
        4  1.0105   -2.2601  0.0     0.0     0.0   0.0   0.0  0.0  3;
	    5  0.9728   -4.0604  0.0     0.0     1.25  0.5   0.0  0.0  3;
	    6  0.9890   -3.7090  0.0     0.0     0.9   0.3   0.0  0.0  3;
	    7  1.0114    3.8471  0.0     0.0     0.0   0.0   0.0  0.0  3;
	    8  0.9972    0.7834  0.0     0.0     1.0   0.35  0.0  0.0  3;
        9  1.0180    2.0381  0.0     0.0     0.0   0.0   0.0  0.0  3];

% line data format
% line: from bus(q), to bus(2), resistance(pu)(3), reactance(pu)(4),
%       line charging(pu)(5) el programa lo divide automaticamente entre 2, verificar, 
%tap ratio(6), tap phase(7)
%       1    2   3         4       5        6    7 
line = [1    4   0.0     0.0576   2*0.00    1.0  0;
        2    7   0.0     0.0625   2*0.00    1.0  0;
        3    9   0.0     0.0586   2*0.00    1.0  0;
        4    5   0.01    0.085    2*0.088   0.0  0;
        4    6   0.017   0.092    2*0.079   0.0  0;
        5    7   0.032   0.161    2*0.153   0.0  0;
        6    9   0.039   0.17     2*0.179   0.0  0;
        7    8   0.0085  0.072    2*0.0745  0.0  0;
        8    9   0.0119  0.1008   2*0.1045  0.0  0];

% Machine data format
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant T"_do(sec),

%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant T"_qo(sec),

%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. bus number
%
% note: all the following machines use sub-transient model
%          1  2  3   4    5     6       7       8    9    10     11     12 
mac_con = [1  1 100 0.0  0.0  0.146   0.0608   0.0  8.96  0.0  0.0969  0.0969  0.0  0.31  0.0 23.64  0.0125  0  1;
           2  2 100 0.0  0.0  0.8958  0.1198   0.0  6.0   0.0  0.8645  0.1969  0.0  0.535 0.0  6.4   6.8e-3  0  2;
           3  3 100 0.0  0.0  1.3125  0.1813   0.0  5.89  0.0  1.2578  0.25    0.0  0.6   0.0  3.01  4.8e-3  0  3];

% Exciter data format
% exciter:  1. exciter type - 3 for ST3
%           2. machine number
%           3. input filter time constant T_R 
%           4. voltage regulator gain K_A
%           5. voltage regulator time constant T_A
%           6. voltage regulator time constant T_B
%           7. voltage regulator time constant T_C
%           8. maximum voltage regulator output V_Rmax
%           9. minimum voltage regulator output V_Rmin
%          10. maximum internal signal V_Imax
%          11. minimum internal signal V_Imin
%          12. first stage regulator gain K_J
%          13. potential circuit gain coefficient K_p
%          14. potential circuit phase angle theta_p
%          15. current circuit gain coefficient K_I
%          16. potential source reactance X_L
%          17. rectifier loading factor K_C
%          18. maximum field voltage E_fdmax 
%          19. inner loop feedback constant K_G
%          20. maximum inner loop voltage feedback V_Gmax
%          1 2  3   4    5    6     7     8     9     10    11    12
exc_con = [0 1 0 20 0.2  20   0.06  0     0    1.0   -0.9  0.0  0.46   3.1   0.33  2.3  0.1   0.1   1.0    0    0   0;
           0 2 0 20 0.2  20   1.0   6.67  1.0 10.0  -10.0  0.2  -0.2 200.0   4.37 20    4.83  0.09  1.1    8.63 1   6.53;
           0 3 0 20 0.2  20   0     0     0    5.0   -5.0  0     0     0     0     0    0     0     0      0    0   0];
           %2 4 0.01 300.0   0.01  0     0    4.95  -4.9  1.0   1.33  3.05  0.279 2.29 0.117 0.1   0.675  0    0   0];

%pss on exciter 2 and 3
pss_con = [1 2 300  20.0  0.1   0.02  0.1   0.02  0.2 -0.05;
           1 3 300  20.0  0.06  0.04  0.08  0.04  0.2 -0.05];

% governor model
% tg_con matrix format
%column	       data			unit
%  1	turbine model number (=1)	
%  2	machine number	
%  3	speed set point   wf		pu
%  4	steady state gain 1/R		pu
%  5	maximum power order  Tmax	pu on generator base
%  6	servo time constant   Ts	sec
%  7	governor time constant  Tc	sec
%  8	transient gain time constant T3	sec
%  9	HP section time constant   T4	sec
% 10	reheater time constant    T5	sec

tg_con = [1  1  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
          1  2  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
          1  3  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
          1  4  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0];


% non-conforming load
% col 1           bus number
% col 2           fraction const active power load
% col 3           fraction const reactive power load
% col 4           fraction const active current load
% col 5           fraction const reactive current load


disp('0.5 constant current load, svc at bus 101')

%svc
% col 1           svc number
% col 2           bus number
% col 3           svc base MVA
% col 4           maximum susceptance Bcvmax(pu)
% col 5           minimum susceptance Bcvmin(pu)
% col 6           regulator gain
% col 7		  regulator time constant (s)

svc_con = [1  101  600  1  0  100  0.05];

%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  - 0 three phase
%                            - 1 line to ground
%                            - 2 line-to-line to ground
%                            - 3 line-to-line
%                            - 4 loss of line with no fault
%                            - 5 loss of load at bus
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)

sw_con = [  0   0    0    0    0    0    0.01;%sets intitial time step
          0.1   3  101    0    0    0    0.005; %apply three phase fault at bus 3, on line 3-101
          0.15  0    0    0    0    0    0.005556; %clear fault at bus 3
          0.20  0    0    0    0    0    0.005556; %clear remote end
          0.50  0    0    0    0    0    0.01; % increase time step 
          1.0   0    0    0    0    0    0.01; % increase time step
         10.0   0    0    0    0    0    0.0]; % end simulation
     


