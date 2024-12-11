function [bus_sol,line_flow] = ...
loadflow(bus,line,tol,iter_max,vmin,vmax,acc,display,flag)
% Syntax:    [bus_sol,line_flow] =
% loadflow(bus,line,tol,iter_max,vmin,vmax,acc,display,flag)
%
% Purpose:   solve the load-flow equations of power systems
%            modified to eliminate do loops and improve the use 
%            sparse matices
% Input:    bus       - bus data
%           line      - line data
%           tol       - tolerance for convergence
%           iter_max  - maximum number of iterations
%           vmin      - voltage minimum limit
%           vmax      - voltage maximum limit
%           acc       - acceleration factor
%           display   - 'y', generate load-flow study report
%                        else, no load-flow study report
%           flag      - 1, form new Jacobian every iteration
%                       2, form new Jacobian every other 
%                           iteration

% Output:   bus_sol   - bus solution (see report for the
%                         solution format)
%           line_flow - line flow solution (see report)
%
% See also:  
%
% Algorithm: Newton-Raphson method using the polar form of 
%   the equations for P(real power) and Q(reactive power).
%
% Calls:     Y_sparse, calc, form_jac
%

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
% Version:   2.0
% Author:     Graham Rogers
% Date:         March 1994
% Version:   1.0
% Authors:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%
% ***********************************************************
global bus_int Y

tt = clock;     % start the total time clock
jay = sqrt(-1);
load_bus = 3;
gen_bus = 2;
swing_bus = 1;
if exist('flag') == 0
  flag = 1;
end
if flag <1 | flag > 2
  error('LOADFLOW: flag not recognized')
end
nline = length(line(:,1));     % number of lines
nbus = length(bus(:,1));     % number of buses
% set maximum and minimum voltage
volt_min = vmin*ones(nbus,1);
volt_max = vmax*ones(nbus,1);
% build admittance matrix Y
%disp('building Y matrix')

[Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line);
%disp('setting up for Jacobian calculation')
% process bus data
bus_no = bus(:,1);
V = bus(:,2);
ang = bus(:,3)*pi/180;
Pg = bus(:,4);
Qg = bus(:,5);
Pl = bus(:,6);
Ql = bus(:,7);
Gb = bus(:,8);
Bb = bus(:,9);
bus_type = round(bus(:,10));
sw_bno=ones(nbus,1);
g_bno=sw_bno;
% set up index for Jacobian calculation
%% form PQV_no and PQ_no
bus_zeros=zeros(nbus,1);
bus_index=[1:1:nbus]';
swing_index=find(bus_type==1); 
sw_bno(swing_index)=bus_zeros(swing_index);
PQV_no=find(bus_type >=2);
PQ_no=find(bus_type==3);
gen_index=find(bus_type==2);
g_bno(gen_index)=bus_zeros(gen_index);    
%sw_bno is a vector having ones everywhere but the swing bus locations
%g_bno is a vector having ones everywhere but the generator bus locations
% construct sparse angle reduction matrix
il = length(PQV_no);
ii = [1:1:il]';
ang_red = sparse(ii,PQV_no,ones(il,1),il,nbus);
% construct sparse voltage reduction matrix
il = length(PQ_no);
ii = [1:1:il]';
volt_red = sparse(ii,PQ_no,ones(il,1),il,nbus);
iter = 0;     % initialize iteration counter
% calculate the power mismatch and check convergence

  [delP,delQ,P,Q,conv_flag] =...
             calc(nbus,V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol);
%keyboard

st = clock;     % start the iteration time clock
%% start iteration process
while (conv_flag == 1 & iter < iter_max)
iter = iter + 1;
% Form the Jacobean matrix
   if flag == 2
      if iter == 2*fix(iter/2) + 1
       clear Jac
       Jac = form_jac(V,ang,Y,ang_red,volt_red);
      end
     else
       clear Jac
	Jac=form_jac(V,ang,Y,ang_red,volt_red);
   end
% reduced mismatch real and reactive power vectors
  red_delP = ang_red*delP;
  red_delQ = volt_red*delQ;
  clear delP delQ
% solve for voltage magnitude and phase angle increments
  temp = Jac\[red_delP; red_delQ];
% expand solution vectors to all buses
  delAng = ang_red'*temp(1:length(PQV_no),:);
  delV = volt_red'*temp(length(PQV_no)+1:length(PQV_no)+length(PQ_no),:);
% update voltage magnitude and phase angle
  V = V + acc*delV;
  V = max(V,volt_min);  % voltage higher than minimum
  V = min(V,volt_max);  % voltage lower than maximum
  ang = ang + acc*delAng;
% calculate the power mismatch and check convergence
%  t_s = clock;
  [delP,delQ,P,Q,conv_flag] =...
             calc(nbus,V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol);
%  etime(clock,t_s)
end; %%
ste = clock;     % end the iteration time clock
gen_index=find(bus_type==2);
load_index = find(bus_type==3);
    Pg(gen_index) = P(gen_index) + Pl(gen_index);
    Qg(gen_index) = Q(gen_index) + Ql(gen_index);
    Pl(load_index) = Pg(load_index) - P(load_index);
    Ql(load_index) = Qg(load_index) - Q(load_index);
Pg(SB) = P(SB) + Pl(SB); Qg(SB) = Q(SB) + Ql(SB);
VV = V.*exp(jay*ang);  % solution voltage 
% calculate the line flows and power losses
tap_index = find(abs(line(:,6))>0);
tap_ratio = ones(nline,1);
tap_ratio(tap_index)=line(tap_index,6);
phase_shift(:,1) = line(:,7);
tps = tap_ratio.*exp(jay*phase_shift*pi/180);
from_bus = line(:,1);
from_int = bus_int(round(from_bus));
to_bus = line(:,2);
to_int = bus_int(round(to_bus));
r = line(:,3);
rx = line(:,4);
chrg = line(:,5);
z = r + jay*rx;
y = ones(nline,1)./z;


MW_s = VV(from_int).*conj((VV(from_int) - tps.*VV(to_int)).*y ...
       + VV(from_int).*(jay*chrg/2))./(tps.*conj(tps));
P_s = real(MW_s);     % active power sent out by from_bus
                      % to to_bus
Q_s = imag(MW_s);     % reactive power sent out by 
                      % from_bus to to_bus
MW_r = VV(to_int).*conj((VV(to_int) ...
       - VV(from_int)./tps).*y ...
       + VV(to_int).*(jay*chrg/2));
P_r = real(MW_r);     % active power received by to_bus 
                      % from from_bus
Q_r = imag(MW_r);     % reactive power received by 
                      % to_bus from from_bus
iline = [1:1:nline]';
  line_ffrom = [iline from_bus to_bus P_s Q_s];
  line_fto   = [iline to_bus from_bus P_r Q_r];
% keyboard
P_loss = sum(P_s) + sum(P_r) ;
Q_loss = sum(Q_s) + sum(Q_r) ;
bus_sol=[bus_no  V  ang*180/pi Pg Qg Pl Ql Gb Bb bus_type];
line_flow = [line_ffrom; line_fto];
% display results
if display == 'y',
  clc
  disp('                             LOAD-FLOW STUDY')
  disp('                    REPORT OF POWER FLOW CALCULATIONS ')
  disp(' ')
  disp(date)
  fprintf('SWING BUS                  : BUS %g \n', SB)
  fprintf('NUMBER OF ITERATIONS       : %g \n', iter)
  fprintf('SOLUTION TIME              : %g sec.\n',etime(ste,st))
  fprintf('TOTAL TIME                 : %g sec.\n',etime(clock,tt))
  fprintf('TOTAL REAL POWER LOSSES    : %g.\n',P_loss)
  fprintf('TOTAL REACTIVE POWER LOSSES: %g.\n\n',Q_loss)
  if conv_flag == 0,
    disp('                                      GENERATION             LOAD')
    disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
    disp(bus_sol(:,1:7))

    disp('                      LINE FLOWS                     ')
    disp('      LINE  FROM BUS    TO BUS      REAL  REACTIVE   ')
    disp(line_ffrom)
    disp(line_fto)
  end
end; %

%iter
if iter > iter_max,
  fprintf('Note: Solution did not converge in %g iterations.\n', iter_max)
end

return

