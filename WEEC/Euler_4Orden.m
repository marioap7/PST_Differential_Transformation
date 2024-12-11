% m.file to simulate the 16 machine, 68 bus system
% using the Matlab Power System Toolbox

clear; clear global; clc

jay = sqrt(-1);
pst_var % set up global variable

%Matriz bus y line
data3_Machine_9Bus;


basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base
syn_ref = 0 ;     % synchronous reference frame

% solve for loadflow - loadflow parameter
  tol = 1e-4;   % tolerance for convergence
  iter_max = 15; % maximum number of iterations
  vmin = 0.5;   % voltage minimum
  vmax = 1.5;   % voltage maximum
  acc = 1.0;   % acceleration factor

  [bus_sol,line_flw] = loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2);
  bus = bus_sol;

disp('Performing simulation.')

% simulation
t_switch(1) = 0;     % all time in second+s, start time
t_switch(2) = 1.0;  % time to apply fault 0.08
t_switch(3) = 1.0 + 3/60; % time to clear fault, 3 cycles 0.11
t_switch(4) = 20;   % end time 
h = 1/60; % integration stepsize

k_switch(1) = round((t_switch(2)-t_switch(1))/h)+1;
k_switch(2) = round((t_switch(3)-t_switch(1))/h)+1;
k_switch(3) = round((t_switch(4)-t_switch(1))/h)+1;

% step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
[Y_red,V_rec] = red_ybus(bus,line); % pre-fault admittance matrix

% create bus matrix with load increased on node 28
bus_f = bus;
bus_f(5,9) = -100;
[Y_red_f,V_rec_f] = red_ybus(bus_f,line);   % fault-on admittance matrix

% step 2: initialization

flag = 0;

  mac_tra(0,1,bus,flag);    % machine model
  smpexc(0,1,bus,flag,h);   % exciter model

% step 3: perform a predictor-corrector integration 
tic
for k = 1:k_switch(3)
  % step 3a: network solution
  mach_ref(k) = 0; % No se toma maquina de referencia.

  flag = 1;

  mac_tra(0,k,bus,flag);    % network-machine interface
  smpexc(0,k,bus,flag,h);
  psi = psi_re(:,k) + jay*psi_im(:,k); % voltajes internos en marco Re-Im

  if k >= k_switch(2)
     cur = Y_red*psi; % network solution corrientes de generadores POSTFALLA
      bus_v(:,k) = V_rec*psi; % bus voltage reconstruction      
  elseif k >= k_switch(1)
     cur = Y_red_f*psi;  % network solution FALLA
      bus_v(:,k) = V_rec_f*psi; % bus voltage reconstruction      
  else

     cur = Y_red*psi; % network solution PREFALLA
      bus_v(:,k) = V_rec*psi; % bus voltage reconstruction      
  end
  
  cur_re(:,k) = real(cur); 
  cur_im(:,k) = imag(cur);
  % step 3b: compute dynamics and integrate
  
  flag = 2;

  pmech(:,k) = pmech(:,1); % constant mechanical input power
  mac_tra(0,k,bus,flag); % dynamics calculation
  [nexc dum] = size(exc_con);
  exc_sig(:,k) = zeros(nexc,1);
  smpexc(0,k,bus,flag,h);

  if k ~=k_switch(3)
    j = k+1;
    % following statements are predictor steps
    mac_ang(:,j) = mac_ang(:,k) + h*dmac_ang(:,k); 
    mac_spd(:,j) = mac_spd(:,k) + h*dmac_spd(:,k);
    edprime(:,j) = edprime(:,k) + h*dedprime(:,k);
    eqprime(:,j) = eqprime(:,k) + h*deqprime(:,k);
    Efd(:,j)     = Efd(:,k)     + h*dEfd(:,k);


    flag = 1;
    mach_ref(j) = 0;
    vex(:,j) = vex(:,k);% field voltage not changed for machines with no exciters
    mac_tra(0,j,bus,flag); 
    smpexc(0,j,bus,flag,h);
    psi = psi_re(:,j) + jay*psi_im(:,j);
    
    if k >= k_switch(2)
      cur = Y_red*psi; % POSTFALLA
    elseif k >= k_switch(1)
      cur = Y_red_f*psi; % FALLA
    else
      cur = Y_red*psi; % PREFALLA
    end
    
    cur_re(:,j) = real(cur); 
    cur_im(:,j) = imag(cur);
    pmech(:,j) = pmech(:,k);

    flag = 2;
    mac_tra(0,j,bus,flag);
    exc_sig(:,j) = zeros(nexc,1);
    smpexc(0,j,bus,flag,h);

    % following statements are corrector steps
    mac_ang(:,j) = mac_ang(:,k) + h*(dmac_ang(:,k)+dmac_ang(:,j))/2.;
    mac_spd(:,j) = mac_spd(:,k) + h*(dmac_spd(:,k)+dmac_spd(:,j))/2.;
    edprime(:,j) = edprime(:,k) + h*(dedprime(:,k)+dedprime(:,j))/2.;
    eqprime(:,j) = eqprime(:,k) + h*(deqprime(:,k)+deqprime(:,j))/2.;
    Efd(:,j)     = Efd(:,k)     +  h*(dEfd(:,k)+dEfd(:,j))/2.;
  end
end
toc
t = [t_switch(1):h:t_switch(4)]'; % time
delta = (mac_ang - mac_ang(1,:))*180/pi;
w = mac_spd;

% plot(t,w)
% xlabel('Tiempo (s)')
% ylabel('Velocidad (pu)')
% grid
