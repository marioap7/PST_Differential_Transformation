% m.file to simulate the 16 machine, 68 bus system
% using the Matlab Power System Toolbox

clear; clear global; clc;

global  basmva basrad syn_ref mach_ref
global  bus_v psi_re psi_im cur_re cur_im exc_con

% synchronous machine variables
global  mac_con mac_ang mac_spd eqprime edprime vex
global  pmech dmac_ang dmac_spd deqprime dedprime Efd dEfd

global Ng


d46_mexico_clas_4o
%data3_Machine_9Bus

mac_con(:,20:21)=0;
mac_con(:,22:23)=1;



% solve for loadflow - loadflow parameter
  tol = 1e-4;   % tolerance for convergence
  iter_max = 15; % maximum number of iterations
  vmin = 0.5;   % voltage minimum
  vmax = 1.5;   % voltage maximum
  acc = 1.0;   % acceleration factor

  [bus_sol,line_flw] = loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2);
  bus = bus_sol;
  
jay = sqrt(-1);

basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base
syn_ref = 0 ;     % synchronous reference frame

Ng = length(mac_con(:,1));
mac_em_idx = 1:Ng; %Maquinas con modelo clasico
Mod = 2*ones(1,Ng);

disp('Performing simulation.')

% simulation
t_switch(1) = 0;     % all time in second+s, start time
t_switch(2) = 1.0;  % time to apply fault
t_switch(3) = 1.0 + 3/60; % time to clear fault, 3 cycles
t_switch(4) = 20.0;   % end time

h = 1/60; % integration stepsize

k_switch(1) = round((t_switch(2) - t_switch(1))/h) + 1;
k_switch(2) = round((t_switch(3) - t_switch(1))/h) + 1;
k_switch(3) = round((t_switch(4) - t_switch(1))/h) + 1;


% step 1: construct reduced Y matrix 
[Y_red,V_rec] = red_ybus(bus,line); %pre-fault admittance matrix

% create bus matrix with load increased on node 28
bus_f = bus;
bus_f(5,9) = -150.0;
%bus_f(5,9) = -150.0;
[Y_red_f,V_rec_f] = red_ybus(bus_f,line);   % fault-on admittance matrix


% line_pf = line;
% nline = length(line(:,1));
% line_pf(nline+1,:)=[32	33	-0.0007	-0.0089	-0.1342	0	0 ];% remove line 7-6
% %line_pf(nline+1,:) = [5 7 -0.032 -0.161 -0.153 0 0]; %remove line 5-7
% [Y_red_pf,V_rec_pf] = red_ybus(bus,line_pf);  % post-fault

%=========================================================
% Calculo de las condiciones iniciales

flag = 0;

for i = 1:Ng
    mac_tra(i,1,bus,flag);
    smpexc(i,1,bus,flag,h);
    
end
%=========================================================

% step 3: perform a predictor-corrector integration
tic
for k = 1:k_switch(3)
    % step 3a: network solution
    mach_ref(k) = 0; % No se toma maquina de referencia.
%============================================
    flag = 1;
    
    mac_tra(0,k,bus,flag);
    smpexc(0,k,bus,flag,h);

%===============================================
    
    % network-machine interface
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
    
    
    mac_tra(0,k,bus,flag);
    smpexc(0,k,bus,flag,h);
    
    if k ~=k_switch(3)
        j = k+1;
        % following statements are predictor steps
        mac_ang(:,j) = mac_ang(:,k) + h*dmac_ang(:,k);
        mac_spd(:,j) = mac_spd(:,k) + h*dmac_spd(:,k);
        edprime(:,j) = edprime(:,k) + h*dedprime(:,k);
        eqprime(:,j) = eqprime(:,k) + h*deqprime(:,k);
        Efd(:,j)     = Efd(:,k)      + h*dEfd(:,k);
        
        
%======================================================
        flag = 1;
        mach_ref(j) = 0;
         %  machines with no exciters    
        vex(:,j) = vex(:,k);    % field voltage not changed for

        
        mac_tra(0,j,bus,flag);
        smpexc(0,j,bus,flag,h);
   %===================================================
        
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
        smpexc(0,j,bus,flag,h);

        
        % following statements are corrector steps
        mac_ang(:,j) = mac_ang(:,k) + h*(dmac_ang(:,k)+dmac_ang(:,j))/2.;
        mac_spd(:,j) = mac_spd(:,k) + h*(dmac_spd(:,k)+dmac_spd(:,j))/2.;
        edprime(:,j) = edprime(:,k) + h*(dedprime(:,k)+dedprime(:,j))/2.;
        eqprime(:,j) = eqprime(:,k) + h*(deqprime(:,k)+deqprime(:,j))/2.;
        Efd(:,j)     = Efd(:,k)     + h*(dEfd(:,k)    +dEfd(:,j))/2;
    end
end
toc
t = 0:h:t_switch(4); % time

delta = (mac_ang - mac_ang(1,:))*180/pi;
w = mac_spd;
% figure(1)
% plot(t,delta(1,:))
% xlabel('tiempo (s)')
% ylabel('\delta (grados)')
% grid

%figure(2)
plot(t,delta(2,:))
xlabel('tiempo (s)')
ylabel('velocidad (pu)')
grid
