clc; clear; clear global;

global k_switch Ng k Y_red Y_red_f modelo


jay = sqrt(-1);
pst_var % set up global variable
basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base

d46_mexico_clas_4o % Datos de entrada
% data16m_4o
Ng = length(mac_con(:,1));
xpd = mac_con(1:Ng,7);
H = mac_con(:,16);
D = mac_con(:,17);


disp('Performing simulation.')

% simulation
t_switch(1) = 0;
t_switch(2) = 1.0;  % time to apply fault
t_switch(3) = 1.0 + 3/60; % time to clear fault, 3 cycles
t_switch(4) = 20; 

h = 1/60; % integration stepsize

t0 = t_switch(1);
tf = t_switch(4);

% solve for loadflow - loadflow parameter
  tol = 1e-4;   % tolerance for convergence
  iter_max = 15; % maximum number of iterations
  vmin = 0.5;   % voltage minimum
  vmax = 1.5;   % voltage maximum
  acc = 1.0;   % acceleration factor

  [bus_sol,line_flw] = loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2);
  bus = bus_sol;

  % step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
[Y_red,V_rec] = red_ybus(bus,line); % pre-fault admittance matrix

% create bus matrix with load increased on node 28
bus_f = bus;
bus_f(5,9) = -150;
[Y_red_f,V_rec_f] = red_ybus(bus_f,line);   % fault-on admittance matrix

 
k_switch(1) = round((t_switch(2)-t_switch(1))/h)+1;
k_switch(2) = round((t_switch(3)-t_switch(1))/h)+1;
k_switch(3) = round((t_switch(4)-t_switch(1))/h)+1;

%Calculo de las condiciones iniciales

[Ng,dum] = size(mac_con);
if dum < 23                    % set power fraction
   mac_con = [mac_con ones(Ng,2)];  % to unity
end

%modelo = [4*ones(1,8) 2 4 4 2 2 2 2 4*ones(1,5) 2 2 2 2 4 4 4 2 2 2 4*ones(1,10) 2 4 2 4 4 4];
modelo = 4*ones(1,Ng);
for i = 1:Ng
    
    busnum = bus_int(mac_con(i,2)); % bus number
    mac_pot(i,1) = basmva/mac_con(i,3); % scaled MVA base
    mac_pot(i,2) = 1.0; % base kv
    
    eterm(i,1) = bus(busnum,2);  % terminal bus voltage
    theta(i,1) = bus(busnum,3)*pi/180; %terminal bus angle in radians 
    
    pelect(i,1) = bus(busnum,4)*mac_con(i,22); %electrical output power, active               
    qelect(i,1) = bus(busnum,5)*mac_con(i,23); % electrical output power, reactive
    
    curr = sqrt(pelect(i,1)^2+qelect(i,1)^2)/eterm(i,1)*mac_pot(i,1);
    phi = atan2(qelect(i,1),pelect(i,1)); %power factor angle
      %current magnitude on generator base
      
     % voltage in real and imaginary parts on system reference frame                               
    v = eterm(i,1)*exp(jay*theta(i,1)); 
    curr = curr*exp(jay*(theta(i,1)-phi)); % complex current  
      
        if mac_con(i,14) == 0
            mac_con(i,14) = 999.0;
        end
        eprime = v + jay*mac_con(i,7)*curr; 
        ei = v + jay*mac_con(i,11)*curr;
      
        mac_ang(i,1) = atan2(imag(ei),real(ei)); % machine angle (delta)
        mac_spd(i,1) = 1; % machine speed at steady state
                            
        rot = sin(mac_ang(i,1)) + jay*cos(mac_ang(i,1));% system reference frame rotation                           
        eprime = eprime*rot;
        edprime(i,1) = real(eprime); 
        eqprime(i,1) = imag(eprime); 
      
        curr = curr.*rot;  
        curdg(i,1) = real(curr); 
        curqg(i,1) = imag(curr);
        curd(i,1) = real(curr)/mac_pot(i,1); 
        curq(i,1) = imag(curr)/mac_pot(i,1);

        v = v*rot;
        ed(i,1) = real(v); 
        eq(i,1) = imag(v);
 
        %E_Isat = eqprime(i,1); % select higher voltage
        vex(i,1) = eqprime(i,1) + (mac_con(i,6) - mac_con(i,7)).*curdg(i,1);
        pmech(i,1) = pelect(i,1)*mac_pot(i,1); % set input   
        efd(i,1) = vex(i,1);
        
        V_A(i,1) = efd(i,1)/exc_con(i,4); % laglead

        exc_pot(i,4) = exc_con(i,7)/exc_con(i,5);
        exc_pot(i,3) = eterm(i,1) + V_A(i,1); % reference voltage
        exc_pot(i,5) = 1;
end




y0 = [edprime eqprime mac_ang mac_spd efd];
v0 = [pelect curd curq ed eq V_A];

tol = 1e-4;


t = t0;
y = y0;
v = v0;
tout = t;
N = tf/h+1;

% The main loop
tic
k = 1;
for k = 1:N
    % Compute the slopes
    [temp, temp1] = feval('fun_Generador_RT',t,y,v);
    f1 = temp;
    v = temp1;
    y1 = y;
    erroRT = 1e2;
    iter = 1;
    while erroRT > tol
        [temp, temp1] = feval('fun_Generador_RT', t+h, y1,v);
        f2 = temp; 
        ya = y1;
        y1 = y + h*(f1 + f2)/2;
        
        erroRT = max(max(abs(y1 - ya)));
        iter = iter + 1;
    end
    % Update the solution
    t = t + h;
    y = y1;
    tout(k) = t;
    epd(:,k) = y1(:,1);
    epq(:,k) = y1(:,2);
    delta(:,k) = y1(:,3);
    Mac_spd(:,k) = y1(:,4);
    efd(:,k) = y1(:,5);

end
toc
w = Mac_spd;
delta = (delta - delta(1,:))*180/pi;
t = tout;
plot(t,delta(2,:))
title('Velocidad')
xlabel('tiempo (s)')
ylabel('Magnitud (pu)')
grid