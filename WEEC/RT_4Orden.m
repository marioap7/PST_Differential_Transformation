clc; clear;

global k_switch Ng k Y_red Y_red_f


jay = sqrt(-1);
pst_var % set up global variable
basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base

data3_Machine_9Bus; % Datos de entrada

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
bus_f(5,9) = -100;
[Y_red_f,V_rec_f] = red_ybus(bus_f,line);   % fault-on admittance matrix


k_switch(1) = round((t_switch(2)-t_switch(1))/h)+1;
k_switch(2) = round((t_switch(3)-t_switch(1))/h)+1;
k_switch(3) = round((t_switch(4)-t_switch(1))/h)+1;

%Calculo de las condiciones iniciales

[Ng,dum] = size(mac_con);
if dum < 23                    % set power fraction
   mac_con = [mac_con ones(Ng,2)];  % to unity
end

busnum = bus_int(mac_con(:,2)); % bus number 
mac_pot(:,1) = basmva*ones(Ng,1)./mac_con(:,3);% scaled MVA base                     
mac_pot(:,2) = ones(Ng,1); % base kv
      
      % extract bus information
eterm = bus(busnum,2);  % terminal bus voltage
theta = bus(busnum,3)*pi/180; % terminal bus angle in radians 
                          
qelect = bus(busnum,5);%.*mac_con(:,23);% electrical output power, reactive                            
pelect = bus(busnum,4); %.*mac_con(:,22); % electrical output power, active 
                                             
curr = sqrt(pelect.^2+qelect.^2)./eterm.*mac_pot(:,1);%current magnitude            
phi = atan2(qelect,pelect); % power factor angle
                                        
v = eterm.*(cos(theta)+jay*sin(theta));%voltage in real and imaginary parts
curr = curr.*(cos(theta-phi) + jay*sin(theta-phi)); 
 
eprime = v + jay*mac_con(:,7).*curr; 
ei = v + jay*mac_con(:,11).*curr;
      
mac_ang = atan2(imag(ei),real(ei)); % machine angle (delta)
mac_spd = ones(Ng,1); % machine speed at steady state
                            
rot = sin(mac_ang)+jay*cos(mac_ang);% system reference frame rotation                           
eprime = eprime.*rot;
edprime = real(eprime); 
eqprime = imag(eprime); 
      
curr = curr.*rot;  
curdg = real(curr); 
curqg = imag(curr);
curd = real(curr)./mac_pot(:,1); 
curq = imag(curr)./mac_pot(:,1);

v = v.*rot;
ed = real(v); 
eq = imag(v);
 
vex = eqprime + (mac_con(:,6) - mac_con(:,7)).*curdg;
pmech = pelect.*mac_pot(:,1); % set input
      
efd = vex;
V_A = efd./exc_con(:,4); % laglead

exc_pot(:,4) = exc_con(:,7)./exc_con(:,5);
exc_pot(:,3) = eterm + V_A; % reference voltage
exc_pot(:,5) = ones(Ng,1);

y0 = [edprime eqprime mac_ang mac_spd efd];
v0 = [pelect curd curq ed eq V_A];

tol = 1e-4;

% tic
t = t0;
y = y0;
v = v0;
tout = t;
N = tf/h;

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
        if iter > 1000
            disp('Demasidas iteraciones')
            break
        end
    end
    % Update the solution
    t = t + h;
    y = y1;
    tout(k) = t;
    epd(:,k) = y1(:,1);
    epq(:,k) = y1(:,2);
    delta(:,k) = y1(:,3);
    Mac_spd(:,k) = y1(:,4);
    Efd(:,k) = y1(:,5);

end
toc
t = tout;
w = Mac_spd;
delta = (delta - delta(1,:))*180/pi;


% plot(t,w)
% xlabel('Tiempo (s)')
% ylabel('Velocidad (pu)')
% grid