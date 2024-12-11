 % m.file to simulate the 3 machine, 9 bus system
% using the Matlab Power System Toolbox

clear; clear global; clc;

jay = sqrt(-1);
global mac_con basrad basmva bus_int
basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base

data3_Machine_9Bus
% d2a_sub
disp('Performing simulation.')

% simulation
t_switch(1) = 0;
t_switch(2) = 1.0;  % time to apply fault
t_switch(3) = 1.0 + 3/60; % time to clear fault, 3 cycles
t_switch(4) = 20;

h = 1/60; % integration stepsize
p = 2;%Numero de coeficientes
hvec = h.^(0:p);

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
  
%  for i = 1:length(mac_con(:,1))
%      if mac_con(i,3) ~= basmva
%          mac_con(i,6) = basmva*mac_con(i,6)/mac_con(i,3);
%          mac_con(i,7) = basmva*mac_con(i,7)/mac_con(i,3);
%          mac_con(i,11) = basmva*mac_con(i,11)/mac_con(i,3);
%          mac_con(i,12) = basmva*mac_con(i,12)/mac_con(i,3);
%          mac_con(i,16) = mac_con(i,16)*mac_con(i,3)/basmva;
%          mac_con(i,17) = mac_con(i,17)*mac_con(i,3)/basmva;
%          mac_con(i,3) = basmva;
%      end
%  end

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
Vt(:,1) = bus(busnum,2);  % terminal bus voltage
theta(:,1) = bus(busnum,3)*pi/180; % terminal bus angle in radians 
                          
qelect(:,1) = bus(busnum,5);%.*mac_con(:,23);% electrical output power, reactive                            
pelect = bus(busnum,4); %.*mac_con(:,22); % electrical output power, active 
                                             
curr = sqrt(pelect(:,1).^2+qelect(:,1).^2)./Vt(:,1).*mac_pot(:,1);%current magnitude            
phi = atan2(qelect(:,1),pelect(:,1)); % power factor angle
                                        
v = Vt(:,1).*(cos(theta(:,1))+jay*sin(theta(:,1)));%voltage in real and imaginary parts
curr = curr.*(cos(theta(:,1)-phi) + jay*sin(theta(:,1)-phi)); 
 
eprime = v + jay*mac_con(:,7).*curr; 
ei = v + jay*mac_con(:,11).*curr;
      
mac_ang = atan2(imag(ei),real(ei)); % machine angle (delta)
mac_spd = ones(Ng,1); % machine speed at steady state
                            
rot = sin(mac_ang(:,1))+jay*cos(mac_ang(:,1));% system reference frame rotation                           
eprime = eprime.*rot;
edprime(:,1) = real(eprime); 
eqprime(:,1) = imag(eprime); 
      
curr = curr.*rot;  
curdg(:,1) = real(curr); 
curqg(:,1) = imag(curr);
curd(:,1) = real(curr)./mac_pot(:,1); 
curq(:,1) = imag(curr)./mac_pot(:,1);

v = v.*rot;
ed(:,1) = real(v); 
eq(:,1) = imag(v);
 
E_Isat = eqprime(:,1); % select higher voltage
vex(:,1) = eqprime(:,1) + (mac_con(:,6) - mac_con(:,7)).*curdg(:,1);
pmech(:,1) = pelect(:,1).*mac_pot(:,1); % set input
      

Efd(:,1) = vex(:,1);
V_A(:,1) = Efd(:,1)./exc_con(:,4); % laglead

exc_pot(:,4) = exc_con(:,7)./exc_con(:,5);
exc_pot(:,3) = Vt(:,1) + V_A(:,1); % reference voltage
exc_pot(:,5) = ones(Ng,1);

fi = sin(mac_ang);
psi = cos(mac_ang);

delta(:,1) = mac_ang;
w(:,1) = mac_spd;
Edp(:,1) = edprime;
Eqp(:,1) = eqprime;
Ef(:,1) = Efd;
  % generator dynamics calculation
tic
for k = 1:k_switch(3)-1

    for n = 1:p
        psi_re = 0;
        psi_im = 0;
        for m = 1:n
            psi_re = psi_re + fi(:,m).*edprime(:,n-m+1) + psi(:,m).*eqprime(:,n-m+1);
            psi_im = psi_im - psi(:,m).*edprime(:,n-m+1) + fi(:,m).*eqprime(:,n-m+1);
        end
        vxy(:,n) = psi_re + jay*psi_im;
    
        % network-machine interface

        if k >= k_switch(2) 
            cur(:,n) = Y_red*vxy(:,n); % network solution corrientes de generadores POSTFALLA
            bus_v(:,n) = V_rec*vxy(:,n); % bus voltage reconstruction
        elseif k >= k_switch(1)  
            cur(:,n) = Y_red_f*vxy(:,n);  % network solution FALLA
            bus_v(:,n) = V_rec_f*vxy(:,n); % bus voltage reconstruction
        else
            cur(:,n) = Y_red*vxy(:,n); % network solution PREFALLA
            bus_v(:,n) = V_rec*vxy(:,n); % bus voltage reconstruction
        end
            
        cur_re(:,n) = real(cur(:,n)); 
        cur_im(:,n) = imag(cur(:,n));
         
        % step 3b: compute dynamics and integrate
  %pmech(:,k) = pmech(:,1); % constant mechanical input power 
 

        auxcurd = 0;
        auxcurq = 0;
        for m = 2:n
            auxcurd = auxcurd + fi(:,m).*curd(:,n-m+1) + psi(:,m).*curq(:,n-m+1);
            auxcurq = auxcurq - psi(:,m).*curd(:,n-m+1) + fi(:,m).*curq(:,n-m+1);
        end
        curd(:,n) = fi(:,1).*(cur_re(:,n)-auxcurd) - psi(:,1).*(cur_im(:,n)-auxcurq);
        curq(:,n) = psi(:,1).*(cur_re(:,n)-auxcurd) + fi(:,1).*(cur_im(:,n)-auxcurq);

        curdg(:,n) = curd(:,n).*mac_pot(:,1);
        curqg(:,n) = curq(:,n).*mac_pot(:,1);
         
        ed(:,n) = edprime(:,n) + mac_con(:,7).*curqg(:,n);
        eq(:,n) = eqprime(:,n) - mac_con(:,7).*curdg(:,n);
         
         
        auxpelect = 0;
        aux = 0;
        for m = 1:n
            auxpelect = auxpelect + ed(:,m).*curdg(:,n-m+1) + eq(:,m).*curqg(:,n-m+1);
            aux = aux + ed(:,m).*ed(:,n-m+1) + eq(:,m).*eq(:,n-m+1);
        end
        x(:,n) = aux;
        pelect(:,n) = auxpelect;
         
        if n == 1 
            Vt(:,1) = sqrt(x(:,1));
        else
            aux = 0;
            for m = 2:n-1
                aux = aux + Vt(:,m).*Vt(:,n-m+1);
            end
            Vt(:,n) = x(:,n)./(2*Vt(:,1)) - aux./(2*Vt(:,1));
        end
         
        if n == 1
            V_A(:,n) = exc_pot(:,3) - Vt(:,n);
        else
            V_A(:,n) = -Vt(:,n);
        end
         
        
         %Ecuaciones diferenciales
        
        edprime(:,n+1) = (-edprime(:,n) + (mac_con(:,11)-mac_con(:,7)).*curqg(:,n))./(n*mac_con(:,14));
        eqprime(:,n+1) = (Efd(:,n) - eqprime(:,n) - (mac_con(:,6)-mac_con(:,7)).*curdg(:,n))./(n*mac_con(:,9));
        Efd(:,n+1) = (-Efd(:,n) + exc_con(:,4).*V_A(:,n))./(n*exc_con(:,5)); 
             
             
        if n == 1
            mac_spd(:,n+1) =(pmech(:,1) - pelect(:,n).*mac_pot(:,1) - mac_con(:,17).*(mac_spd(:,n)-ones(Ng,1)))./(n*2*mac_con(:,16));
            mac_ang(:,n+1) = basrad*(mac_spd(:,n)-ones(Ng,1))/n;
        else
            mac_spd(:,n+1) =( -pelect(:,n).*mac_pot(:,1) - mac_con(:,17).*mac_spd(:,n))./(n*2*mac_con(:,16));
            mac_ang(:,n+1) = basrad*(mac_spd(:,n))/n;
        end

        auxfi = 0;
        auxpsi = 0;
        for m = 1:n
            auxfi = auxfi + (n-m+1)/n*psi(:,m).*mac_ang(:,n-m+2);
            auxpsi = auxpsi + (n-m+1)/n*fi(:,m).*mac_ang(:,n-m+2);
        end
        fi(:,n+1) = auxfi;
        psi(:,n+1) = -auxpsi;
    end

    mac_ang = mac_ang*hvec';
    mac_spd = mac_spd*hvec';
    eqprime = eqprime*hvec';
    edprime = edprime*hvec';
    Efd = Efd*hvec';
    
    delta(1:Ng,k+1) = mac_ang;
    w(1:Ng,k+1) = mac_spd;
    Edp(1:Ng,k+1) = edprime;
    Eqp(1:Ng,k+1) = eqprime;
    Ef(1:Ng,k+1) = Efd;
    
    fi = sin(mac_ang);
    psi = cos(mac_ang);
end 
toc
t = 0:h:t_switch(4); % time
delta = (delta - delta(1,:))*180/pi;


% plot(t(1:end-1),w)
% xlabel('Tiempo (s)')
% ylabel('Velocidad (pu)')
% grid

