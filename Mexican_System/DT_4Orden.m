% m.file to simulate the 46 machine, 140 bus system
% using the Matlab Power System Toolbox

clear all; clear global; clc;

global  basmva basrad mac_con bus_int

d46_mexico_clas_4o

h = 1/60; % integration stepsize
p = 2;%Numero de coeficientes
hvec = h.^(0:p);

% simulation
t_switch(1) = 0;     % all time in second+s, start time
t_switch(2) = 1.0;  % time to apply fault
t_switch(3) = 1.0 + 3/60; % time to clear fault, 3 cycles
t_switch(4) = 20;   % end time

jay = sqrt(-1);
mac_con(:,20:21)=0;
mac_con(:,22:23)=1;


  tol = 1e-4;   % tolerance for convergence
  iter_max = 15; % maximum number of iterations
  vmin = 0.5;   % voltage minimum
  vmax = 1.5;   % voltage maximum
  acc = 1.0;   % acceleration factor

  [bus_sol,line_flw] = loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2);
  bus = bus_sol;
  
%mac_con(:,17) = mac_con(:,16);
Ng = length(mac_con(:,1));
xpd = mac_con(1:Ng,7);
H = mac_con(:,16);
D = mac_con(:,17);


basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base


disp('Performing simulation.')



k_switch(1) = round((t_switch(2) - t_switch(1))/h) + 1;
k_switch(2) = round((t_switch(3) - t_switch(1))/h) + 1;
k_switch(3) = round((t_switch(4) - t_switch(1))/h) + 1;

% step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
[Y_red,V_rec] = red_ybus(bus,line); %pre-fault admittance matrix

% create bus matrix with load increased on node 9
bus_f = bus;
bus_f(5,9) = -150.0;
[Y_red_f,V_rec_f] = red_ybus(bus_f,line); %fault-on admittance matrix


% line_pf = line;
% nline = length(line(:,1));
% line_pf(nline+1,:)=[32	33	-0.0007	-0.0089	-0.1342	0	0];% line with negative impedance
% [Y_red_pf,V_rec_pf] = red_ybus(bus,line_pf); %post-fault 

%step 2: initialization
%modelo = [4*ones(1,8) 2 4 4 2 2 2 2 4*ones(1,5) 2 2 2 2 4 4 4 2 2 2 4*ones(1,10) 2 4 2 4 4 4];

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
    mac_spd(i,1) = 0; % machine speed at steady state
                            
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
    Efd(i,1) = vex(i,1);
        
    V_A(i,1) = Efd(i,1)/exc_con(i,4); % laglead

    exc_pot(i,4) = exc_con(i,7)/exc_con(i,5);
    exc_pot(i,3) = eterm(i,1) + V_A(i,1); % reference voltage
    exc_pot(i,5) = 1;
end

delta = mac_ang;
w = mac_spd;
Edp = edprime;
Eqp = eqprime;
Ef = Efd;
fi = sin(mac_ang);
psi = cos(mac_ang);


  % generator dynamics calculation
tic
for k = 1:k_switch(3)-1

    for n = 1:p 
        for i = 1:Ng
            vx = 0;
            vy = 0;
            for m = 1:n
                vx = vx + fi(i,m)*edprime(i,n-m+1) + psi(i,m)*eqprime(i,n-m+1);
                vy = vy - psi(i,m)*edprime(i,n-m+1) + fi(i,m)*eqprime(i,n-m+1);
            end
            vxy(i,n) = vx + vy*1i;
        end
        % network-machine interface

        if k >= k_switch(2) 
            Ixy(:,n) = Y_red*vxy(:,n); % network solution corrientes de generadores POSTFALLA
            bus_v(:,n) = V_rec*vxy(:,n); % bus voltage reconstruction
        elseif k >= k_switch(1)
            Ixy(:,n) = Y_red_f*vxy(:,n);  % network solution FALLA
            bus_v(:,n) = V_rec_f*vxy(:,n); % bus voltage reconstruction
        else
            Ixy(:,n) = Y_red*vxy(:,n); % network solution PREFALLA
            bus_v(:,n) = V_rec*vxy(:,n); % bus voltage reconstruction
        end
            
        Ix(:,n) = real(Ixy(:,n)); 
        Iy(:,n) = imag(Ixy(:,n));
         
        auxcurd = 0;
        auxcurq = 0;
        for m = 2:n
            auxcurd = auxcurd + fi(:,m).*curd(:,n-m+1) + psi(:,m).*curq(:,n-m+1);
            auxcurq = auxcurq - psi(:,m).*curd(:,n-m+1) + fi(:,m).*curq(:,n-m+1);
        end
        curd(:,n) = fi(:,1).*(Ix(:,n)-auxcurd) - psi(:,1).*(Iy(:,n)-auxcurq);
        curq(:,n) = psi(:,1).*(Ix(:,n)-auxcurd) + fi(:,1).*(Iy(:,n)-auxcurq);
            
        curdg(:,n) = curd(:,n).*mac_pot(:,1);
        curqg(:,n) = curq(:,n).*mac_pot(:,1);

        ed(:,n) = edprime(:,n) + mac_con(:,7).*curqg(:,n);
        eq(:,n) = eqprime(:,n) - mac_con(:,7).*curdg(:,n);
        

            auxpelect = 0;
            auxed = 0;
            auxeq = 0;
            for m = 1:n
                auxpelect = auxpelect + ed(:,m).*curdg(:,n-m+1) + eq(:,m).*curqg(:,n-m+1);
                auxed = auxed + ed(:,m).*ed(:,n-m+1);
                auxeq = auxeq + eq(:,m).*eq(:,n-m+1);
            end
            pelect(:,n) = auxpelect;
            x(:,n) = auxed + auxeq;
           
            if n == 1 
                eterm(:,1) = sqrt(x(:,1));
            else
                auxy = 0;
                for m = 2:n-1
                    auxy = auxy + eterm(:,m).*eterm(:,n-m+1);
                end
                eterm(:,n) = x(:,n)./(2*eterm(:,1)) - auxy./(2*eterm(:,1));
            end
              
            if n == 1
                V_A(:,n) = exc_pot(:,3) - eterm(:,n);
            else
                V_A(:,n) = -eterm(:,n);
            end
         

%=================================================
         %Ecuaciones diferenciales
            mac_ang(:,n+1) = basrad./n*(mac_spd(:,n));
            
            if n == 1
                mac_spd(:,n+1) = (pmech - pelect(:,n) - D.*(mac_spd(:,n)))./(2*H*n); 
            else
                mac_spd(:,n+1) = ( -pelect(:,n)*mac_pot(i,1) - D.*(mac_spd(:,n)))./(2*H*n); 
            end
            
            edprime(:,n+1) = (-edprime(:,n) + (mac_con(:,11)-mac_con(:,7)).*curqg(:,n))./(n*mac_con(:,14));
            eqprime(:,n+1) = (Efd(:,n) - eqprime(:,n) - (mac_con(:,6)-mac_con(:,7)).*curdg(:,n))./(n*mac_con(:,9));
            Efd(:,n+1) = (-Efd(:,n) + exc_con(:,4).*V_A(:,n))./(n*exc_con(:,5)); 
 
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

    delta(:,k+1) = mac_ang;
    w(:,k+1) = mac_spd;
    Edp(:,k+1) = edprime;
    Eqp(:,k+1) = eqprime;
    Ef(:,k+1) = Efd;
    
    fi = sin(mac_ang);
    psi = cos(mac_ang);
end 
toc
t = 0:h:t_switch(4); % time

delta(1:Ng,:) = (delta(1:Ng,:) - delta(1,:))*180/pi;
w = w+1;

%figure(1)
% plot(t,delta(2,:))
% xlabel('Tiempo (s)')
% ylabel('E_d^{\prime} (pu)')
% grid

% figure(2)
% plot(t,Eqp)
% xlabel('Tiempo (s)')
% ylabel('E_q^{\prime} (pu)')
% grid
% 
% figure(3)
% plot(t,Ef)
% xlabel('Tiempo (s)')
% ylabel('E_{fd} (pu)')
% grid


