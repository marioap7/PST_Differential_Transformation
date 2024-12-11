
function [f] = mac_tra(i,k,bus,flag)
% Syntax: [f] = mac_tra(i,k,bus,flag)
%
% Purpose: voltage-behind-transient-reactance generator
%          model, with vectorized computation option
%
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation 
%
% 
pst_var
global  basmva basrad mach_ref
global  psi_re psi_im cur_re cur_im bus_int Ng

% synchronous machine variables
global  mac_con mac_pot busnum mac_ang mac_spd eqprime edprime 
global  curd curq curdg curqg vex eterm theta ed eq 
global  pmech pelect qelect dmac_ang dmac_spd deqprime dedprime 



jay = sqrt(-1);
if flag == 0 % initialization
    % vectorized computation
      [nmach,dum] = size(mac_con);
      if dum < 23                    % set power fraction
        mac_con = [mac_con ones(nmach,2)];  % to unity
      end
      if mac_con(i,14) == 0
            mac_con(i,14) = 999.0;
      end
      busnum = bus_int(mac_con(i,2)); % bus number 
      mac_pot(i,1) = basmva/mac_con(i,3); % scaled MVA base                  
      mac_pot(i,2) = 1; % base kv
      
      % extract bus information
      eterm(i,1) = bus(busnum,2);  % terminal bus voltage
      theta(i,1) = bus(busnum,3)*pi/180;  
                          % terminal bus angle in radians
      qelect(i,1) = bus(busnum,5)*mac_con(i,23);                            
      pelect(i,1) = bus(busnum,4); %.*mac_con(:,22);  
                        % electrical output power, active
                        % electrical output power, reactive
                        
      curr = sqrt(pelect(i,1)^2 + qelect(i,1)^2) ...
            /eterm(i,1)*mac_pot(i,1);  % current magnitude
      phi = atan2(qelect(i,1),pelect(i,1)); 
                                        % power factor angle
      v = eterm(i,1)*(cos(theta(i,1))+jay*sin(theta(i,1))); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
      curr = curr*(cos(theta(i,1)-phi) + jay*sin(theta(i,1)-phi)); 
 
      eprime = v + jay*mac_con(i,7)*curr; 
      ei = v + jay*mac_con(i,11)*curr;
      
      mac_ang(i,1) = atan2(imag(ei),real(ei)); % machine angle (delta)
      mac_spd(i,1) = 1; 
                            % machine speed at steady state
      rot = sin(mac_ang(i,1))+jay*cos(mac_ang(i,1)); 
                          % system reference frame rotation
      eprime = eprime*rot;
      edprime(i,1) = real(eprime); 
      eqprime(i,1) = imag(eprime); 
      curr = curr*rot;
      
      curdg(i,1) = real(curr); 
      curqg(i,1) = imag(curr);
      curd(i,1) = real(curr)/mac_pot(i,1); 
      curq(i,1) = imag(curr)/mac_pot(i,1);
      v = v*rot;
      
      ed(i,1) = real(v); 
      eq(i,1) = imag(v);
 
      %E_Isat = eqprime(i,1); % select higher voltage
      vex(i,1) = eqprime(i,1) + (mac_con(i,6) - mac_con(i,7))*curdg(i,1);
      pmech(i,1) = pelect(i,1)*mac_pot(i,1); % set input
end

if flag == 1 % network interface computation 
    % vectorized computation
      [nmach,dum] = size(mac_con);
      mac_ang(:,k) = mac_ang(:,k);
                     % wrt machine reference
      psi_re(:,k) = sin(mac_ang(:,k)).*edprime(:,k) + ...
                    cos(mac_ang(:,k)).*eqprime(:,k); % real part of psi
      psi_im(:,k) = -cos(mac_ang(:,k)).*edprime(:,k) + ...
                     sin(mac_ang(:,k)).*eqprime(:,k); % imag part of psi
end
if flag == 2 % generator dynamics calculation

    % vectorized computation
      [nmach,dum] = size(mac_con);
      curd(:,k) = sin(mac_ang(:,k)).*cur_re(:,k) - ...
                  cos(mac_ang(:,k)).*cur_im(:,k); % d-axis current
      curq(:,k) = cos(mac_ang(:,k)).*cur_re(:,k) + ...
                  sin(mac_ang(:,k)).*cur_im(:,k); % q-axis current
      curdg(:,k) = curd(:,k).*mac_pot(:,1);
      curqg(:,k) = curq(:,k).*mac_pot(:,1);

      dedprime(:,k) = (-edprime(:,k) + (mac_con(:,11)-...
                  mac_con(:,7)).*curqg(:,k))./mac_con(:,14);
      deqprime(:,k) = (vex(:,k) - eqprime(:,k) - (mac_con(:,6)-...
                  mac_con(:,7)).*curdg(:,k))./mac_con(:,9);
      
      ed(:,k) = edprime(:,k) + mac_con(:,7).*curqg(:,k);
      eq(:,k) = eqprime(:,k) - mac_con(:,7).*curdg(:,k);
      eterm(:,k) = sqrt(ed(:,k).^2  + eq(:,k).^2);

      pelect(:,k) = eq(:,k).*curq(:,k) + ed(:,k).*curd(:,k);
      qelect(:,k) = eq(:,k).*curd(:,k) - ed(:,k).*curq(:,k);
      
      dmac_ang(:,k) = basrad*(mac_spd(:,k)- ones(1,Ng)');
      dmac_spd(:,k) =(pmech(:,k) - pelect(:,k).*mac_pot(:,1) ...
                      -mac_con(:,17).*(mac_spd(:,k) - ones(1,Ng)'))./(2*mac_con(:,16));
end
