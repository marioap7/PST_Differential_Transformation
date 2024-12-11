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

global  basmva basrad mach_ref
global  psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot busnum mac_ang mac_spd eqprime edprime 
global  curd curq curdg curqg vex eterm theta 
global  pmech pelect qelect dmac_ang dmac_spd deqprime dedprime 



jay = sqrt(-1);
if flag == 0 % initialization
    % vectorized computation
      [nmach,dum] = size(mac_con);
      if dum < 23                    % set power fraction
        mac_con = [mac_con ones(nmach,2)];  % to unity
      end
      busnum = bus_int(mac_con(:,2)); % bus number 
      mac_pot(:,1) = basmva*ones(nmach,1)./mac_con(:,3); 
                          % scaled MVA base
      mac_pot(:,2) = ones(nmach,1); % base kv
      % extract bus information
      eterm(:,1) = bus(busnum,2);  % terminal bus voltage
      theta(:,1) = bus(busnum,3)*pi/180;  
                          % terminal bus angle in radians
      qelect(:,1) = bus(busnum,5);%.*mac_con(:,23);                            
      pelect = bus(busnum,4); %.*mac_con(:,22);  
                        % electrical output power, active
                        % electrical output power, reactive
      curr = sqrt(pelect(:,1).^2+qelect(:,1).^2) ...
            ./eterm(:,1).*mac_pot(:,1);  % current magnitude
      phi = atan2(qelect(:,1),pelect(:,1)); 
                                        % power factor angle
      v = eterm(:,1).*(cos(theta(:,1))+jay*sin(theta(:,1))); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
      curr = curr.*(cos(theta(:,1)-phi) + jay*sin(theta(:,1)-phi)); 
 
      eprime = v + jay*mac_con(:,7).*curr; 
      ei = v + jay*mac_con(:,11).*curr;
      
      mac_ang = atan2(imag(ei),real(ei)); % machine angle (delta)
      mac_spd = ones(nmach,1); 
                            % machine speed at steady state
      rot = sin(mac_ang(:,1))+jay*cos(mac_ang(:,1)); 
                          % system reference frame rotation
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
end

if flag == 1 % network interface computation 
    % vectorized computation
      [nmach,dum] = size(mac_con);
      mac_ang(:,k) = mac_ang(:,k)-mach_ref(k)*ones(nmach,1);
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

      for i = 1:nmach
         dedprime(i,k) = (-edprime(i,k) + (mac_con(i,11)-...
                  mac_con(i,7)).*curqg(i,k))./mac_con(i,14);
         deqprime(i,k) = (vex(i,k) - eqprime(i,k) - (mac_con(i,6)-...
                  mac_con(i,7)).*curdg(i,k))./mac_con(i,9);
      end
      ed(:,k) = edprime(:,k) + mac_con(:,7).*curqg(:,k);
      eq(:,k) = eqprime(:,k) - mac_con(:,7).*curdg(:,k);
      eterm(:,k) = sqrt(ed(:,k).^2+eq(:,k).^2);

      pelect(:,k) = eq(:,k).*curq(:,k) + ed(:,k).*curd(:,k);
      qelect(:,k) = eq(:,k).*curd(:,k) - ed(:,k).*curq(:,k);
      dmac_ang(:,k) = basrad*(mac_spd(:,k)-ones(nmach,1));
      dmac_spd(:,k) =(pmech(:,k) - pelect(:,k).*mac_pot(:,1) ...
                      -mac_con(:,17).*(mac_spd(:,k) - ones(nmach,1)))./(2*mac_con(:,16));
end
