function [y1, y2] = fun_Generador_RT(t0,y0,v0)

    global k k_switch Y_red Y_red_f mac_pot
    jay = sqrt(-1);
    pst_var % set up global variable
    
    edprime = y0(:,1);
    eqprime = y0(:,2);
    mac_ang = y0(:,3);
    mac_spd = y0(:,4);
    efd = y0(:,5);
    
    % network interface computation 

  psi_re = sin(mac_ang).*edprime + cos(mac_ang).*eqprime; % real part of psi
  psi_im = -cos(mac_ang).*edprime + sin(mac_ang).*eqprime; % imag part of psi


  psi = psi_re + jay*psi_im; % voltajes internos en marco Re-Im

  if k >= k_switch(2)
     cur = Y_red*psi; % network solution corrientes de generadores POSTFALLA
      %bus_v(:,k) = V_rec*psi; % bus voltage reconstruction      
  elseif k >= k_switch(1)
      cur = Y_red_f*psi;  % network solution FALLA
      %bus_v(:,k) = V_rec_f*psi; % bus voltage reconstruction      
  else
     cur = Y_red*psi; % network solution PREFALLA
      %bus_v(:,k) = V_rec*psi; % bus voltage reconstruction      
  end
  
  cur_re = real(cur); 
  cur_im = imag(cur);
  
  % step 3b: compute dynamics and integrate
  %pmech(:,k) = pmech(:,1); % constant mechanical input power 
  
  % generator dynamics calculation
      curd = sin(mac_ang).*cur_re - cos(mac_ang).*cur_im; % d-axis current
      curq = cos(mac_ang).*cur_re + sin(mac_ang).*cur_im; % q-axis current
      
      curdg = curd.*mac_pot(:,1);
      curqg = curq.*mac_pot(:,1);
      
      ed = edprime + mac_con(:,7).*curqg;
      eq = eqprime - mac_con(:,7).*curdg;
      eterm = sqrt(ed.^2 + eq.^2);

      pelect = eq.*curq + ed.*curd;
      %qelect(:,k) = eq.*curd - ed.*curq;
      
      V_A = exc_pot(:,3) - eterm;
    
    %Solucion de las ecuaciones diferenciales
    
        mac_ang = basrad*(mac_spd - 1);
        mac_spd =(pmech - pelect.*mac_pot(:,1) - mac_con(:,17).*(mac_spd - 1))./(2*mac_con(:,16));
        edprime = (-edprime + (mac_con(:,11) - mac_con(:,7)).*curqg)./mac_con(:,14);
        eqprime = (efd - eqprime - (mac_con(:,6) - mac_con(:,7)).*curdg)./mac_con(:,9);
                %for i = 1:Ng
        efd = (-efd + exc_con(:,4).*V_A)./exc_con(:,5);
        %end
    y1 = [edprime eqprime mac_ang mac_spd efd];
    y2 = [pelect curd curq ed eq V_A];   
end
