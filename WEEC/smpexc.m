function [f] = smpexc(i,k,bus,flag,h)
% Syntax: [f] = smpexc(i,k,bus,flag,h)
%
% Purpose: simple excitation system, (exc_con(i,1)=0)
%            with vectorized computation option
%           
% Input: i - generator number
%            0, vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation
%        h - integration stepsize, for anti-windup reset

% synchronous machine variables
global  vex eterm mac_int

% excitation system variables
global  exc_con exc_pot Efd V_R V_A V_As dEfd dV_R dV_As exc_sig



[nexc dum] =size(exc_con);
jay = sqrt(-1);

if flag == 0 % initialization
    n = mac_int(exc_con(:,2)); % machine number
    Efd(:,1) = vex(n,1);
    V_A(:,1) = Efd(:,1)./exc_con(:,4); % laglead
    V_As(:,1) = V_A(:,1); % leadlag state variable 
    V_R(:,1)=zeros(nexc,1);
    exc_pot(:,4) = exc_con(:,7)./exc_con(:,5);
    err = V_A(:,1); % summing junction error
    exc_pot(:,5)=ones(nexc,1);
    exc_pot(:,3) = eterm(n,1)+err; % reference voltage
end


if flag == 1 % network interface computation
   n = mac_int(exc_con(:,2)); % machine number
   vex(n,k) = Efd(:,k); % set field voltage for machines
end

if flag == 2 % exciter dynamics calculation

  n = mac_int(exc_con(:,2)); % machine number
  [nexc dum] =size(exc_con);
  
  err = exc_sig(n,k) + exc_pot(:,3) - eterm(n,k);
  dV_R(:,k)=zeros(nexc,1);
  V_R(:,k)=eterm(n,k);
  err=exc_sig(n,k)+exc_pot(:,3)-V_R(:,k);
  dV_As(:,k) = zeros(nexc,1);
  V_As(:,k) = err;
  V_A(:,k) = err;
  dEfd(:,k) = (-Efd(:,k)+exc_con(:,4).*V_A(:,k))./exc_con(:,5);
end
