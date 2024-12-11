function [f] = pss(i,k,bus,flag)
% Syntax: [f] = pss(i,k,bus,flag)
%
% Purpose: power system stabilization model
%           
% Input: i - generator number
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - system dynamics computation
%
% Output: f - dummy variable 
%
% Files:
%
% See Also: exc_dc12, exc_st3

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1992

% pss variables
global  pss_con pss_pot 
global  pss1 pss2 pss3 dpss1 dpss2 dpss3

global  exc_sig mac_spd pelect basmva mac_int mac_con

if i ~= 0
  if pss_con(i,1) ~= 1 & pss_con(i,1) ~= 2 
    disp('PSS: inappropriate power system stablizer model')
    stop
  end
 else
  if pss_con(1,1) ~= 1 & pss_con(1,1) ~= 2 
    disp('PSS: inappropriate power system stabilizer model')
    stop
  end
end
f = 0;

[npss dum] =size(pss_con);
jay = sqrt(-1);
num_mac=length(mac_con(:,1));
if i == 0 & num_mac~=npss
  disp('pss data inconsistent with vector computation')
  stop
end
if flag == 0; % initialization
  if i ~= 0  % scalar computation
  n = mac_int(pss_con(i,2)); % machine number
    if pss_con(i,1) == 1
      pss1(i,1) = mac_spd(n,1);
     else
      pss1(i,1) = pelect(n,1)*basmva/mac_con(n,3);
    end
    pss2(i,1) = 0.;
    pss3(i,1) = 0.;
    exc_sig(i,1) = 0.;
    pss_pot(i,1) = pss_con(i,5)/pss_con(i,6);
    if pss_con(i,8) ~= 0
      pss_pot(i,2) = pss_con(i,7)/pss_con(i,8);
    end
  else
    %vector computation
    pss_pot=ones(num_mac,2);
    n=mac_int(pss_con(:,2));
    if pss_con(1,1)==1
     pss1(:,1) = zeros(num_mac,1);
    else
     pss1(:,1)=pelect(n,1)*basmva./mac_con(n,3);
    end
    pss2(:,1)=zeros(num_mac,1);
    pss3(:,1)=zeros(num_mac,1);
    exc_sig=zeros(num_mac,1);
    pss_pot(:,1)=pss_con(:,5)./pss_con(:,6);
    if min(abs(pss_con(:,8)))>0
       pss_pot(:,2)=pss_con(:,7)./pss_con(:,8);
    end
  end
end

if flag == 1 % network interface computation
 if i ~= 0 % scalar computation
    n = mac_int(pss_con(i,2)); % machine number
    if pss_con(i,1) == 1
      var1 = (mac_spd(i,k)-pss1(i,k))/pss_con(i,4);
     else
      n = mac_int(pss_con(i,2)); % machine number 
      var1 = (pelect(i,k)*basmva/mac_con(n,3)-pss1(i,k))/pss_con(i,4);
    end

    var2 = pss_pot(i,1)*pss_con(i,3)*var1 + pss2(i,k);

    if pss_con(i,8) == 0
      var3 = var2;
     else
      var3 = pss_pot(i,2)*var2 + pss3(i,k);
    end
    exc_sig(i,k) = min(pss_con(i,9),max(var3,-pss_con(i,9)));
  else
     % vector computation
    n = mac_int(pss_con(:,2)); % machine number vector
    if pss_con(1,1) == 1
      var1 = mac_spd(:,k)-pss1(:,k) - ones(num_mac,1);
     else
      var1 = (pelect(:,k)*basmva./mac_con(n,3)-pss1(:,k))./pss_con(:,4);
    end

    var2 = pss_pot(:,1).*(pss_con(:,3).*var1) + pss2(:,k);
   

    if min(abs( pss_con(:,8))) == 0
      var3 = var2;
     else
      var3 = pss_pot(:,2).*var2 + pss3(:,k);
     
    end
    exc_sig(:,k) = min(pss_con(:,9),max(var3,-pss_con(:,9)));    
  end
end

if flag == 2 % pss dynamics calculation
  if i ~= 0 % scalar computation
    n = mac_int(pss_con(i,2)); % machine number
    if pss_con(i,1) == 1
      var1 = (mac_spd(i,k)-pss1(i,k))/pss_con(i,4);
     else
      n = mac_int(pss_con(i,2)); % machine number 
      var1 = (pelect(i,k)*basmva./mac_con(n,3)-pss1(i,k))/pss_con(i,4);
    end
    dpss1(i,k) = var1;

    var2 = pss_pot(i,1)*pss_con(i,3)*var1 + pss2(i,k);
    dpss2(i,k) = ((1-pss_pot(i,1))*pss_con(i,3)*var1 - pss2(i,k))/pss_con(i,6);

    if pss_con(i,8) == 0
      var3 = var2;
      dpss3(i,k) = dpss2(i,k);
     else
      var3 = pss_pot(i,2)*var2 + pss3(i,k);
      dpss3(i,k) = ((1-pss_pot(i,2))*var2 - pss3(i,k))/pss_con(i,8);
    end
    exc_sig(i,k) = min(pss_con(i,9),max(var3,-pss_con(i,9)));
  else
     % vector computation
    n = mac_int(pss_con(:,2)); % machine number vector
    if pss_con(1,1) == 1
      var1 = mac_spd(:,k)-pss1(:,k) - ones(num_mac,1);
     else
      var1 = (pelect(:,k)*basmva./mac_con(n,3)-pss1(:,k))./pss_con(:,4);
    end
    dpss1(:,k) = var1./pss_con(:,4);

    var2 = pss_pot(:,1).*(pss_con(:,3).*var1) + pss2(:,k);
    dpss2(:,k) = ((ones(num_mac,1)-pss_pot(:,1)).*(pss_con(:,3).*var1 )- pss2(:,k))./pss_con(:,6);

    if min(abs( pss_con(:,8))) == 0
      var3 = var2;
      dpss3(:,k) = dpss2(:,k);
     else
      var3 = pss_pot(:,2).*var2 + pss3(:,k);
      dpss3(:,k) = ((ones(num_mac,1)-pss_pot(:,2)).*var2 - pss3(:,k))./pss_con(:,8);
    end
    exc_sig(:,k) = min(pss_con(:,9),max(var3,-pss_con(:,9)));    
  end
end
