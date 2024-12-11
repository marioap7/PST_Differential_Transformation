function [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi)
% Syntax:   [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi) 
%
% Purpose:  Compute line flows. Inputs can be vectors.
%
% Input:    V1        - from bus complex voltage
%           V2        - to bus complex voltage
%           R         - line resistance
%           X         - line reactance
%           B         - line charging
%           tap       - tap ratio
%           phi       - phase shifter angle in degrees
% Output:   S1        - complex power injection at from bus
%           S2        - complex power injection at to bus
%
% See also:  
%
% Algorithm: 
%
% Calls:     
%

% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   1.0
% Authors:   Joe H. Chow
% Date:      March 1992
%
% ***********************************************************
% jay = sqrt(-1);
% nline = length(V1);
% for i = 1:nline
%   if tap(i) == 0
%     tap(i) = 1;
%   end
% end
% tps = tap.*exp(jay*phi*pi/180);
% 
% z = R + jay*X;
% y = ones(nline,1)./z;
% 
% S1 = V1.*conj((V1 - tps.*V2).*y ...
%        + V1.*(jay*B/2))./(tps.*conj(tps));
% 
% S2 = V2.*conj((V2 - V1./tps).*y ...
%        + V2.*(jay*B/2));
jay = sqrt(-1);
[nline,dummy] = size(V1);
for i = 1:nline
  if tap(i) == 0
    tap(i) = 1;
  end
end
tps = tap.*exp(jay*phi*pi/180);
tpsi = diag(ones(nline,1)./tps);
tps = diag(tps);
z = R + jay*X;
y = diag(ones(nline,1)./z);
chg = diag(jay*B/2);
cur1 = tps*(y*(tpsi*V1-V2) + chg*V1);
cur2 = y*(V2 - tpsi*V1) + chg*V2;
S1 = V1.*conj(cur1);
S2 = V2.*conj(cur2);