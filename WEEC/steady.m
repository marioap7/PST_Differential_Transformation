% disp('Calculo del estado permanente de Generadores')
% disp('Se suponen En Memoria los siguientes datos :')
% disp(' Ng= Numero de Generadores')
% disp(' P = Vector de Potencia Real entregada por Cada Maquina  [pu]')
% disp(' Q= Vector de Potencia Reactiva entregada por Cada Maquina  [pu]')
% disp(' V = Vector de Magnitudes de Voltaje Terminal  [pu]')
% disp(' Theta = Vector de Angulos del Voltaje Terminal  [Grados]')
% disp(' Ra = Vector de resistencias de armadura  [pu]')
% disp(' Xq = Vector de Reactancias en eje-q  [pu]')
% disp(' Xd = Vector de Reactancias en eje-d  [pu]')
% disp(' Xpq = Vector de Reactancias Transitorias en eje-q  [pu]')
% disp(' Xpd = Vector de Reactancias Transitorias en eje-d  [pu]')
% disp(' X2q = Vector de Reactancias Subtransitorias en eje-q  [pu]')
% disp(' X2d = Vector de Reactancias Subtransitorias en eje-d  [pu]')

% Ng = 1; P = 1.0; Q = 0.47975;
%  D = 0; Xd = 0.05; X2d = 0.012; Tpq0 = 0.75; V =1.0298; Theta = 3.8280; % ORIGINAL
%    Ra = 0.0; Xpq = 0.010; X2q = 0.012; H = 3.5; Vinf = 1.0;
%       Xq = 0.05; Xpd = 0.015; Tpd0 = 4.5; Re = 0.0; Xe = 0.06875;
%          Ka = 75; Ta=0.05; Ke=-0.017; Aex=0; Bex=0.1; Mod=0;

%Ng = 1; P = 2.25; Q = 1.7219;
% D = 0; Xd = 0.05; X2d = 0.012; Tpq0 = 0.75; V =1.2292 ; Theta = 24.8972;  %Xpd = Xpq
%   Ra = 0.0; Xpq = 0.0125; X2q = 0.012; H = 3.5; Vinf = 1.0;
%      Xq = 0.05; Xpd = 0.0125; Tpd0 = 4.5; Re = 0.0; Xe = 0.23;
%         Ka = 75; Ta=0.05; Ke=-0.017; Aex=0; Bex=0.1; Mod=0;


i=sqrt(-1);
for k=1:Ng,
Theta(k) = Theta(k)*pi/180.0;
Vt(k)=V(k)*(cos(Theta(k))+i*sin(Theta(k)));
%disp('Corriente  [pu]  :')
Cor=(P(k)-i*Q(k))/conj(Vt(k));
Ix=real(Cor);
Iy=imag(Cor);
Vx=real(Vt(k));
Vy=imag(Vt(k));
VEq(k,1)=Vt(k)+(Ra(k)+i*Xq(k))*Cor;
Eprime(k,1)=Vt(k)+(Ra(k)+i*Xq(k))*Cor;

Eprime2(k,1) = Eprime(k,1);
%disp('Voltaje  Eq [pu]  y  Angulo Delta [Grados]')
Eq(k)=abs(VEq(k,1));
Del(k)=atan2(imag(VEq(k,1)),real(VEq(k,1)));
Delta=Del(k)*180.0/pi;

rot = sin(Del(k))+i*cos(Del(k)); 
                          % system reference frame rotation
      Eprime2(k,1) = Eprime2(k,1)*rot;
%disp(' Epd  y  Epq  [pu]')      
      Epd(1,k) = real(Eprime2(k,1)); 
      Epq(1,k) = imag(Eprime2(k,1));

Id(k)=Ix*sin(Del(k))-Iy*cos(Del(k));
Iq(k)=Ix*cos(Del(k))+Iy*sin(Del(k));
Vd(k)=Vx*sin(Del(k))-Vy*cos(Del(k));
Vq(k)=Vx*cos(Del(k))+Vy*sin(Del(k));
%disp(' Efd [pu]  :')
Efd(k)=Eq(k)-(Xq(k)-Xd(k))*Id(k);
%pause

%disp('Fluko-d  y  Fluko-q  [pu]')
Fd(k)=Vq(k)+Ra(k)*Iq(k);
Fq(k)=-Vd(k)-Ra(k)*Id(k);
%pause
%disp('Par Electrico  [pu]')
Pelec(k)=Fd(k)*Iq(k)-Fq(k)*Id(k);
%disp(' E2d  y  E2q  [pu]')
E2d(k)=(Xq(k)-X2q(k))*Iq(k);
E2q(k)=Efd(k)-(Xd(k)-X2d(k))*Id(k);
end

Epq
Del
Pelec