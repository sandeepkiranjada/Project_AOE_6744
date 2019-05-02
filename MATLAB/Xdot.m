function [ Xd ] = Xdot( t,X,A,U )
%XDOT Summary of this function goes here
%   Detailed explanation goes here

O = 2*pi*0.1;
% O = 1.2130;
O = 1;
a = 5;
% W = @(t) a*sin(O.*t)+a; % Acceleration of the leader.
W = @(t) t*0;
G = zeros(length(U),1);
G(end) = 1;


maxV = 44.704; % 100 mph in m/s
minV = -6.7056; % -15 mph in m/s

maxA = 13.4112; % 30 mph/s in m/s^2
minA = -13.4112; % -30 mph/s m/s^2

N = length(X);

Xd_max  = ((mod(1:N,2)==0)*maxV + (mod(1:N,2)==1)*maxA)';
Xd_min  = ((mod(1:N,2)==0)*minV + (mod(1:N,2)==1)*minA)';


Xd = A * X + U + G*W(t);

% Xd(Xd_max<Xd)=Xd_max(Xd_max<Xd);
% 
% Xd(Xd_min>Xd)=Xd_min(Xd_min>Xd);

end

