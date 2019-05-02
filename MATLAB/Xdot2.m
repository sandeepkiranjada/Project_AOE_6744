function [ Xd ] = Xdot2( t,X,A )
%XDOT Summary of this function goes here
%   Detailed explanation goes here


O = 2*pi*0.1;
% O = 1.2130;
O = 0.38;
a = 5;
W = @(t) a*sin(O.*t)+a; % Acceleration of the leader.
% W = @(t) t*0;
G = zeros(length(X),1);
G(end) = 1;

Xd = A * X + G*W(t);

end

