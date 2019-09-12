%% Vehicle Platoon

close all; clc; clear;

%% Functions

d = 10;

N = 7; % Platoon size with Leader

B = zeros(2*N,N);

for n=2:2:2*N
    B(n,n/2) = 1;
end

A = diag(mod((1:N*2-1),2)==1,1)*1;


%% LQR
Q = diag([10 1]);
R = 1;


Q = Q./R;
R = 1;
[Klq,S,E]=lqr([0 1;0 0],[0 1]',Q,R);

%% Gains

kp = Klq(1);
kd = Klq(2);

for n=0:N-2
    Kb(n+1,n*2+1:n*2+2) = [kp kd];
end

% K  = K(:,3:end-2);
Kb(end-1,end-1) = 0;
Kb(end,:) = Kb(end,:).*0;
Kb(end,end) = kd;

TT = triu(-1*ones(N-1));
TP = [zeros(2*N-2,2),diag(ones(2*N-2,1))]-[diag(ones(2*N-2,1)),zeros(2*N-2,2)];
K = TT*Kb;
K = [K;zeros(1,2*N-2)];
Sep = [mod((1:N*2-2),2)==1]'*d;
Sep(end-2) = 0;
% K = K;

%% Attacque

K_a = Kb;

K_a(3,3:6) = [-1 2 1 -2];

%% simulation

% X0  = zeros(1,N*2);
% X0(1,end-1:end) = [0 50];

X0  = [0 0 10 0 20 0 30 0 40 0 50 0 50 22.3]';

tfin = 50;
dt = 0.01;

clear X t
t=0:dt:tfin;

maxA = 13.4112; % 30 mph/s in m/s^2
minA = -13.4112; % -30 mph/s m/s^2

X(1,:) = X0';

a = 15;
f= 0.05;
w = 1.53;

%% Simulation Euler-Cauchy
tf = 50;
dt = 0.01;

clear X t
t=0:dt:tf;

maxA = 13.4112; % 30 mph/s in m/s^2
minA = -13.4112; % -30 mph/s m/s^2

X(1,:) = X0';

for n=2:length(t)
%     E(n-1,:) = (R - C*X(n-1,:)')';
    
    U(n-1,:) = (-K*(TP*X(n-1,:)'-Sep)) ;
%     U(n-1,2) = U(n-1,2) - kd*10*sin(-0*pi/2+0.1*t(n));
%     U(n-1,4) = U(n-1,4) + kd*10*sin(0.1*t(n));
%     U(n-1,:) = (K_a*E(n-1,:)');

    U(maxA<U)=maxA;
    U(minA>U)=minA;   

    Xd = A*X(n-1,:)' + B*U(n-1,:)';
    X(n,:) = X(n-1,:) + dt*Xd';
end



%% Plots
% U = -K*C*X'+K*R;
%U = -K_a*C*X'+K_a*R;

maxA = 13.4112; % 30 mph/s in m/s^2
minA = -13.4112; % -30 mph/s m/s^2

% U(maxA<U)=maxA;
% U(-maxA>U)=maxA;

for n=1:N-1
figure(1); plot(t,X(:,n*2-1)); hold on;
end
title('Positions');
xlabel('Time (s)');
ylabel('Positions (m)');
legend('show')
grid on



for n=1:N-1
figure(2); plot(t,X(:,n*2).*2.23694); hold on
end
title('Velocities');
xlabel('Time (s)');
ylabel('Velocities (mph)');
legend('show')
grid on



for n=1:N-2
figure(3); plot(t,X(:,(n+1)*2-1)-X(:,(n)*2-1)); hold on;
end
title('Seperation');
xlabel('Time (s)');
ylabel('dP (m)');
legend('show')
grid on



for n=1:N-1
figure(4); plot(t,X(:,(n+1)*2).*2.23694-X(:,n*2).*2.23694); hold on
end
title('Error in Velocities');
xlabel('Time (s)');
ylabel('dV (mph)');
legend('show')
grid on


for n=1:N-1
figure(5); plot(t(2:end),U(:,n)./9.806); hold on
end
title('Accelerations');
xlabel('Time (s)');
ylabel('Accelerations (g)');
legend('show')
grid on


%% Energy Calculation


