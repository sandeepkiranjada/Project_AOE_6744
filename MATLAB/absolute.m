%% Vehicle Platoon

close all; clc; clear;

%% Functions

kpi = 1;
kdi = 12;
d = 10;

N = 5; % Platoon size with Leader

R = or(mod((3:N*5-2),5)==1,mod((3:N*5-2),5)==3)' .* d; % Reference Signal
R(end) = 22.352;

W = @(t) t*0; % Acceleration of the leader.
G = zeros(2*N,1);
G(end) = 1;

A = diag(mod((1:N*2-1),2)==1,1)*1;
B = zeros(2*N,N);
C = zeros(N*5-4,2*N);
for n=2:2:2*N
    B(n,n/2) = 1;
end

Ci = [-1  0  1  0  0  0; ...
       0 -1  0  1  0  0; ...
       0  0 -1  0  1  0; ...
       0  0  0 -1  0  1; ...
       0  0  0  1  0  0];
   
C1 = [ -1  0  1  0; ...
        0 -1  0  1; ...
        0  1  0  0];

CN = [ -1  0  1  0; ...
        0 -1  0  1; ...
        0  0  0  1];

C(1:3,1:4) = C1;
C(end-2:end,end-3:end) = CN;

K = zeros(N,length(R));
Ki = [kpi kdi -kpi -kdi 0];
K1 = [-kpi -kdi 0];
KN = [kpi kdi kdi];

K(1,1:3) = K1;
K(end,end-2:end) = KN;

for n=0:N-3
    C(((n)*5+4):((n)*5+8),((n)*2 + 1):((n)*2 + 6)) = Ci;
    
    K(n+2,n*5+4:n*5+8) = Ki;
end

% simulate instability
K_a = K;
% %kps
% K_a(4,14) = -3;
% K_a(4,16) = 3;
%kds
K_a(3,10) = -1.5;
K_a(3,12) = -K_a(3,10);
% K_a(4,14) = 1.2;
% K_a(4,16) = -K_a(4,14);

% attacker modifying C
% C(7,:) = -0.7*C(7,:);
% C(15,:) = -0.7*C(15,:);

% B(end,end) = 0;

A_tilda = A-B*K*C;
%A_tilda = A-B*K_a*C;
% U_tilda = B*K*R+G*W(1);

U_tilda = B*K*R;
%U_tilda = B*K_a*R;

% Eigen_A_tilda = eig(A_tilda);

%% simulation

X0  = [0 0 10 0 20 0 30 0 40 0]';
% X0  = zeros(N*2,1);
% X0(end) = 22.352; % 50 mph in m/s


% Xdot = @(t,X) (A_tilda * X + U_tilda);

% [t,X] = ode45(@(t,X) Xdot(t,X,A_tilda,U_tilda),[0 20],X0);

%% Simulation Euler-Cauchy
tf = 50;
dt = 0.01;

clear X t
t=0:dt:tf;

maxA = 13.4112; % 30 mph/s in m/s^2
minA = -13.4112; % -30 mph/s m/s^2

X(1,:) = X0';

for n=2:length(t)
    E(n-1,:) = (R - C*X(n-1,:)')';
    
    U(n-1,:) = (K*E(n-1,:)');
    U(n-1,2) = U(n-1,2) - kdi*10*sin(-0*pi/2+0.1*t(n));
    U(n-1,4) = U(n-1,4) + kdi*10*sin(0.1*t(n));
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

for n=1:N
figure(1); plot(t,X(:,n*2-1)); hold on;
end
title('Positions');
xlabel('Time (s)');
ylabel('Positions (m)');
legend('show')
grid on



for n=1:N
figure(2); plot(t,X(:,n*2).*2.23694); hold on
end
title('Velocities');
xlabel('Time (s)');
ylabel('Velocities (mph)');
legend('show')
grid on



for n=1:N-1
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
figure(5); plot(t(1:end-1),U(:,n)./9.806); hold on
end
title('Accelerations');
xlabel('Time (s)');
ylabel('Accelerations (g)');
legend('show')
grid on


