%% Vehicle Platoon

close all; clc; clear;

%% Functions

kp = 1;
kd = 3.3;
d = 10;

N = 5; % Platoon size with Leader

for n = 1:N-1
    B(n*2,n:n+1) = [-1 1];
    
end

for n=0:N-1
    K(n+1,n*2+1:n*2+4) = [-kp -kd kp kd];
end

K  = K(:,3:end-2);
A = diag(mod((1:N*2-3),2)==1,1)*1;

%% LQR
% Q = eye(N*2);
% R = 50*eye(N);
% [K,S,E]=lqr(A,B,Q,R);

%% Attacque

K_a = K;

K_a(3,3:6) = [-1 2 1 -2];

%% simulation

% X0  = zeros(1,N*2);
% X0(1,end-1:end) = [0 50];

X0 = [-0 5 -0 5 -0 5 -0 5];

tfin = 100;
dt = 0.01;

clear X t
t=0:dt:tfin;

maxA = 13.4112; % 30 mph/s in m/s^2
minA = -13.4112; % -30 mph/s m/s^2

X(1,:) = X0';

a = 15;
f= 0.05;
w = 1.53;

for n=2:length(t)
%     E(n-1,:) = (R - C*X(n-1,:)')';
%     U(n-1,:) = ((K)*X(n-1,:)');
%     U(n-1,:) = ((K+K_a*a*sin(w*t(n)))*X(n-1,:)');
    U(n-1,:) = (K_a*X(n-1,:)');

    U(maxA<U)=maxA;
    U(minA>U)=minA;   

    Xd = A*X(n-1,:)' + B*U(n-1,:)';
    X(n,:) = X(n-1,:) + dt*Xd';
end

%% Plots
for n=1:N-1
figure(1); plot(t,X(:,(n)*2-1)+d); hold on;
n
end
title('Seperation');
xlabel('Time (s)');
ylabel('dP (m)');
legend('show')
grid on



for n=1:N-1
figure(2); plot(t,X(:,n*2).*2.23694); hold on
n
end
title('Error in Velocities');
xlabel('Time (s)');
ylabel('dV (mph)');
legend('show')
grid on


for n=1:N-1
figure(3); plot(t(1:end-1),U(:,n)./9.806); hold on
end
title('Accelerations');
xlabel('Time (s)');
ylabel('Accelerations (g)');
legend('show')
grid on


%% Robustness
% 
% Ks = -10:2:10;
% K_a = K;
% figure;
% for n=1:length(Ks)
%     K_a(3,6) = Ks(n);
%     K_a(3,4) = -Ks(n);
%     [b,a] = ss2tf(A,B(:,3),K_a(3,:),0);
%     RD_u3 = tf(b,a);
%     bode(RD_u3+1);hold on
% end
% 
% for n=1:length(Ks)
%     K_a(3,6) = Ks(n);
%     K_a(3,4) = -Ks(n)
%     [b,a] = ss2tf(A,B(:,3),K_a(3,:),0);
%     RD_u3 = tf(b,a);
%     figure;margin(RD_u3+1);legend(num2str(Ks(n)))
% end
