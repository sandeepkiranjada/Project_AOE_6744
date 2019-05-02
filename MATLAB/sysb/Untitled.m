clear; clc; close all;

N  = 5;
U = sym('U',[N 1]);

Temp = eye(N)-diag(ones(1,N-1),1);

T1 = Temp(1:end-1,:);

X = sym('x',[N 1]);
V = sym('v',[N 1]);

Dx = T1*X;
Dv = T1*V;

for n=0:N-2
    Xh(n*2+1,1) = Dx(n+1);
    Xh(n*2+2,1) = Dv(n+1);
end

syms a b kp kd l;

Ub = a*Dx+b*Dv;

T2 = inv(triu(-1*ones(N-1)));

Un = inv(T2)*Ub;

for n = 1:N-1
    B(n*2,n:n+1) = [-1 1];
    
end

for n=0:N-1
    K(n+1,n*2+1:n*2+4) = [-kp -kd kp kd];
end

K  = K(:,3:end-2);
A = diag(mod((1:N*2-3),2)==1,1)*1;
sI = eye(length(A)).*l;

char = sI-A+B*K;



