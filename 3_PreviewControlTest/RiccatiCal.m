function [P,K,f,t] = RiccatiCal(dt)
close all;
zc = 0.32;
g  = 9.8;
A = [1 dt dt^2/2; 0 1 dt; 0 0 1];
b = [dt^3/6 ; dt^2/2 ; dt];
c = [1 0 -zc/g];
Q =1000;R=1;
E =eye(3);
S =[0;0;0];
syms P

[P,L,G] = dare(A,b,c'*Q*c,1,S,E);
% [P] = solve(P == A'*P*A+c'*Q*c-A'*P*b*(R+b'*P*b)^(-1)*b'*P*A);

K =inv(R+b'*P*b)*b'*P*A;
N=1000;
f = zeros(1,N);
for i = 1:1:N
    f(i) = (R+b'*P*b)^(-1)*b'*((A-b*K)')^(i-1)*c'*Q;
end

t =0.02:0.02:0.02*N;
figure(1)
plot(t,f);
end