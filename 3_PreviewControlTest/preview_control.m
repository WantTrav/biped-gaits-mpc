function [comx,comy] = preview_control()
Tr= 1;
dt = 0.02;
num = 11;
Tn = num*round(Tr/dt);

Stepx = 0.1;
Stepwidth = 0.1;

Sx = Stepx*ones(1,num);
Sx(1) = 0;
Sx(2) = 0;
Sy = Stepwidth*ones(1,num);
Sy(1) = 0;
Sy(2) = 0;
Sy(3) = Stepwidth/2;

Np = 100;
n=10;
[zmpx,zmpy] = ZmpPlan(Tr,dt,n,Sx,Sy);
[P,K,f,A,b] = RiccatiCal(dt,Np);

comx = zeros(3,Tn);
comy = zeros(3,Tn);
ux = zeros(Tn,1);
uy = zeros(Tn,1);
for k =1:Tn
    zmpx_ref = zmpx(k+1:k+Np);
    zmpy_ref = zmpy(k+1:k+Np);
    
    ux(k) = -K*comx(:,k)+ f*zmpx_ref;
    comx(:,k+1) = A*comx(:,k)+b*ux(k);
    
    uy(k) = -K*comy(:,k)+ f*zmpy_ref;
    comy(:,k+1) = A*comy(:,k)+b*uy(k);
    
end
close all;
figure(1)
plot(comx(1,:),comy(1,:));

hold on
plot(zmpx,zmpy);
legend('COM¹ì¼£','ZMP')

figure(2)
plot(comx(1,:));
hold on
plot(zmpx);
legend('COM-x','ZMP-x');
figure(3)
plot(comy(1,:));
hold on
plot(zmpy);
axis([0 600 -0.06 0.06])
legend('COM-y','ZMP-y');
end
function [P,K,f,A,b] = RiccatiCal(dt,N)
close all;
zc = 1.2;
g  = 9.8;
A = [1 dt dt^2/2; 0 1 dt; 0 0 1];
b = [dt^3/6 ; dt^2/2 ; dt];
c = [1 0 -zc/g];
Q =1000; R=1;
E =eye(3);
S =[0;0;0];
syms P

[P,L,G] = dare(A,b,c'*Q*c,1,S,E);
% [P] = solve(P == A'*P*A+c'*Q*c-A'*P*b*(R+b'*P*b)^(-1)*b'*P*A);

K =(R+b'*P*b)\(b'*P*A);

f = zeros(1,N);
for i = 1:1:N
    f(i) = (R+b'*P*b)^(-1)*b'*((A-b*K)')^(i-1)*c'*Q;
end

t =0.02:0.02:0.02*N;
% figure(1)
% plot(t,f);
end

function [zmpx,zmpy] = ZmpPlan(Tr,dt,num1,Sx,Sy)
% dt = 0.02;
% num=10;
n_t_r = round(Tr/dt);
num = num1+3;
Footx = zeros(1,num);
Footy = zeros(1,num);
Sx(num1+1) = 0;
Sx(num1+2) = 0;
Sx(num1+3) = 0;
Sy(num1+1) = Sy(num1)/2;
Sy(num1+2) = 0;
Sy(num1+3) = 0;
for i=2:num
    Footx(i)=Footx(i-1)+Sx(i);
    Footy(i)=Footy(i-1)+Sy(i)*(-1)^(i+1);
end

zmpx = zeros(num*n_t_r,1);
zmpy = zeros(num*n_t_r,1);

for j =1:num
    zmpx((j-1)*n_t_r+1:j*n_t_r) = Footx(j);
    zmpy((j-1)*n_t_r+1:j*n_t_r) = Footy(j);
end
figure(1)
subplot(2,1,1);
plot(zmpx);
subplot(2,1,2);
plot(zmpy);
figure(2);
plot(zmpx,zmpy);
end