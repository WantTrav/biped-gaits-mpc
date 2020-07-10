function [zmpx,zmpy] = ZmpPlan(Tr,Stepx,Stepwidth)
dt = 0.02;
num=10;
n_t_r = round(Tr/dt);
Sx = Stepx*ones(1,num);
Sx(1) = 0;
Sx(2) = 0;
Sy = Stepwidth*ones(1,num);
Sy(1) = 0;
Sy(2) = 0;
Sy(3) = Stepwidth/2;
Footx = zeros(1,num);
Footy = zeros(1,num);

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