function [yyy] = nao_20180118_energy_optimal_thress_mass_dot1

alpha =0.25; %%%%双足相占比；
belta =0.3;%%%双足相期间ZMP唯一量（对于点式，只是一个固定点，没有移动）
T_ref =2;   %%%%周期
stepx =60; %%%%%步长
yyy=PDC_energy_calculation(alpha,belta,T_ref,stepx);


end
function [yyy]=PDC_energy_calculation(alpha,belta, T_ref, stepx)
v_gama = 0;
palse = 1;  stepwidth=100;
M = 2.841+0.141+0.390+0.301+0.134+0.172+0.141+0.390+0.301+0.134+0.172;     %%%5.117kg
%单足相ZMP容许变动范围，
Px_max = 70; Px_min = 30; Py_max = 30; Py_min = 23; 
% Px_max = 20; Px_min = 20; Py_max = 20; Py_min = 20; 
%% 约束条件
%约束条件1：稳定性条件，首要条件
Px_max1 = 102; Px_min1 = 30; Py_max1 = 30; Py_min1 =30; 
%约束条件2:关节角度限制
angle_max = [0.379,0.484,2.11,0.92,0.397,0.79,0.484,2.11,0.92,0.768]; 
angle_min = [-0.79,-1.535,-0.1,-1.186,-0.768,-0.379,-1.535,-0.1,-1.186,-0.397];
%几何约束条件3：
V_angle_max = [4.2,6.5,6.5,6.5,4.2,4.2,6.5,6.5,6.5,4.2];
%约束条件4：
ficion_co = 0.32; %地面摩擦系数；                    

%%约束条件判断和循环求解
nx = 16;
[t_goal,yx1,y1,bodyp21,com_multi1,zmp1,totZMP1,Rp1,Lp1,Lpp_ground1,Rpp_ground1,vel_angle1,acc_angle1,J_R1,J_L1,Co_R1,Co_L1]=nao_20170918_azr(palse,nx,stepwidth,stepx,T_ref,alpha,Px_max,Px_min,Py_max,Py_min,belta);
y11111 = y1';
bodyp21111= bodyp21'; 
[axx, bxx]= size(bodyp21111);
[axx1, bxx1]= size(zmp1);
det_bxxx = bxx-bxx1;
zmp1111 = bodyp21111;
zmp1111(1:2,det_bxxx+1:end)=zmp1;
%%%%角度保存
save nao_20180118_energy_optimal_thress_mass_dot1_angle.txt y11111 -ascii;
save nao_20180118_energy_optimal_thress_mass_dot1_bodyp2.txt bodyp21111 -ascii;
save nao_20180118_energy_optimal_thress_mass_dot1_zmp.txt zmp1111 -ascii;


dt = t_goal(2)-t_goal(1);
nn_des = 2*round(T_ref/dt);
zmp_ref_x = zmp1(1,end-nn_des+1:end); zmp_ref_y = zmp1(2,end-nn_des+1:end);
zmp_real_x = totZMP1(1,end-nn_des+1:end); zmp_real_y = totZMP1(2,end-nn_des+1:end);

xfoot3= 10*stepx;
yfoot3= -stepwidth/2;
xfoot4= 11*stepx;
yfoot4= stepwidth/2;
Td = round(alpha*T_ref/dt);
zmp_det_x1=zeros(1,nn_des/2-Td);  zmp_det_x2=zeros(1,nn_des/2-Td);
zmp_det_y1=zeros(1,nn_des/2-Td);  zmp_det_y2=zeros(1,nn_des/2-Td);
for i=Td+1:nn_des/2                       %%左脚支撑
    zmp_det_x1(i)=zmp_real_x(i)-xfoot3;
    zmp_det_y1(i)=zmp_real_y(i)-yfoot3;   % 
end
for i=nn_des/2+1+Td:nn_des                   %%右脚支撑
    zmp_det_x2(i)=zmp_real_x(i)-xfoot4;
    zmp_det_y2(i)=zmp_real_y(i)-yfoot4;   % 
end
max_zmp_det_x1 = max(zmp_det_x1);  min_zmp_det_x1 = min(zmp_det_x1); max_zmp_det_x2 = max(zmp_det_x2); min_zmp_det_x2 = min(zmp_det_x2);
max_zmp_det_y1 = max(zmp_det_y1);  min_zmp_det_y1 = min(zmp_det_y1); max_zmp_det_y2 = max(zmp_det_y2); min_zmp_det_y2 = min(zmp_det_y2);

yxx = y1(end-nn_des+1:end,:);
y_max = max(yxx);
y_min = min(yxx);                
vangl = vel_angle1(:,end-nn_des+1:end)';
velo_max =max(abs(vangl));

det_angle = zeros(nn_des,10);
for i=1:10
    for j=1:nn_des-1
        det_angle(nn_des-i+1,i)=(y1(end-j+1,i)-y1(end-j,i))/dt;
    end
end
det_angle_max = max(det_angle);

com = bodyp21';
vel_com = Diff(t_goal',com);
acc_com = Diff(t_goal',vel_com)/1000; %%%质心加速度：3*b: m/s^2
F = zeros(3,nn_des);
G = [0;0;-9.8];
for i = 1:nn_des
    F(:,nn_des-i+1) = acc_com(:,end-i+1)-G;
end
%%%%注意，切向力为x,y向的矢量和，法向力为z方向（已经独立了）:
real_co =zeros(1,nn_des);
for i = 1:nn_des
    real_co(i) = sqrt(F(1,i)^2+F(2,i)^2)/F(3,i);   
end
real_co_max = max(real_co);

if (((max_zmp_det_x1-Px_max1)<=0) && ((min_zmp_det_x1+Px_min1)>=0)&&((max_zmp_det_y1-Py_max1)<=0)&&((min_zmp_det_y1+Py_min1)>=0)...
        &&((max_zmp_det_x2-Px_max1)<=0) && ((min_zmp_det_x2+Px_min1)>=0)&&((max_zmp_det_y2-Py_min1)<=0)&&((min_zmp_det_y2+Py_max1)>=0)) 
    if (max(y_max-angle_max)<=0) && (min(y_min-angle_min)>=0)
        if (max(velo_max-V_angle_max)<=0)
            if (real_co_max - ficion_co)<=0
                flag_energy=0;
            else
%                 fprintf('%s\n','不满足条件4');
                flag_energy=4;
                               
            end
        else
%             fprintf('%s\n','不满足条件3');
            flag_energy=3;
                      
        end
    else
%         fprintf('%s\n','不满足条件2');
        flag_energy=2;
                     
    end    
else
%     fprintf('%s\n','不满足条件1');
    flag_energy=1;                  
end

%^%%%%无论是否满足边界条件都进行计算，都进行计算
com = bodyp21';
[a,b]=size(zmp1);
[a1,b1]=size(bodyp21');
n_t = b1-b+1;

%%%% 角度、角速度与角加速度
yx=yx1(n_t:end,:);    y=y1(n_t:end,:);    vel_angle = vel_angle1(:,n_t:end);   acc_angle = acc_angle1(:,n_t:end);
%%%% 左右脚的全局坐标系位置:3*1*b
Lfoot = Lp1(:,:,n_t:end);   Rfoot = Rp1(:,:,n_t:end); 
Lpp_ground = Lpp_ground1(:,:,n_t:end);
Rpp_ground = Rpp_ground1(:,:,n_t:end);
%%%% 质心位置: 3*b
com = com(:,n_t:end);
%%%% 左右脚的雅可比矩阵： 6*5*b
J_R = J_R1(:,:,n_t:end)/1000;    J_L = J_L1(:,:,n_t:end)/1000;  %以m为单位

%%%%这里主要计算每个时刻左右脚的力分配，
t_goal = t_goal(n_t:end);
% %%%% 规划的质心轨迹
vel_com = Diff(t_goal',com);
acc_com = Diff(t_goal',vel_com)/1000; %%%质心加速度：3*b: m/s^2
%%%%%%这里求解质心移动轨迹:x相移动到移动到500mm为界限
picste=find(30>=com(1,:));
n_start = picste(end);
picene = find(530>com(1,:));
n_goal = picene(end);
%%%计算全局反力矢量：
F=zeros(3,n_goal);
G= [0;0;-9.8];
for i=n_start+1:n_goal
    F(:,i) = M*(acc_com(:,i)-G);
end
%%%%%左、右脚力和力矩分配
F_R = zeros(6,n_goal); F_L = zeros(6,n_goal);
for i=n_start+1:n_goal
    F_R(1:3,i) = Co_R1(:,:,i)*F(:,i);
    F_L(1:3,i) = Co_L1(:,:,i)*F(:,i);
end
%%%%%%%%这里左右脚所受地面反力矩
%%%%%计算全局地面反力矩
GRm = zeros(3,n_goal);
for i=n_start:n_goal
    GRm(:,i) = cross((Rfoot(:,:,i)-com(:,i))/1000,F_R(1:3,i))+cross((Lfoot(:,:,i)-com(:,i))/1000,F_L(1:3,i));
    F_R(4:6,i) = Co_R1(:,:,i)*GRm(:,i);
    F_L(4:6,i) = Co_L1(:,:,i)*GRm(:,i);
end

Torque = zeros(10,n_goal);
energy = 0;
dt = t_goal(2)-t_goal(1);
vv_joint = zeros(10,n_goal);
detz= zeros(3,1);
for i=2:n_goal
    det_Rp =  (Rfoot(:,:,i)-com(:,i))-(Rfoot(:,:,i-1)-com(:,i-1));
    vv_joint(1:5,i)= J_R(:,:,i)\[det_Rp;detz];
    det_Lp =  (Lfoot(:,:,i)-com(:,i))-(Lfoot(:,:,i-1)-com(:,i-1));
    vv_joint(6:10,i)= J_L(:,:,i)\[det_Lp;detz];                                        
end                                
for i=n_start+1:n_goal
    Torque(1:5,i) = J_R(:,:,i)'*F_R(:,i);
    Torque(6:10,i) = J_L(:,:,i)'*F_L(:,i);
    Power = abs(vv_joint(:,i)')* abs(Torque(:,i));  
    energy =energy+abs(Power)*dt; 
end
W = energy/(t_goal(n_goal)-t_goal(n_start));                                                

yyy=[belta;T_ref;stepx;real_co_max;v_gama;flag_energy;energy;W];

%%
function dy=Diff(x,y)
    if length(size(y))==3
        [a,c,b]=size(y);
        y=reshape(y,[a,b]);
        cs=csapi(x,y);
        pp1=fnder(cs);
        dy=fnval(pp1,x);
        dy=reshape(dy,[a,1,b]);
    else 
        cs=csapi(x,y);
        pp1=fnder(cs);
        dy=fnval(pp1,x);
    end
end
end

function [t_goal,yx,y,bodyp2,com,zmp,bodyZMP,Rp,Lp,Lpp_ground,Rpp_ground,vel_angle,acc_angle,JR,JL,Co_R,Co_L] = nao_20170918_azr(palse,num,stepwidth,stepx,T_ref,alpha,Px_max,Px_min,Py_max,Py_min,belta)
global Link N1;
global bodypppp;
%global comx;
%间隔取为两秒，if else使得低电平输入时，不用求解，从而提高效率
dt=0.02;
[t_goal,zmp,com1,Lp,LR,Lv,Lw,Rp,RR,Rv,Rw,Lpp_ground,Rpp_ground,Co_R,Co_L]=plan_zmp_com_AZR(num,stepx,stepwidth,dt,T_ref,alpha,Px_max,Px_min,Py_max,Py_min,belta);
% zmp = [px,py,pz];
bodyp1=com1;
T_ini = 6;
n_t_ini =round(T_ini/dt);
t=t_goal;
N1=length(t);
if palse==1
    Initial(N1);
    for pii=1:2 %这个for循环用来修正误差ZMP，循环次数为修正次数
        if isempty(bodypppp)
            bodyp2=bodyp1';
            bodypppp=reshape(bodyp1,[3,1,N1]);            
        end
        [bp,bR,bv,bw]=trans_bp(bodypppp,t,N1);%trans输出的bodyp依然是三维数组
        Link(1).p=bp;Link(1).R=bR;Link(1).v=bv;Link(1).w=bw;
        R_RefpR.v=Rv;R_RefpR.w=Rw;R_RefpR.p=Rp;R_RefpR.R=RR;
        L_RefpR.v=Lv;L_RefpR.w=Lw;L_RefpR.p=Lp;L_RefpR.R=LR;
        Forward_Kinematics(2,N1);
        [JR]=InverseKinematics(6,R_RefpR);
        [JL]=InverseKinematics(11,L_RefpR);
        com=calCom();%计算总体质心
        [totZMP,bodyZMP]=cal_linkZMP(com,R_RefpR,L_RefpR,t);
        
        %%%-------------反馈修正ZMP反馈量-------
        zmp = zmp(1:2,:);
        n_ini = round(T_ref/dt);
        [zmp_x,zpm_y]=size(bodyZMP);
        zmp_ref = zeros(zmp_x,zpm_y);
        zmp_ref(1,1:n_t_ini) = bodyZMP(1,1:n_t_ini); 
        zmp_ref(1,n_t_ini+1:end) = zmp(1,n_ini+1:end);
        zmp_ref(2,1:n_t_ini) = bodyZMP(2,1:n_t_ini); 
        zmp_ref(2,n_t_ini+1:end) = zmp(2,n_ini+1:end);
        [del_com,del_p]=cal_del_x(bodyZMP,t,zmp_ref);
        bodypppp=bodypppp+reshape(del_com,3,1,[]); 
    end 
    c=zeros(N1,10);
    for j=2:11
        c(:,j-1)=Link(j).q(:);%Nx10
    end
    y=c;
else 
    y=zeros(N1,10);
end
y=y(:);
clear global bodypppp;
y=reshape(y,[],10);

yx=[t',y];
y1= reshape(y,[10,length(t)]);   %%%10*n
y1x =y';
vel_angle = Diff(t',y1x);         %%%t=n*1,
acc_angle = Diff(t',vel_angle);
function dy=Diff(x,y)
    if length(size(y))==3
        [a,c,b]=size(y);
        y=reshape(y,[a,b]);
        cs=csapi(x,y);
        pp1=fnder(cs);
        dy=fnval(pp1,x);
        dy=reshape(dy,[a,1,b]);
    else 
        cs=csapi(x,y);
        pp1=fnder(cs);
        dy=fnval(pp1,x);
    end
end

% figure (8)
% hold on;
% plot(com1(1,:),com1(2,:),'--');
% plot(com(1,:),com(2,:),'r');
% legend('design_com','multilink_com');
% figure (9)
% hold on;
% plot(zmp_ref(2,:));
% plot(totZMP(2,:),'--');
% plot(bodyZMP(2,:),'r');
% legend('zmp_ref','totZMPy','bodyZMPy');
end
function Initial(tN)
global Link
Link(1).name='Body';
Link(1).brotherID=0;
Link(1).childID=2;
Link(1).motherID=0;
Name={'Body','Rleg_hip_r','Rleg_hip_p','Rleg_knee','Rleg_ank_p','Rleg_ank_r',...
      'Lleg_hip_r','Lleg_hip_p','Lleg_knee','Lleg_ank_p','Lleg_ank_r'};
BrotherID=[0 7 0 0 0 0 0 0 0 0 0];
ChildID=[2 3 4 5 6 0 8 9 10 11 0];
MotherID=[0 1 2 3 4 5 1 7 8 9 10];
a_Init={[1 0 0],[0 1 0],[0 1 0],[0 1 0],[1 0 0],...
   [1 0 0],[0 1 0],[0 1 0],[0 1 0],[1 0 0]};
b_Init={[0 -50 -85],[0 0 0],[0 0 -100],[0 0 -102.9],[0 0 0],...
    [0 50 -85],[0 0 0],[0 0 -100],[0 0 -102.9],[0 0 0]};    %mm
c_Init={[-1.44,-1.29,54.92],[-15.49,-0.29,-5.15],[1.38,-2.21,-53.73],...
        [4.53,-2.25,-49.36],[0.45,-0.29,6.85],[25.42,-3.3,-32.39],...
       [-15.49,0.29,-5.15],[1.38,2.21,-53.73],[4.53,2.25,-49.36],...
       [0.45,0.29,6.85],[25.42,3.3,-32.39]};  %mm
m_Init=[2.841,0.141,0.390,0.301,0.134,0.172,0.141,0.390,0.301,0.134,0.172];%kg
dq_Init=zeros(1,11);%link(1)的dq没有意义。

I1=[13953.66,6.19,-198.09;6.19,13318.88,-196.16;-198.09,-196.16,2682.42];%kg*mm^2
I2=[2.76,-0.02,-4.11;-0.02,9.83,0.00;-4.11,0.00,8.81];
I3=[1637.48,-0.92,85.88;-0.92,1592.21,-39.18;85.88,-39.18,303.98];
I4=[1182.83,-0.90,28.00,;-0.90,1128.28,-38.48;28.00,-38.48,191.45];
I5=[38.51,-0.06,3.87;-0.06,74.31,-0.00;3.87,-0.00,54.91];
I6=[269.30,5.88,139.13;5.88,643.47,-18.85;139.13,-18.85,525.03];
I7=[2.76,-0.02,-4.08;-0.02,9.83,-0.00;-4.08,-0.00,8.81];
I8=[1637.20,-0.92,85.31;-0.92,1591.07,-38.36;85.31,-38.36,303.74];
I9=[1182.08,0.63,36.50,;0.63,1128.65,39.50;36.50,39.50,193.22];
I10=[38.51,-0.03,3.86;-0.03,74.27,-0.02;3.86,-0.02,54.87];
I11=[269.44,-5.70,139.38;-5.70,644.43,18.74;139.38,18.74,525.76];
I_Init={I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11};
v_Init=zeros(11,1);
w_Init=zeros(11,1);
q=[0 0 0 0.2 0 0 0 0 0.2 0 0];%第一个是没有用的
%这里q和I都在第三维复制。
for j=1:1:11
    Link(j).name=Name(j);
    Link(j).brotherID=BrotherID(j);
    Link(j).childID=ChildID(j);
    Link(j).motherID=MotherID(j);
    Link(j).c=repmat(c_Init{1,j},[1,1,tN]);
    Link(j).m=m_Init(j);
    Link(j).dq=dq_Init(j);
    Link(j).I=repmat(I_Init{1,j},[1,1,tN]);%把矩阵放到cell里面，要注意cell的元素取用
    Link(j).v=v_Init(j);
    Link(j).w=w_Init(j);
    Link(j).q=repmat(q(j),[1,1],tN);
    if j>1
    Link(j).a=a_Init{1,j-1}';
    Link(j).b=b_Init{1,j-1}';
    end
end
end

%%%%生成Zmp，质心,足部位姿
function [txxx,zmp,com,Lpx,LRx,Lvx,Lwx,Rpx,RRx,Rvx,Rwx,Lpp_ground,Rpp_ground,Co_R,Co_L]=plan_zmp_com_AZR(num,steplength,stepwidth,dt,T_ref,alpha,Px_max,Px_min,Py_max,Py_min,belta)
T=T_ref; 
Tdd = alpha*T_ref;
n_td = round(Tdd/dt);
Tdou = n_td*dt;
%---静止状态到稳定状态的过渡
T_ini=6;
t0=dt:dt:T_ini;
COMy0=55;

[px,py,pz,comx1,comy1,Co_R,Co_L,n_jieduan]=nao_ZMP2CoM(num,dt,steplength,stepwidth,T_ref,alpha,Px_max,Px_min,Py_max,Py_min,belta);

tn=length(comx1);
N1=fix(T/dt);

Dest1=[dt  T_ini/4 T_ini/2 3*T_ini/4  T_ini-Tdou T_ini+dt T_ini+2*dt];       %初始阶段中只完成右脚向前跨一步，并且质心也移动到右脚上
x1=[0  -5 -10 -15  (comx1(N1+1)-15)/2  comx1(N1+1) comx1(N1+2)];
y1=[0 50 COMy0 COMy0   50 comy1(N1+1) comy1(N1+2)]; 
pp=pchip(Dest1,x1);
comx0=ppval(pp,t0);
pp2=pchip(Dest1,y1);
comy0=ppval(pp2,t0);
comx1=comx1(:,(N1+1):tn);
comy1=comy1(:,(N1+1):tn);
comx=[comx0,comx1];
comy=[comy0,comy1];

Destx=[0 0.01 0.02 0.03 T_ini-2-2*dt T_ini-2-dt T_ini-2 T_ini-1 T_ini];
z1=[333.09 333.0 332.95 332.85 310.31 310.1 310 310 310];
pp3=pchip(Destx,z1);
comz0=ppval(pp3,t0);
comz1=310*ones(1,tn-N1);
comz=[comz0,comz1];

%%%%%初始阶段利用CoM轨迹求解第1s的力分配系数。
for i=1:N1
    Co_L(:,:,i) = eye(3);
end

zmp=[px;py;pz];
com=[comx;comy;comz];
[Lpx,LRx,Lvx,Lwx,Rpx,RRx,Rvx,Rwx,Lpp_ground,Rpp_ground]=FootpR_polynomial_mod2(num,n_jieduan,dt,steplength,stepwidth,T_ref,alpha);


%%%%%%作图分析
tnnn =length(comx);
txxx=dt:dt:tnnn*dt;
Lpp=reshape(Lpx,3,tnnn);
Rpp=reshape(Rpx,3,tnnn);

% figure (6)
% hold on;
% plot(comx,comy,'--');
% plot(Lpp(1,:),Lpp(2,:),'--');
% plot(Rpp(1,:),Rpp(2,:));
% legend('com','lfoot','rfoot'); 
% 
% figure (7)
% subplot(2,1,1);
% hold on;
% plot(txxx,comx,'r');
% plot(txxx,Lpp(1,:),'--');
% plot(txxx,Rpp(1,:));
% legend('comx','rfootx','lfootx'); 
% subplot(2,1,2);
% hold on;
% plot(txxx,comy,'r');
% plot(txxx,Lpp(2,:),'--');
% plot(txxx,Rpp(2,:));
% legend('comy','rfooty','lfooty'); 

end
function [px,py,pz,comx,comy,Co_R,Co_L,n_jieduan]=nao_ZMP2CoM(n,dt,steplength,stepwidth,T_ref,alpha,Px_max,Px_min,Py_max,Py_min,belta)
global Ts;
global Td;
global T1;

%%
n_jieduan = 3;
Zc=310;   g=9800;   w=sqrt(g/Zc);
nt = n_jieduan*T_ref;
[AN]=round(nt/dt);
%周期时间
Tdd = alpha*T_ref;
n_td = round(Tdd/dt);
Tdd = n_td*dt;
Ts=(T_ref-Tdd)*ones(1,n);    Td=Tdd*ones(1,n);

%%%%%利用FSR位置容许裕度和最优速度
d2=Px_min*belta;   w1=Py_min*belta;  

T1 = Ts+Td;  
Kt=zeros(1,n);
for k=2:n
    Kt(k)=Kt(k-1)+round(T1(k-1)/dt);
end
N_sum = Kt(end)+round(T1(end)/dt);

%% 
tx_ini=Kt(2)+1;
[zmp,com] = nao_20171123_optimal_ZMP_COM(w,dt,T_ref,Tdd,n,n_jieduan,steplength,stepwidth,d2,w1);
px1 = zmp(1,:);
py1 = zmp(2,:);
pz = zeros(1,N_sum-AN);
comx1 = com(1,:);
comy1 = com(2,:);


%逆矩阵求解，相当于离线生成
t_des=dt:dt:(N_sum-AN)*dt;
px=px1; py=py1; comx= comx1; comy=comy1;

%%%%%力分配的矩阵分析：
Co_R =  zeros(3,3,N_sum-AN); Co_L =  zeros(3,3,N_sum-AN);
%%%%%注意，第一个周期也就是前1S并没有设计ZMP轨迹，所以无法计算力分配轨迹，这个要放到CoM外层计算
%%%%%注意，从1s中开始，ZMP开始于双足相，从左脚移动到右脚，然后开始右脚的单足支撑相。
for i=tx_ini:(N_sum-AN)
    b = find(i>Kt);
    i_zhouqi = b(end);    
    ki = Kt(i_zhouqi);   Kd = round(Td(i_zhouqi)/dt);   %Ks = round(Ts(i_zhouqi)/dt);    
    if   i<=(ki+Kd)                        %%%%双足相
        if rem(i_zhouqi,2)==0                     %右足支撑
            %%%%选择改进的方法
            if px(ki+Kd)-px(ki)==0                %双足相期间ZMP理想轨迹与y轴平行                
                Co_R(2,2,i) = abs((py(i)-py(ki+1))/(py(ki+Kd)-py(ki+1)));
                if Co_R(2,2,i) >1;
                    Co_R(2,2,i)=1;
                end                
                Co_R(1,1,i) = Co_R(2,2,i);                
                Co_R(3,3,i) = sqrt((Co_R(1,1,i)^2+Co_R(2,2,i)^2)/2);                  
            else
                %%%%%%结合参考文献，利用ZMP投影求解
                Co_R(1,1,i) = abs(((py(ki+Kd)-py(ki+1))*((py(i)-py(ki+1)))+(px(ki+Kd)-px(ki))*(px(i)-px(ki)))/((py(ki+Kd)-py(ki+1))^2+(px(ki+Kd)-px(ki))^2));
                if Co_R(1,1,i) >1;
                    Co_R(1,1,i)=1;
                end
                Co_R(2,2,i) = Co_R(1,1,i);
                Co_R(3,3,i) = sqrt((Co_R(1,1,i)^2+Co_R(2,2,i)^2)/2);                
            end               
            Co_L(1,1,i) = 1-Co_R(1,1,i);            
            Co_L(2,2,i) = 1-Co_R(2,2,i);
            Co_L(3,3,i) = 1-Co_R(3,3,i);
        else

            if px(ki+Kd)-px(ki)==0                %双足相期间ZMP理想轨迹与y轴平行                
                Co_L(2,2,i) = abs((py(i)-py(ki+1))/(py(ki+Kd)-py(ki+1)));
                if Co_L(2,2,i) >1;
                    Co_L(2,2,i)=1;
                end                
                Co_L(1,1,i) = Co_L(2,2,i);                
                Co_L(3,3,i) = sqrt((Co_L(1,1,i)^2+Co_L(2,2,i)^2)/2);                  
            else
                %%%%%%结合参考文献，利用ZMP投影求解
                Co_L(1,1,i) = abs(((py(ki+Kd)-py(ki+1))*((py(i)-py(ki+1)))+(px(ki+Kd)-px(ki))*(px(i)-px(ki)))/((py(ki+Kd)-py(ki+1))^2+(px(ki+Kd)-px(ki))^2));
                if Co_L(1,1,i) >1;
                    Co_L(1,1,i)=1;
                end
                Co_L(2,2,i) = Co_L(1,1,i);
                Co_L(3,3,i) = sqrt((Co_L(1,1,i)^2+Co_L(2,2,i)^2)/2);                
                     
            end              
            Co_R(1,1,i) = 1-Co_L(1,1,i);     
            Co_R(2,2,i) = 1-Co_L(2,2,i);             
            Co_R(3,3,i) = 1-Co_L(3,3,i);            
        end
    else                                   %%%%单足相
        if rem(i_zhouqi,2)==0                     %右足支撑
            Co_R(:,:,i) = eye(3);
        else
            Co_L(:,:,i) = eye(3);
        end
                
    end
end


figure (2)
subplot(2,1,1);
hold on;
plot(t_des,comx);
plot(t_des,px,'r.');
legend('comx','zmpx');
subplot(2,1,2);
hold on;
plot(t_des,comy);
plot(t_des,py,'r.');
legend('comy','zmpy');


end
function [zmp,com,comv,comacc] = nao_20171123_optimal_ZMP_COM(w,dt,T_ref,det_t,nxx,n_jieduan,steplength,stepwidth,a_n,b_n)
stepwidth = stepwidth/2;
nt = n_jieduan*T_ref;
n_t_n = round(T_ref/dt);
[AN]=round(nt/dt);
N_sum = nxx*n_t_n;

[zmp1,zmp_eq1,com1,comv1,comacc1,comaccx1,comxx1,comvx1,E1,Exx1,footx,Sx,comx,zmpx]=nao_20180111_zmp_comx_three_mass(dt,steplength,nxx,a_n,T_ref,det_t);
[zmp2,zmp_eq2,com2,comv2,comacc2,comaccx2,comxx2,comvx,E2,Exx2,footy,Sy,comy,zmpy]=nao_20180111_zmp_comy_three_mass(dt,stepwidth,nxx,b_n,T_ref,det_t);
zmp = [zmpx(1,1:N_sum-AN);zmpy(1,1:N_sum-AN)];
com = [comx(1,1:N_sum-AN);comy(1,1:N_sum-AN)];
%%%%%针对单足相中间时刻的质心不连续：（前后取四个点，进行拟合）：
n_det_dsp = round(det_t/dt);
n_det_ssp = n_t_n-n_det_dsp;
n_det_ssp_half = round(n_det_ssp/2);
for xx = 2:nxx-n_jieduan
    xxn = (xx-1)*n_t_n + n_det_dsp;
    if rem(n_det_ssp,2)==0    %%%%单足相为偶数
        x1 = com(1,xxn+n_det_ssp_half-2);   
        x2 = com(1,xxn+n_det_ssp_half+3);
        vx1 = (com(1,xxn+n_det_ssp_half-2) - com(1,xxn+n_det_ssp_half-3))/dt;   
        vx2 = (com(1,xxn+n_det_ssp_half+4) - com(1,xxn+n_det_ssp_half+3))/dt;
        xxx = [x1';x2';vx1';vx2'];
        A = [0,0,       0,    1;...
            (5*dt)^3,(5*dt)^2,(5*dt),1;...
            0,0,        1,    0;...
            3*(5*dt)^2,2*(5*dt), 1,    0];
        B = A\xxx;
        for j = 1:6
            AA = [((j-1)*dt)^3,((j-1)*dt)^2,((j-1)*dt),1];
            com(1,xxn+n_det_ssp_half-2+j-1) = (AA*B)';    
        end
        

        x1 = com(2,xxn+n_det_ssp_half-2);   
        x2 = com(2,xxn+n_det_ssp_half+3);
        vx1 = (com(2,xxn+n_det_ssp_half-2) - com(2,xxn+n_det_ssp_half-3))/dt;   
        vx2 = (com(2,xxn+n_det_ssp_half+4) - com(2,xxn+n_det_ssp_half+3))/dt;
        xxx = [x1';x2';vx1';vx2'];
        A = [0,0,       0,    1;...
            (5*dt)^3,(5*dt)^2,(5*dt),1;...
            0,0,        1,    0;...
            3*(5*dt)^2,2*(5*dt), 1,    0];
        B = A\xxx;
        for j = 1:6
            AA = [((j-1)*dt)^3,((j-1)*dt)^2,((j-1)*dt),1];
            com(2,xxn+n_det_ssp_half-2+j-1) = (AA*B)';    
        end        
               
    else                      %%%%单足相为奇数
        x1 = com(:,xxn+n_det_ssp_half-2);   
        x2 = com(:,xxn+n_det_ssp_half+2);
        vx1 = (com(:,xxn+n_det_ssp_half-2) - com(:,xxn+n_det_ssp_half-3))/dt;   
        vx2 = (com(:,xxn+n_det_ssp_half+3) - com(:,xxn+n_det_ssp_half+2))/dt;
        xxx = [x1';x2';vx1';vx2'];
        A = [0,0,       0,    1;...
            (4*dt)^3,(4*dt)^2,(4*dt),1;...
            0, 0,        1,    0;...
            3*(4*dt)^2,2*(4*dt), 1,    0];
        B = A\xxx;
        for j = 1:5
            AA = [((j-1)*dt)^3,((j-1)*dt)^2,((j-1)*dt),1];
            com(:,xxn+n_det_ssp_half-2+j-1) = (AA*B)';    
        end  
        
    end

end
%%%zmp
for xx = 2:nxx-n_jieduan
    xxn = (xx-1)*n_t_n + n_det_dsp;
    if rem(n_det_ssp,2)==0    %%%%单足相为偶数
        x1 = zmp(1,xxn+n_det_ssp_half-2);   
        x2 = zmp(1,xxn+n_det_ssp_half+3);
        vx1 = (zmp(1,xxn+n_det_ssp_half-2) - zmp(1,xxn+n_det_ssp_half-3))/dt;   
        vx2 = (zmp(1,xxn+n_det_ssp_half+4) - zmp(1,xxn+n_det_ssp_half+3))/dt;
        xxx = [x1';x2';vx1';vx2'];
        A = [0,0,       0,    1;...
            (5*dt)^3,(5*dt)^2,(5*dt),1;...
            0,0,        1,    0;...
            3*(5*dt)^2,2*(5*dt), 1,    0];
        B = A\xxx;
        for j = 1:6
            AA = [((j-1)*dt)^3,((j-1)*dt)^2,((j-1)*dt),1];
            zmp(1,xxn+n_det_ssp_half-2+j-1) = (AA*B)';    
        end
        

        x1 = zmp(2,xxn+n_det_ssp_half-2);   
        x2 = zmp(2,xxn+n_det_ssp_half+3);
        vx1 = (zmp(2,xxn+n_det_ssp_half-2) - zmp(2,xxn+n_det_ssp_half-3))/dt;   
        vx2 = (zmp(2,xxn+n_det_ssp_half+4) - zmp(2,xxn+n_det_ssp_half+3))/dt;
        xxx = [x1';x2';vx1';vx2'];
        A = [0,0,       0,    1;...
            (5*dt)^3,(5*dt)^2,(5*dt),1;...
            0,0,        1,    0;...
            3*(5*dt)^2,2*(5*dt), 1,    0];
        B = A\xxx;
        for j = 1:6
            AA = [((j-1)*dt)^3,((j-1)*dt)^2,((j-1)*dt),1];
            zmp(2,xxn+n_det_ssp_half-2+j-1) = (AA*B)';    
        end        
               
    else                      %%%%单足相为奇数
        x1 = zmp(:,xxn+n_det_ssp_half-2);   
        x2 = zmp(:,xxn+n_det_ssp_half+2);
        vx1 = (zmp(:,xxn+n_det_ssp_half-2) - zmp(:,xxn+n_det_ssp_half-3))/dt;   
        vx2 = (zmp(:,xxn+n_det_ssp_half+3) - zmp(:,xxn+n_det_ssp_half+2))/dt;
        xxx = [x1';x2';vx1';vx2'];
        A = [0,0,       0,    1;...
            (4*dt)^3,(4*dt)^2,(4*dt),1;...
            0, 0,        1,    0;...
            3*(4*dt)^2,2*(4*dt), 1,    0];
        B = A\xxx;
        for j = 1:5
            AA = [((j-1)*dt)^3,((j-1)*dt)^2,((j-1)*dt),1];
            zmp(:,xxn+n_det_ssp_half-2+j-1) = (AA*B)';    
        end  
        
    end

end
comv = [comv1(4,:);comv2(4,:)];
comacc = [comacc1(4,:);comacc2(4,:)];

end
function [zmp,zmp_eq,com,comv,comacc,comaccx,comxx,comvx,E,Exx,footx,Sx,comx,zmpx]=nao_20180111_zmp_comx_three_mass(dt,steplength,nxx,a_n,T_ref,det_t)
zmp = zeros(4,10);zmp_eq = zeros(4,10);com = zeros(4,10);comv = zeros(4,10);comacc = zeros(4,10);
comaccx = zeros(4,1);comxx = zeros(4,1);comvx = zeros(4,1);
E = zeros(4,1);Exx = zeros(4,1);

%%%%摆动足和支撑足部运动轨迹；
r_a = zeros(4,1);
l_a = zeros(4,1);
r_v = zeros(4,1);
l_v = zeros(4,1);
r_acc = zeros(4,1);
l_acc = zeros(4,1);
Er=zeros(4,1);
El=zeros(4,1);

%%%%左、右腿部摆动腿足部轨迹；
R_foot = zeros(4,10);
L_foot = zeros(4,10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%三质点问题求解初始条件：
g = 9800; Zc = 310; 
MB =2.841; ML=1.138; MM=MB+2*ML;
fb = 45.2;
a = [0 0 0;0 0 0;0 0 -100;0 0 -202.9;0 0 -202.9];
b = [-15.49,-0.29,-5.15;1.38,-2.21,-53.73;4.53,-2.25,-49.36;4.53,-2.25,-49.36;25.42,-3.3,-32.39];
c = [0.141,0.390,0.301,0.134,0.172];
ss = (c(1)*(a(1,:)+b(1,:))+c(2)*(a(2,:)+b(2,:))+c(3)*(a(3,:)+b(3,:))+c(4)*(a(4,:)+b(4,:))+c(5)*(a(5,:)+b(5,:)))/(ML);
Cxb = 0; Cyb = b(1,2); Czb = Zc;


%%周期
t_ref = zeros(nxx,2);  t_ref(:,1) = T_ref/2*ones(nxx,1); t_ref(:,2) = T_ref*ones(nxx,1);  
   
%%%步长与步宽
% steplength=steplength-2;
Sx = steplength*ones(1,nxx);
Sy = 100*zeros(1,nxx);

%%%%%足部支撑点交替改变
footx = zeros(1,nxx);
footy = zeros(1,nxx);
for j=1:nxx
    if j == 1
        footx(1) = 0;
    else
        footx(j) = footx(j-1) + Sx(j-1);
    end
end

%%%%等效ZMP位移
ax = zeros(4,nxx);ay=zeros(4,nxx);
cx=zeros(2,nxx); cy=zeros(2,nxx);

%%%%%参考ZMP前半单足相终止点和后半双足相终止点
xfoot =zeros(1,2);
for j = 1:nxx
    %%%第一步是右脚开始支撑（右脚支撑前半单足相+过渡相+左脚支撑后半单足相）：左右支撑出现对称性，所需没有必要区分偶数还是奇数
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================无双足相：t_i ==t_f===========
    % =============================================
    t_i = t_ref(j,1); t_f = t_ref(j,1);t_fn = t_ref(j,2); 
    n_t_i = round(t_i/dt); n_t_f = round(t_f/dt);

    n_t_n = round(t_fn/dt);
    if j == 1
        k_n_t_n = 0;
    else
        k_n_t_n = k_n_t_n + n_t_n;
    end
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 0;
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sx(j)-a_n;    

    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = 0;     x_if = a_n;   x_velo_i = 0;                              %相对于当前坐标footx(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：2次多项式： 
    xfoot_ii = 0;     xfoot_if = Sx(j);   xfoot_velo_ii = 0;               %相对于当前坐标footx(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1); foota_i3 = 0;
    
    %%%%%%等效ZMP轨迹公式求解=============
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else 
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end
    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=-foota_i3*ML/(2*MM)*(MM/(MB+ML));   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-(foota_i1-Czs/g*6*foota_i3)*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ex)*(MM/(MB+ML));
    B = [exp(-w*0),exp(w*0);exp(-w*t_i),exp(w*t_i)];
    COMx_is = 10;COMx_es = Sx(j)/2+10;
    comx_condi1=COMx_is-(ax(1,j)*(0)^3+ax(2,j)*(0)^2+ax(3,j)*(0)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*0+2*ax(2,j)));
    comx_condi2=COMx_es-(ax(1,j)*(t_i)^3+ax(2,j)*(t_i)^2+ax(3,j)*(t_i)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*t_i+2*ax(2,j)));
    cx(:,j) = B\[comx_condi1;comx_condi2];  %cy(:,i) = B\comy_condi;    
    for i=1:n_t_i
        zmp(1,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footx(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(1,k_n_t_n+i) = footx(j);
            R_foot(1,k_n_t_n+i) = foota_i3*(dt*i)^3+foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);        
        else
            R_foot(1,k_n_t_n+i) = footx(j);
            L_foot(1,k_n_t_n+i) = foota_i3*(dt*i)^3+foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);            
        end 
        tj = dt*i;
        com(1,k_n_t_n+i)=cx(1,j)*exp(-w*tj)+cx(2,j)*exp(w*tj)+ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*tj+2*ax(2,j))+footx(j);            
        comv(1,k_n_t_n+i) = -w*cx(1,j)*exp(-w*tj)+w*cx(2,j)*exp(w*tj)+3*ax(1,j)*(tj)^2+2*ax(2,j)*(tj)^1+ax(3,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j));
        comacc(1,k_n_t_n+i) = w^2*cx(1,j)*exp(-w*tj)+w^2*cx(2,j)*exp(w*tj)+6*ax(1,j)*(tj)^1+2*ax(2,j);
        zmp_eq(1,k_n_t_n+i) = ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+footx(j);
        comaccx(1) = comaccx(1)+ comacc(1,k_n_t_n+i)^2*dt;
        comxx(1) = comxx(1)+ com(1,k_n_t_n+i)^2;
        comvx(1) = comvx(1)+ comv(1,k_n_t_n+i)^2;
        E(1) = E(1) + abs(comv(1,k_n_t_n+i)* comacc(1,k_n_t_n+i)*dt);            

        r_a(1,k_n_t_n+i) = (R_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        l_a(1,k_n_t_n+i) = (L_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 3*foota_i3*(dt*i)^2+2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            r_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 6*foota_i3*(dt*i)^1+2*foota_i2)/2;
            l_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;  
        else     
            l_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 3*foota_i3*(dt*i)^2+2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            l_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 6*foota_i3*(dt*i)^1+2*foota_i2)/2;
            r_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;             
        end
        Er(1)=Er(1)+abs(r_v(1,k_n_t_n+i)* r_acc(1,k_n_t_n+i)*dt);
        El(1)=El(1)+abs(l_v(1,k_n_t_n+i)* l_acc(1,k_n_t_n+i)*dt);   
        
    end
    
    %%%%==========================================质心轨迹后半段
    % =============================================
    x_fi = -a_n;  x_ff = 0;     x_velo_f = 0;                              %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi =-Sx(j);xfoot_ff =0;       xfoot_velo_f = 0;                 %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);       

    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数            
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ex)*(MM/(MB+ML))+Sx(j);  %%%加上Sx(j)，转化为相对 坐标式footx(j) 

    B = [exp(-w*t_f),exp(w*t_f);exp(-w*t_fn),exp(w*t_fn)];
    COMx_is = Sx(j)/2+10;COMx_es = Sx(j)+10;
    comx_condi1=COMx_is-(ax(1,j)*(t_f)^3+ax(2,j)*(t_f)^2+ax(3,j)*(t_f)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*t_f+2*ax(2,j)));
    comx_condi2=COMx_es-(ax(1,j)*(t_fn)^3+ax(2,j)*(t_fn)^2+ax(3,j)*(t_fn)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*t_fn+2*ax(2,j)));
    cx(:,j) = B\[comx_condi1;comx_condi2];  %cy(:,i) = B\comy_condi;      

    for i=n_t_f+1:n_t_n      
        zmp(1,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footx(j)+Sx(j); 
        if rem(j,2)==0 %右脚支撑
            L_foot(1,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            R_foot(1,k_n_t_n+i) = footx(j)+Sx(j);          
        else
            R_foot(1,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            L_foot(1,k_n_t_n+i) = footx(j)+Sx(j);   
        end        
        tj = dt*i;
        com(1,k_n_t_n+i)=cx(1,j)*exp(-w*tj)+cx(2,j)*exp(w*tj)+ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*tj+2*ax(2,j))+footx(j);            
        comv(1,k_n_t_n+i) = -w*cx(1,j)*exp(-w*tj)+w*cx(2,j)*exp(w*tj)+3*ax(1,j)*(tj)^2+2*ax(2,j)*(tj)^1+ax(3,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j));
        comacc(1,k_n_t_n+i) = w^2*cx(1,j)*exp(-w*tj)+w^2*cx(2,j)*exp(w*tj)+6*ax(1,j)*(tj)^1+2*ax(2,j);
        zmp_eq(1,k_n_t_n+i) = ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+footx(j);
        comaccx(1) = comaccx(1)+ comacc(1,k_n_t_n+i)^2*dt;
        comxx(1) = comxx(1)+ com(1,k_n_t_n+i)^2;
        comvx(1) = comvx(1)+ comv(1,k_n_t_n+i)^2;
        E(1) = E(1) + abs(comv(1,k_n_t_n+i)* comacc(1,k_n_t_n+i)*dt); 
        
        r_a(1,k_n_t_n+i) = (R_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        l_a(1,k_n_t_n+i) = (L_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            r_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;  
        else     
            l_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            l_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;             
        end
        Er(1)=Er(1)+abs(r_v(1,k_n_t_n+i)* r_acc(1,k_n_t_n+i)*dt);
        El(1)=El(1)+abs(l_v(1,k_n_t_n+i)* l_acc(1,k_n_t_n+i)*dt);         
        
    end
    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================有双足相,点式参考轨迹===========
    % =============================================
    %%%从中间到两边拓展:det_t 为过渡相（也就是双足相时间）
    det_t_n = round(det_t/dt);
    t_i = t_ref(j,1)-round(det_t_n/2)*dt; 
    n_t_i = round(t_i/dt);
    t_i = n_t_i*dt;
    t_f = t_i+det_t_n*dt;
    n_t_f = round(t_f/dt);    
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 6.5;
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sx(j)-a_n;    


    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = 0;     x_if = a_n;   x_velo_i = 0;                              %相对于当前坐标footx(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：（利用对称性），修正摆动腿足部轨迹为二次多项式
    xfoot_ii = 0;     xfoot_if = Sx(j);   xfoot_velo_ii = 0;               %相对于当前坐标footx(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1);
    
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end


    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-foota_i1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ex)*(MM/(MB+ML));

    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_i00 = ax(4,j); a_i10 = ax(3,j);  a_i20 = ax(2,j); 
    for i=1:n_t_i
        zmp(2,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footx(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(2,k_n_t_n+i) = footx(j);
            R_foot(2,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);               
        else
            R_foot(2,k_n_t_n+i) = footx(j);
            L_foot(2,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);              
        end        
    end
    
    %%%%========================================== 质心轨迹后半段
    % =============================================
    x_fi = -a_n;  x_ff = 0;     x_velo_f = 0;                              %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi =-Sx(j);xfoot_ff =0;       xfoot_velo_f = 0;                  %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);              
    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数  
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ex)*(MM/(MB+ML))+Sx(j);  %%%加上Sx(j)，转化为相对 坐标式footx(j) 
    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_f00 = ax(4,j); a_f10 = ax(3,j);  a_f20 = ax(2,j); 
    for i=n_t_f:n_t_n    
        zmp(2,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footx(j)+Sx(j); 
        if rem(j,2)==0 %右脚支撑
            L_foot(2,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            R_foot(2,k_n_t_n+i) = footx(j)+Sx(j);           
        else
            R_foot(2,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            L_foot(2,k_n_t_n+i) = footx(j)+Sx(j);   
        end               
    end     

    %%%% zmp轨迹过渡中间相
    for i=n_t_i+1:n_t_f-1
        %%% 3次多项式求解参考ZMP和质心轨迹。
        if rem(j,2)==0 %右脚支撑
            L_foot(2,k_n_t_n+i) = footx(j);
            R_foot(2,k_n_t_n+i) = footx(j)+Sx(j);           
        else
            R_foot(2,k_n_t_n+i) = footx(j);
            L_foot(2,k_n_t_n+i) = footx(j)+Sx(j);
        end

        %%%%ZMP轨迹初始和末端时刻位置和速度
        zmp_condi = [a_n+footx(j);-a_n+footx(j)+Sx(j);0;0];
        ttt_f = [(t_i)^3,(t_i)^2,(t_i)^1,1;(t_f)^3,(t_f)^2,t_f,1;3*(t_i)^2,2*(t_i),1,0;3*(t_f)^2,2*(t_f),1,0];    
        a_zmp = ttt_f\zmp_condi;
        tj = dt*i;
        zmp(2,k_n_t_n+i) = a_zmp(4) + a_zmp(3)*(tj) + a_zmp(2)*(tj)^2 + a_zmp(1)*(tj)^3;           
    end       
    
    %%%%%%%%%%======================================全局质心加速度最小=======================求解    
    a_i0 = a_i00; a_i1 = a_i10;  a_i2 = a_i20; 
    a_f0 = a_f00; a_f1 = a_f10;  a_f2 = a_f20; 
    xx_u_t0 = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_s_t0 = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_u_tn = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_tn = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);   
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    Ah = xx_u_t0 - xx_s_t0;
    Bh = xx_s_tn - xx_u_tn;
    % h_add1 = -(w^2)/4*[Ah*exp(-w*t_i);Bh*exp(w*t_f-w*t_fn)];
    h_add1 = 0;
    h_add2 = -w*[a_i2*(1-exp(-w*t_i));a_f2*(1-exp(w*t_f-w*t_fn))];
    h_add = h_add1 + h_add2 ;
    %%%%%%%%矩阵求解
    % Wpre1 = w^2/4*exp(-2*w*t_i);   
    % Wpost1 = w^2/4*exp(2*w*t_f-2*w*t_fn); 
    Wpre1 = 0;   
    Wpost1 = 0; 
    Wpre2 = w^3/8*(1-exp(-2*w*t_i));   
    Wpost2 = w^3/8*(1-exp(2*w*t_f-2*w*t_fn));   
    Wpre = Wpre1 + Wpre2 ;
    Wpost = Wpost1 + Wpost2;


    W = [Wpre,0;0,Wpost];
    G = [det_t^3/3,det_t^2/2;det_t^2/2,det_t]; 
    G_inv = inv(G);
    phi_u = [1/2;w/2];  phi_s = [1/2;-w/2];
    M = [1,det_t;0,1];
    H11 = -M*phi_s;  H1 = [H11,phi_u];
    H21 = -M*phi_u;  H2 = [H21,phi_s];

    %%%计算已知分量：
    %%%xx_s_t = a0 + (2*a2 - w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%xx_u_t = a0 + (2*a2 + w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%其中：xx_c_t = a0 + a1*t + a2*t^2 + (2*a2)/w^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    F = [xx_s_ti;xx_u_tf];
    A = W + H2'*G_inv*H2;
    h = (W - H2'*G_inv*H1)*F + h_add;
    Phi = A\h;

    x_u_ti = Phi(1);
    x_s_tf = Phi(2);


    %%%% 质心轨迹前半段；
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);    
    for i=1:n_t_i        
        com(2,k_n_t_n+i) = exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 + (2*a_i2)/w^2);
        comv(2,k_n_t_n+i) = w*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i1 + 2*a_i2*(dt*i));
        comacc(2,k_n_t_n+i) = w^2*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + 2*a_i2;
        zmp_eq(2,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2;
        comaccx(2) = comaccx(2)+ comacc(2,k_n_t_n+i)^2*dt;
        comxx(2) = comxx(2)+ com(2,k_n_t_n+i)^2;
        comvx(2) = comvx(2)+ comv(2,k_n_t_n+i)^2;
        E(2) = E(2) + abs(comv(2,k_n_t_n+i)* comacc(2,k_n_t_n+i)*dt);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i);
        
        r_a(2,k_n_t_n+i) = (R_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(2,k_n_t_n+i) = (L_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footx(j))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            r_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) +2*foota_i2)/2;
            l_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;  
        else     
            l_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            l_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) +2*foota_i2)/2;
            r_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;             
        end
        Er(2)=Er(2)+abs(r_v(2,k_n_t_n+i)* r_acc(2,k_n_t_n+i)*dt);
        El(2)=El(2)+abs(l_v(2,k_n_t_n+i)* l_acc(2,k_n_t_n+i)*dt);          
    end

    %%%% 质心轨迹后半段
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    for i=n_t_f+1:n_t_n
        com(2,k_n_t_n+i) = exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2+ (a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 + (2*a_f2)/w^2);
        comv(2,k_n_t_n+i) = -w*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + (a_f1 + 2*a_f2*(dt*i));
        comacc(2,k_n_t_n+i) = w^2*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + 2*a_f2;
        zmp_eq(2,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2;
        comaccx(2) = comaccx(2)+ comacc(2,k_n_t_n+i)^2*dt;
        comxx(2) = comxx(2)+ com(2,k_n_t_n+i)^2; 
        comvx(2) = comvx(2)+ comv(2,k_n_t_n+i)^2;
        E(2) = E(2) + abs(comv(2,k_n_t_n+i)* comacc(2,k_n_t_n+i)*dt);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i);
        
        r_a(2,k_n_t_n+i) = (R_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(2,k_n_t_n+i) = (L_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footx(j))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            r_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;  
        else     
            l_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            l_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;             
        end
        Er(2)=Er(2)+abs(r_v(2,k_n_t_n+i)* r_acc(2,k_n_t_n+i)*dt);
        El(2)=El(2)+abs(l_v(2,k_n_t_n+i)* l_acc(2,k_n_t_n+i)*dt);          
    end

    %%%%求解中间阶段的质心最优加速度和对应的质心轨迹和ZMP轨迹
    d_ti_tf = H2*Phi+H1*F;
    B = [0;1];
    for i=n_t_i+1:n_t_f  
        E_A_ti = [1,0;t_f-i*dt,1];
        comacc(2,k_n_t_n+i) = B'* E_A_ti * G_inv * d_ti_tf;
        %%%这里使用v_t1 = v_t2+(t1-t2)*a_t2;
        comv(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)*dt + comv(2,k_n_t_n+i-1);
        com(2,k_n_t_n+i) = comv(2,k_n_t_n+i)*dt + com(2,k_n_t_n+i-1);
        zmp_eq(2,k_n_t_n+i) = com(2,k_n_t_n+i) - comacc(2,k_n_t_n+i)/(w^2);
        comaccx(2) = comaccx(2)+ comacc(2,k_n_t_n+i)^2*dt;
        comxx(2) = comxx(2)+ com(2,k_n_t_n+i)^2;
        comvx(2) = comvx(2)+ comv(2,k_n_t_n+i)^2;
        E(2) = E(2) + abs(comv(2,k_n_t_n+i)* comacc(2,k_n_t_n+i)*dt);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i);
        
        r_a(2,k_n_t_n+i) = (R_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(2,k_n_t_n+i) = (L_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footx(j))/2+Cxs;
        r_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
        l_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
        r_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;
        l_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;  
        Er(2)=Er(2)+abs(r_v(2,k_n_t_n+i)* r_acc(2,k_n_t_n+i)*dt);
        El(2)=El(2)+abs(l_v(2,k_n_t_n+i)* l_acc(2,k_n_t_n+i)*dt);         
    end
    for i=1:n_t_n
        zmp_eq(2,k_n_t_n+i) = zmp_eq(2,k_n_t_n+i) +footx(j);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i) + footx(j);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================有双足相,直线式参考轨迹===========
    % =============================================
    %%%从中间到两边拓展:det_t 为过渡相（也就是双足相时间）
    det_t_n = round(det_t/dt);
    t_i = t_ref(j,1)-round(det_t_n/2)*dt; 
    n_t_i = round(t_i/dt);
    t_i = n_t_i*dt;
    t_f = t_i+det_t_n*dt;
    n_t_f = round(t_f/dt);    
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 10;
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sx(j)-a_n;    


    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = 0;     x_if = a_n;   x_velo_i = (x_if-x_ii)/t_i;                %相对于当前坐标footx(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：（利用对称性），修正摆动腿足部轨迹为二次多项式
    xfoot_ii = 0;     xfoot_if = Sx(j);   xfoot_velo_ii = 0;               %相对于当前坐标footx(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1);
    
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end


    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-foota_i1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ex)*(MM/(MB+ML));

    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_i00 = ax(4,j); a_i10 = ax(3,j);  a_i20 = ax(2,j); 
    for i=1:n_t_i
        zmp(3,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footx(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(3,k_n_t_n+i) = footx(j);
            R_foot(3,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);               
        else
            R_foot(3,k_n_t_n+i) = footx(j);
            L_foot(3,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);              
        end        
    end
    
    %%%%========================================== 质心轨迹后半段
    % =============================================
    x_fi = -a_n;  x_ff = 0;     x_velo_f = (x_ff-x_fi)/(t_fn-t_f);         %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi =-Sx(j);xfoot_ff =0;       xfoot_velo_f = 0;                  %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);              
    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数  
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ex)*(MM/(MB+ML))+Sx(j);  %%%加上Sx(j)，转化为相对 坐标式footx(j) 
    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_f00 = ax(4,j); a_f10 = ax(3,j);  a_f20 = ax(2,j); 
    for i=n_t_f:n_t_n    
        zmp(3,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footx(j)+Sx(j); 
        if rem(j,2)==0 %右脚支撑
            L_foot(3,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            R_foot(3,k_n_t_n+i) = footx(j)+Sx(j);           
        else
            R_foot(3,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            L_foot(3,k_n_t_n+i) = footx(j)+Sx(j);   
        end               
    end     

    %%%% zmp轨迹过渡中间相
    for i=n_t_i+1:n_t_f-1
        %%% 3次多项式求解参考ZMP和质心轨迹。
        if rem(j,2)==0 %右脚支撑
            L_foot(3,k_n_t_n+i) = footx(j);
            R_foot(3,k_n_t_n+i) = footx(j)+Sx(j);           
        else
            R_foot(3,k_n_t_n+i) = footx(j);
            L_foot(3,k_n_t_n+i) = footx(j)+Sx(j);
        end

        %%%%ZMP轨迹初始和末端时刻位置和速度
        zmp_condi = [a_n+footx(j);-a_n+footx(j)+Sx(j);0;0];
        ttt_f = [(t_i)^3,(t_i)^2,(t_i)^1,1;(t_f)^3,(t_f)^2,t_f,1;3*(t_i)^2,2*(t_i),1,0;3*(t_f)^2,2*(t_f),1,0];    
        a_zmp = ttt_f\zmp_condi;
        tj = dt*i;
        zmp(3,k_n_t_n+i) = a_zmp(4) + a_zmp(3)*(tj) + a_zmp(2)*(tj)^2 + a_zmp(1)*(tj)^3;           
    end       
    
    %%%%%%%%%%======================================全局质心加速度最小=======================求解    
    a_i0 = a_i00; a_i1 = a_i10;  a_i2 = a_i20; 
    a_f0 = a_f00; a_f1 = a_f10;  a_f2 = a_f20; 
    xx_u_t0 = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_s_t0 = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_u_tn = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_tn = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);   
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    Ah = xx_u_t0 - xx_s_t0;
    Bh = xx_s_tn - xx_u_tn;
    % h_add1 = -(w^2)/4*[Ah*exp(-w*t_i);Bh*exp(w*t_f-w*t_fn)];
    h_add1 = 0;
    h_add2 = -w*[a_i2*(1-exp(-w*t_i));a_f2*(1-exp(w*t_f-w*t_fn))];
    h_add = h_add1 + h_add2 ;
    %%%%%%%%矩阵求解
    % Wpre1 = w^2/4*exp(-2*w*t_i);   
    % Wpost1 = w^2/4*exp(2*w*t_f-2*w*t_fn); 
    Wpre1 = 0;   
    Wpost1 = 0; 
    Wpre2 = w^3/8*(1-exp(-2*w*t_i));   
    Wpost2 = w^3/8*(1-exp(2*w*t_f-2*w*t_fn));   
    Wpre = Wpre1 + Wpre2 ;
    Wpost = Wpost1 + Wpost2;


    W = [Wpre,0;0,Wpost];
    G = [det_t^3/3,det_t^2/2;det_t^2/2,det_t]; 
    G_inv = inv(G);
    phi_u = [1/2;w/2];  phi_s = [1/2;-w/2];
    M = [1,det_t;0,1];
    H11 = -M*phi_s;  H1 = [H11,phi_u];
    H21 = -M*phi_u;  H2 = [H21,phi_s];

    %%%计算已知分量：
    %%%xx_s_t = a0 + (2*a2 - w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%xx_u_t = a0 + (2*a2 + w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%其中：xx_c_t = a0 + a1*t + a2*t^2 + (2*a2)/w^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    F = [xx_s_ti;xx_u_tf];
    A = W + H2'*G_inv*H2;
    h = (W - H2'*G_inv*H1)*F + h_add;
    Phi = A\h;

    x_u_ti = Phi(1);
    x_s_tf = Phi(2);


    %%%% 质心轨迹前半段；
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);    
    for i=1:n_t_i        
        com(3,k_n_t_n+i) = exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 + (2*a_i2)/w^2);
        comv(3,k_n_t_n+i) = w*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i1 + 2*a_i2*(dt*i));
        comacc(3,k_n_t_n+i) = w^2*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + 2*a_i2;
        zmp_eq(3,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2;
        comaccx(3) = comaccx(3)+ comacc(3,k_n_t_n+i)^2*dt;
        comxx(3) = comxx(3)+ com(3,k_n_t_n+i)^2;
        comvx(3) = comvx(3)+ comv(3,k_n_t_n+i)^2;
        E(3) = E(3) + abs(comv(3,k_n_t_n+i)* comacc(3,k_n_t_n+i)*dt);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i);
        
        r_a(3,k_n_t_n+i) = (R_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(3,k_n_t_n+i) = (L_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footx(j))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            r_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) +2*foota_i2)/2;
            l_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;  
        else     
            l_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            l_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) +2*foota_i2)/2;
            r_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;             
        end
        Er(3)=Er(3)+abs(r_v(3,k_n_t_n+i)* r_acc(3,k_n_t_n+i)*dt);
        El(3)=El(3)+abs(l_v(3,k_n_t_n+i)* l_acc(3,k_n_t_n+i)*dt);          
    end

    %%%% 质心轨迹后半段
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    for i=n_t_f+1:n_t_n
        com(3,k_n_t_n+i) = exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2+ (a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 + (2*a_f2)/w^2);
        comv(3,k_n_t_n+i) = -w*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + (a_f1 + 2*a_f2*(dt*i));
        comacc(3,k_n_t_n+i) = w^2*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + 2*a_f2;
        zmp_eq(3,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2;
        comaccx(3) = comaccx(3)+ comacc(3,k_n_t_n+i)^2*dt;
        comxx(3) = comxx(3)+ com(3,k_n_t_n+i)^2; 
        comvx(3) = comvx(3)+ comv(3,k_n_t_n+i)^2;
        E(3) = E(3) + abs(comv(3,k_n_t_n+i)* comacc(3,k_n_t_n+i)*dt);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i);
        
        r_a(3,k_n_t_n+i) = (R_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(3,k_n_t_n+i) = (L_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footx(j))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            r_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;  
        else     
            l_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            l_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;             
        end
        Er(3)=Er(3)+abs(r_v(3,k_n_t_n+i)* r_acc(3,k_n_t_n+i)*dt);
        El(3)=El(3)+abs(l_v(3,k_n_t_n+i)* l_acc(3,k_n_t_n+i)*dt);          
    end

    %%%%求解中间阶段的质心最优加速度和对应的质心轨迹和ZMP轨迹
    d_ti_tf = H2*Phi+H1*F;
    B = [0;1];
    for i=n_t_i+1:n_t_f  
        E_A_ti = [1,0;t_f-i*dt,1];
        comacc(3,k_n_t_n+i) = B'* E_A_ti * G_inv * d_ti_tf;
        %%%这里使用v_t1 = v_t2+(t1-t2)*a_t2;
        comv(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)*dt + comv(3,k_n_t_n+i-1);
        com(3,k_n_t_n+i) = comv(3,k_n_t_n+i)*dt + com(3,k_n_t_n+i-1);
        zmp_eq(3,k_n_t_n+i) = com(3,k_n_t_n+i) - comacc(3,k_n_t_n+i)/(w^2);
        comaccx(3) = comaccx(3)+ comacc(3,k_n_t_n+i)^2*dt;
        comxx(3) = comxx(3)+ com(3,k_n_t_n+i)^2;
        comvx(3) = comvx(3)+ comv(3,k_n_t_n+i)^2;
        E(3) = E(3) + abs(comv(3,k_n_t_n+i)* comacc(3,k_n_t_n+i)*dt);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i);
        
        r_a(3,k_n_t_n+i) = (R_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(3,k_n_t_n+i) = (L_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footx(j))/2+Cxs;
        r_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
        l_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
        r_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;
        l_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;  
        Er(3)=Er(3)+abs(r_v(3,k_n_t_n+i)* r_acc(3,k_n_t_n+i)*dt);
        El(3)=El(3)+abs(l_v(3,k_n_t_n+i)* l_acc(3,k_n_t_n+i)*dt);          
    end
    for i=1:n_t_n
        zmp_eq(3,k_n_t_n+i) = zmp_eq(3,k_n_t_n+i) +footx(j);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i) + footx(j);
    end       
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================有双足相,抛物线式参考轨迹===========
    % =============================================
    %%%从中间到两边拓展:det_t 为过渡相（也就是双足相时间）
    det_t_n = round(det_t/dt);
    t_i = t_ref(j,1)-round(det_t_n/2)*dt; 
    n_t_i = round(t_i/dt);
    t_i = n_t_i*dt;
    t_f = t_i+det_t_n*dt;
    n_t_f = round(t_f/dt);    
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 10;
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sx(j)-a_n;    


    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = 0;     x_if = a_n;   x_velo_i = 0;                              %相对于当前坐标footx(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：（利用对称性），修正摆动腿足部轨迹为二次多项式
    xfoot_ii = 0;     xfoot_if = Sx(j);   xfoot_velo_ii = 0;               %相对于当前坐标footx(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1);
    
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end


    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-foota_i1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ex)*(MM/(MB+ML));

    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_i00 = ax(4,j); a_i10 = ax(3,j);  a_i20 = ax(2,j); 
    for i=1:n_t_i
        zmp(4,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footx(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(4,k_n_t_n+i) = footx(j);
            R_foot(4,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);               
        else
            R_foot(4,k_n_t_n+i) = footx(j);
            L_foot(4,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footx(j);              
        end        
    end
    
    %%%%========================================== 质心轨迹后半段
    % =============================================
    x_fi = -a_n;  x_ff = 0;     x_velo_f = 0;                              %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi =-Sx(j);xfoot_ff =0;       xfoot_velo_f = 0;                  %相对于全局坐标式footx(j+1)=footx(j)+Sx(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);              
    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数  
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ex)*(MM/(MB+ML))+Sx(j);  %%%加上Sx(j)，转化为相对 坐标式footx(j) 
    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_f00 = ax(4,j); a_f10 = ax(3,j);  a_f20 = ax(2,j); 
    for i=n_t_f:n_t_n    
        zmp(4,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footx(j)+Sx(j); 
        if rem(j,2)==0 %右脚支撑
            L_foot(4,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            R_foot(4,k_n_t_n+i) = footx(j)+Sx(j);           
        else
            R_foot(4,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footx(j)+Sx(j);
            L_foot(4,k_n_t_n+i) = footx(j)+Sx(j);   
        end               
    end     

    %%%% zmp轨迹过渡中间相
    for i=n_t_i+1:n_t_f-1
        %%% 3次多项式求解参考ZMP和质心轨迹。
        if rem(j,2)==0 %右脚支撑
            L_foot(4,k_n_t_n+i) = footx(j);
            R_foot(4,k_n_t_n+i) = footx(j)+Sx(j);           
        else
            R_foot(4,k_n_t_n+i) = footx(j);
            L_foot(4,k_n_t_n+i) = footx(j)+Sx(j);
        end

        %%%%ZMP轨迹初始和末端时刻位置和速度
        zmp_condi = [a_n+footx(j);-a_n+footx(j)+Sx(j);0;0];
        ttt_f = [(t_i)^3,(t_i)^2,(t_i)^1,1;(t_f)^3,(t_f)^2,t_f,1;3*(t_i)^2,2*(t_i),1,0;3*(t_f)^2,2*(t_f),1,0];    
        a_zmp = ttt_f\zmp_condi;
        tj = dt*i;
        zmp(4,k_n_t_n+i) = a_zmp(4) + a_zmp(3)*(tj) + a_zmp(2)*(tj)^2 + a_zmp(1)*(tj)^3;           
    end       
    
    %%%%%%%%%%======================================全局质心加速度最小=======================求解    
    a_i0 = a_i00; a_i1 = a_i10;  a_i2 = a_i20; 
    a_f0 = a_f00; a_f1 = a_f10;  a_f2 = a_f20; 
    xx_u_t0 = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_s_t0 = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_u_tn = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_tn = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);   
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    Ah = xx_u_t0 - xx_s_t0;
    Bh = xx_s_tn - xx_u_tn;
    % h_add1 = -(w^2)/4*[Ah*exp(-w*t_i);Bh*exp(w*t_f-w*t_fn)];
    h_add1 = 0;
    h_add2 = -w*[a_i2*(1-exp(-w*t_i));a_f2*(1-exp(w*t_f-w*t_fn))];
    h_add = h_add1 + h_add2 ;
    %%%%%%%%矩阵求解
    % Wpre1 = w^2/4*exp(-2*w*t_i);   
    % Wpost1 = w^2/4*exp(2*w*t_f-2*w*t_fn); 
    Wpre1 = 0;   
    Wpost1 = 0; 
    Wpre2 = w^3/8*(1-exp(-2*w*t_i));   
    Wpost2 = w^3/8*(1-exp(2*w*t_f-2*w*t_fn));   
    Wpre = Wpre1 + Wpre2 ;
    Wpost = Wpost1 + Wpost2;


    W = [Wpre,0;0,Wpost];
    G = [det_t^3/3,det_t^2/2;det_t^2/2,det_t]; 
    G_inv = inv(G);
    phi_u = [1/2;w/2];  phi_s = [1/2;-w/2];
    M = [1,det_t;0,1];
    H11 = -M*phi_s;  H1 = [H11,phi_u];
    H21 = -M*phi_u;  H2 = [H21,phi_s];

    %%%计算已知分量：
    %%%xx_s_t = a0 + (2*a2 - w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%xx_u_t = a0 + (2*a2 + w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%其中：xx_c_t = a0 + a1*t + a2*t^2 + (2*a2)/w^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    F = [xx_s_ti;xx_u_tf];
    A = W + H2'*G_inv*H2;
    h = (W - H2'*G_inv*H1)*F + h_add;
    Phi = A\h;

    x_u_ti = Phi(1);
    x_s_tf = Phi(2);


    %%%% 质心轨迹前半段；
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);    
    for i=1:n_t_i        
        com(4,k_n_t_n+i) = exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 + (2*a_i2)/w^2);
        comv(4,k_n_t_n+i) = w*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i1 + 2*a_i2*(dt*i));
        comacc(4,k_n_t_n+i) = w^2*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + 2*a_i2;
        zmp_eq(4,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2;
        comaccx(4) = comaccx(4)+ comacc(4,k_n_t_n+i)^2*dt;
        comxx(4) = comxx(4)+ com(4,k_n_t_n+i)^2;
        comvx(4) = comvx(4)+ comv(4,k_n_t_n+i)^2;
        E(4) = E(4) + abs(comv(4,k_n_t_n+i)* comacc(4,k_n_t_n+i)*dt);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i);

        r_a(4,k_n_t_n+i) = (R_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(4,k_n_t_n+i) = (L_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footx(j))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            r_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) +2*foota_i2)/2;
            l_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;  
        else     
            l_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            l_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) +2*foota_i2)/2;
            r_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;             
        end
        Er(4)=Er(4)+abs(r_v(4,k_n_t_n+i)* r_acc(4,k_n_t_n+i)*dt);
        El(4)=El(4)+abs(l_v(4,k_n_t_n+i)* l_acc(4,k_n_t_n+i)*dt);         
    end

    %%%% 质心轨迹后半段
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    for i=n_t_f+1:n_t_n
        com(4,k_n_t_n+i) = exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2+ (a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 + (2*a_f2)/w^2);
        comv(4,k_n_t_n+i) = -w*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + (a_f1 + 2*a_f2*(dt*i));
        comacc(4,k_n_t_n+i) = w^2*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + 2*a_f2;
        zmp_eq(4,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2;
        comaccx(4) = comaccx(4)+ comacc(4,k_n_t_n+i)^2*dt;
        comxx(4) = comxx(4)+ com(4,k_n_t_n+i)^2; 
        comvx(4) = comvx(4)+ comv(4,k_n_t_n+i)^2;
        E(4) = E(4) + abs(comv(4,k_n_t_n+i)* comacc(4,k_n_t_n+i)*dt);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i);

        r_a(4,k_n_t_n+i) = (R_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(4,k_n_t_n+i) = (L_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footx(j))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            r_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;  
        else     
            l_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            l_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;             
        end
        Er(4)=Er(4)+abs(r_v(4,k_n_t_n+i)* r_acc(4,k_n_t_n+i)*dt);
        El(4)=El(4)+abs(l_v(4,k_n_t_n+i)* l_acc(4,k_n_t_n+i)*dt);         
    end

    %%%%求解中间阶段的质心最优加速度和对应的质心轨迹和ZMP轨迹
    d_ti_tf = H2*Phi+H1*F;
    B = [0;1];
    for i=n_t_i+1:n_t_f  
        E_A_ti = [1,0;t_f-i*dt,1];
        comacc(4,k_n_t_n+i) = B'* E_A_ti * G_inv * d_ti_tf;
        %%%这里使用v_t1 = v_t2+(t1-t2)*a_t2;
        comv(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)*dt + comv(4,k_n_t_n+i-1);
        com(4,k_n_t_n+i) = comv(4,k_n_t_n+i)*dt + com(4,k_n_t_n+i-1);
        zmp_eq(4,k_n_t_n+i) = com(4,k_n_t_n+i) - comacc(4,k_n_t_n+i)/(w^2);
        comaccx(4) = comaccx(4)+ comacc(4,k_n_t_n+i)^2*dt;
        comxx(4) = comxx(4)+ com(4,k_n_t_n+i)^2;
        comvx(4) = comvx(4)+ comv(4,k_n_t_n+i)^2;
        E(4) = E(4) + abs(comv(4,k_n_t_n+i)* comacc(4,k_n_t_n+i)*dt);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i);
        
        r_a(4,k_n_t_n+i) = (R_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footx(j))/2+Cxs;
        l_a(4,k_n_t_n+i) = (L_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footx(j))/2+Cxs;
        r_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
        l_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
        r_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;
        l_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;  
        Er(4)=Er(4)+abs(r_v(4,k_n_t_n+i)* r_acc(4,k_n_t_n+i)*dt);
        El(4)=El(4)+abs(l_v(4,k_n_t_n+i)* l_acc(4,k_n_t_n+i)*dt);             
    end
    for i=1:n_t_n
        zmp_eq(4,k_n_t_n+i) = zmp_eq(4,k_n_t_n+i) +footx(j);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i) + footx(j);
    end       
    
    
n_t_h = round((n_t_i + n_t_f)/2)+k_n_t_n;
Exx(1) = Exx(1)+(comv(1,n_t_h))^2 - 1/2*(comv(1,k_n_t_n+1))^2 - 1/2*(comv(1,end))^2; 
Exx(2) = Exx(2)+(comv(2,n_t_h))^2 - 1/2*(comv(2,k_n_t_n+1))^2 - 1/2*(comv(2,end))^2; 
Exx(3) = Exx(3)+(comv(3,n_t_h))^2 - 1/2*(comv(3,k_n_t_n+1))^2 - 1/2*(comv(3,end))^2; 
Exx(4) = Exx(4)+(comv(4,n_t_h))^2 - 1/2*(comv(4,k_n_t_n+1))^2 - 1/2*(comv(4,end))^2;     
    
end

det_t_n = round(det_t/dt);
t_i = t_ref(1)-round(det_t_n/2)*dt; 
n_t_i = round(t_i/dt);
t_i = n_t_i*dt;
t_f = t_i+det_t_n*dt;
n_t_f = round(t_f/dt);

det_comx = com(2,n_t_n-n_t_f+1)-com(2,1);
com_i = -com(2,1:(n_t_n-n_t_f));
for i =1:n_t_n-n_t_f
    com_i(i) = -com(2,n_t_n-n_t_f+1-i)-(com(2,2)-com(2,1));
end
com_f =  com(2,end-n_t_i+1:end);
com1 = com(2,:);
com_f = com_f+det_comx;
comx = [com_i,com1,com_f];
%%%%前面要加上一个双足相。
comx = comx(1,(n_t_n-det_t_n+1:end));
xxx = comx(1,round(det_t_n/2)+1);
comx = comx-xxx;
com_start = zeros(1,n_t_n);
comx = [com_start,comx];
% comx = comx+3.3;
zmp=zmp_eq;
det_zmpx = zmp(2,n_t_n-n_t_f+1)-zmp(2,1);
zmp_i = -zmp(2,1:(n_t_n-n_t_f));
for i =1:n_t_n-n_t_f
    zmp_i(i) = -zmp(2,n_t_n-n_t_f+1-i)-(zmp(2,2)-zmp(2,1));
end
zmp_f =  zmp(2,end-n_t_i+1:end);
zmp1 = zmp(2,:);
zmp_f = zmp_f+det_zmpx;
zmpx = [zmp_i,zmp1,zmp_f];
%%%%前面要加上一个双足相。
zmpx = zmpx(1,(n_t_n-det_t_n+1:end));
xxx = zmpx(1,round(det_t_n/2)+1);
zmpx = zmpx-xxx;
zmp_start = zeros(1,n_t_n);
zmpx = [zmp_start,zmpx];
end
function [zmp,zmp_eq,com,comv,comacc,comaccx,comxx,comvx,E,Exx,footx,Sx,comy,zmpy]=nao_20180111_zmp_comy_three_mass(dt,stepwidthx,nxx,b_n,T_ref,det_t)
stepwidth=2*(stepwidthx-11);

zmp = zeros(4,10);zmp_eq = zeros(4,10);com = zeros(4,10);comv = zeros(4,10);comacc = zeros(4,10);
comaccx = zeros(4,1);comxx = zeros(4,1);comvx = zeros(4,1);
E = zeros(4,1);Exx = zeros(4,1);

%%%%摆动足和支撑足部运动轨迹；
r_a = zeros(4,1);
l_a = zeros(4,1);
r_v = zeros(4,1);
l_v = zeros(4,1);
r_acc = zeros(4,1);
l_acc = zeros(4,1);
Er=zeros(4,1);
El=zeros(4,1);
%%%%左、右腿部摆动腿足部轨迹；
R_foot = zeros(4,10);
L_foot = zeros(4,10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%三质点问题求解初始条件：
g = 9800; Zc = 310; 
MB =2.841; ML=1.138; MM=MB+2*ML;
fb = 45.2;
a = [0 0 0;0 0 0;0 0 -100;0 0 -202.9;0 0 -202.9];
b = [-15.49,-0.29,-5.15;1.38,-2.21,-53.73;4.53,-2.25,-49.36;4.53,-2.25,-49.36;25.42,-3.3,-32.39];
c = [0.141,0.390,0.301,0.134,0.172];
ss = (c(1)*(a(1,:)+b(1,:))+c(2)*(a(2,:)+b(2,:))+c(3)*(a(3,:)+b(3,:))+c(4)*(a(4,:)+b(4,:))+c(5)*(a(5,:)+b(5,:)))/(ML);
Cxb = 0; Cyb = 0; Czb = Zc;
ss=[0,0,0];

%%周期
t_ref = zeros(nxx,2);  t_ref(:,1) = T_ref/2*ones(nxx,1); t_ref(:,2) = T_ref*ones(nxx,1);  
   
%%%步长与步宽
Sx = 60*ones(1,nxx);
Sy = stepwidth*ones(1,nxx+1);

%%%%%足部支撑点交替改变
footx = zeros(1,nxx);
footy = zeros(1,nxx);
for j=1:nxx
    if j == 1
        footy(1) = -stepwidth/2;
    else
        footy(j) = footy(j-1) + Sy(j)*((-1)^(j));
    end
end

%%%%等效ZMP位移
ax = zeros(4,nxx);ay=zeros(4,nxx);
cx=zeros(2,nxx); cy=zeros(2,nxx);

%%%%%参考ZMP前半单足相终止点和后半双足相终止点
xfoot =zeros(1,2);


for j = 1:nxx
    %%%第一步是右脚开始支撑（右脚支撑前半单足相+过渡相+左脚支撑后半单足相）：左右支撑出现对称性，所需没有必要区分偶数还是奇数
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================无双足相：t_i ==t_f===========
    % =============================================
    t_i = t_ref(j,1); t_f = t_ref(j,1);t_fn = t_ref(j,2); 
    n_t_i = round(t_i/dt); n_t_f = round(t_f/dt);

    n_t_n = round(t_fn/dt);
    if j == 1
        k_n_t_n = 0;
    else
        k_n_t_n = k_n_t_n + n_t_n;
    end
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 0;
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sy(j)-a_n;    

    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = 0;     x_if = a_n;   x_velo_i = 0;                              %相对于当前坐标footy(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：2次多项式： 
    if j==1
        xfoot_ii = Sy(1); 
    else
        xfoot_ii = (Sy(j-1)*((-1)^(j+1))+Sy(j)*((-1)^(j+1)))/2;
    end    
    xfoot_if = Sy(j)*((-1)^(j+1));   xfoot_velo_ii = 0;               %相对于当前坐标footy(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1); foota_i3 = 0;
    
    %%%%%%等效ZMP轨迹公式求解=============
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else 
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end
    Ey = (MB*Cyb+ML*Cyt+ML*Cys)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=-foota_i3*ML/(2*MM)*(MM/(MB+ML));   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-(foota_i1-Czs/g*6*foota_i3)*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ey)*(MM/(MB+ML));
    B = [exp(-w*0),exp(w*0);exp(-w*t_i),exp(w*t_i)];

    COMx_is = 10*(-1)^(j);
    COMx_es = Sy(j)/2*(-1)^(j+1);
    
    comx_condi1=COMx_is-(ax(1,j)*(0)^3+ax(2,j)*(0)^2+ax(3,j)*(0)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*0+2*ax(2,j)));
    comx_condi2=COMx_es-(ax(1,j)*(t_i)^3+ax(2,j)*(t_i)^2+ax(3,j)*(t_i)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*t_i+2*ax(2,j)));
    cx(:,j) = B\[comx_condi1;comx_condi2];  %cy(:,i) = B\comy_condi;    
    for i=1:n_t_i
        zmp(1,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footy(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(1,k_n_t_n+i) = footy(j);
            R_foot(1,k_n_t_n+i) = foota_i3*(dt*i)^3+foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);        
        else
            R_foot(1,k_n_t_n+i) = footy(j);
            L_foot(1,k_n_t_n+i) = foota_i3*(dt*i)^3+foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);            
        end 
        tj = dt*i;
        com(1,k_n_t_n+i)=cx(1,j)*exp(-w*tj)+cx(2,j)*exp(w*tj)+ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*tj+2*ax(2,j))+footy(j);            
        comv(1,k_n_t_n+i) = -w*cx(1,j)*exp(-w*tj)+w*cx(2,j)*exp(w*tj)+3*ax(1,j)*(tj)^2+2*ax(2,j)*(tj)^1+ax(3,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j));
        comacc(1,k_n_t_n+i) = w^2*cx(1,j)*exp(-w*tj)+w^2*cx(2,j)*exp(w*tj)+6*ax(1,j)*(tj)^1+2*ax(2,j);
        zmp_eq(1,k_n_t_n+i) = ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+footy(j);
        comaccx(1) = comaccx(1)+ comacc(1,k_n_t_n+i)^2*dt;
        comxx(1) = comxx(1)+ com(1,k_n_t_n+i)^2;
        comvx(1) = comvx(1)+ comv(1,k_n_t_n+i)^2;
        E(1) = E(1) + abs(comv(1,k_n_t_n+i)* comacc(1,k_n_t_n+i)*dt);            

        r_a(1,k_n_t_n+i) = (R_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        l_a(1,k_n_t_n+i) = (L_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 3*foota_i3*(dt*i)^2+2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            r_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 6*foota_i3*(dt*i)^1+2*foota_i2)/2;
            l_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;  
        else     
            l_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 3*foota_i3*(dt*i)^2+2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            l_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 6*foota_i3*(dt*i)^1+2*foota_i2)/2;
            r_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;             
        end
        Er(1)=Er(1)+abs(r_v(1,k_n_t_n+i)* r_acc(1,k_n_t_n+i)*dt);
        El(1)=El(1)+abs(l_v(1,k_n_t_n+i)* l_acc(1,k_n_t_n+i)*dt);           
        
        
    end
    
    %%%%==========================================质心轨迹后半段
    % =============================================
    x_fi = -a_n;  x_ff = 0;     x_velo_f = 0;                              %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi = Sy(j)*(-1)^j;    xfoot_ff =(Sy(j)*((-1)^(j))+Sy(j+1)*((-1)^(j)))/2;       xfoot_velo_f = 0;   %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);       

    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数            
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ey = (MB*Cyb+ML*Cyt+ML*Cys)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ey)*(MM/(MB+ML))+Sy(j)*(-1)^(j+1);  %%%加上Sy(j)，转化为相对 坐标式footy(j) 

    B = [exp(-w*t_f),exp(w*t_f);exp(-w*t_fn),exp(w*t_fn)];
    COMx_is = Sy(j)/2*(-1)^(j+1);
    COMx_es = (Sy(j)+10)*(-1)^(j+1);
    comx_condi1=COMx_is-(ax(1,j)*(t_f)^3+ax(2,j)*(t_f)^2+ax(3,j)*(t_f)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*t_f+2*ax(2,j)));
    comx_condi2=COMx_es-(ax(1,j)*(t_fn)^3+ax(2,j)*(t_fn)^2+ax(3,j)*(t_fn)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*t_fn+2*ax(2,j)));
    cx(:,j) = B\[comx_condi1;comx_condi2];  %cy(:,i) = B\comy_condi;      

    for i=n_t_f+1:n_t_n      
        zmp(1,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footy(j)+Sy(j)*(-1)^(j+1); 
        if rem(j,2)==0 %右脚支撑
            L_foot(1,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            R_foot(1,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);          
        else
            R_foot(1,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            L_foot(1,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);   
        end        
        tj = dt*i;
        com(1,k_n_t_n+i)=cx(1,j)*exp(-w*tj)+cx(2,j)*exp(w*tj)+ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j)*tj+2*ax(2,j))+footy(j);            
        comv(1,k_n_t_n+i) = -w*cx(1,j)*exp(-w*tj)+w*cx(2,j)*exp(w*tj)+3*ax(1,j)*(tj)^2+2*ax(2,j)*(tj)^1+ax(3,j)+(Ez*MM)/(g*(MB+ML))*(6*ax(1,j));
        comacc(1,k_n_t_n+i) = w^2*cx(1,j)*exp(-w*tj)+w^2*cx(2,j)*exp(w*tj)+6*ax(1,j)*(tj)^1+2*ax(2,j);
        zmp_eq(1,k_n_t_n+i) = ax(1,j)*(tj)^3+ax(2,j)*(tj)^2+ax(3,j)*(tj)^1+ax(4,j)+footy(j);
        comaccx(1) = comaccx(1)+ comacc(1,k_n_t_n+i)^2*dt;
        comxx(1) = comxx(1)+ com(1,k_n_t_n+i)^2;
        comvx(1) = comvx(1)+ comv(1,k_n_t_n+i)^2;
        E(1) = E(1) + abs(comv(1,k_n_t_n+i)* comacc(1,k_n_t_n+i)*dt);   
        
        r_a(1,k_n_t_n+i) = (R_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        l_a(1,k_n_t_n+i) = (L_foot(1,k_n_t_n+i)+com(1,k_n_t_n+i))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            r_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;  
        else     
            l_v(1,k_n_t_n+i) = ( comv(1,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(1,k_n_t_n+i) = comv(1,k_n_t_n+i)/2;
            l_acc(1,k_n_t_n+i) = ( comacc(1,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(1,k_n_t_n+i) = comacc(1,k_n_t_n+i)/2;             
        end
        Er(1)=Er(1)+abs(r_v(1,k_n_t_n+i)* r_acc(1,k_n_t_n+i)*dt);
        El(1)=El(1)+abs(l_v(1,k_n_t_n+i)* l_acc(1,k_n_t_n+i)*dt);          
        
        
    end
    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================有双足相,点式参考轨迹===========
    % =============================================
    %%%从中间到两边拓展:det_t 为过渡相（也就是双足相时间）
    det_t_n = round(det_t/dt);
    t_i = t_ref(j,1)-round(det_t_n/2)*dt; 
    n_t_i = round(t_i/dt);
    t_i = n_t_i*dt;
    t_f = t_i+det_t_n*dt;
    n_t_f = round(t_f/dt);    
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 0;
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sy(j)-a_n;    


    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = 0;     x_if = a_n;   x_velo_i = 0;                              %相对于当前坐标footy(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：（利用对称性），修正摆动腿足部轨迹为二次多项式
    if j==1
        xfoot_ii = Sy(1); 
    else
        xfoot_ii = (Sy(j-1)*((-1)^(j+1))+Sy(j)*((-1)^(j+1)))/2;
    end    
    xfoot_if = Sy(j)*((-1)^(j+1));   xfoot_velo_ii = 0;               %相对于当前坐标footy(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1);
    
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end


    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-foota_i1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ex)*(MM/(MB+ML));

    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_i00 = ax(4,j); a_i10 = ax(3,j);  a_i20 = ax(2,j); 
    for i=1:n_t_i
        zmp(2,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footy(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(2,k_n_t_n+i) = footy(j);
            R_foot(2,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);               
        else
            R_foot(2,k_n_t_n+i) = footy(j);
            L_foot(2,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);              
        end        
    end
    
    %%%%========================================== 质心轨迹后半段
    % =============================================
    x_fi = -a_n;  x_ff = 0;     x_velo_f = 0;                              %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi = Sy(j)*(-1)^j;    xfoot_ff =(Sy(j)*((-1)^(j))+Sy(j+1)*((-1)^(j)))/2;     xfoot_velo_f = 0;                  %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);              
    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数  
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ex)*(MM/(MB+ML))+Sy(j)*(-1)^(j+1);  %%%加上Sy(j)，转化为相对 坐标式footy(j) 
    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_f00 = ax(4,j); a_f10 = ax(3,j);  a_f20 = ax(2,j); 
    for i=n_t_f:n_t_n    
        zmp(2,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footy(j)+Sy(j)*(-1)^(j+1); 
        if rem(j,2)==0 %右脚支撑
            L_foot(2,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            R_foot(2,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);           
        else
            R_foot(2,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            L_foot(2,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);   
        end               
    end     

    %%%% zmp轨迹过渡中间相
    for i=n_t_i+1:n_t_f-1
        %%% 3次多项式求解参考ZMP和质心轨迹。
        if rem(j,2)==0 %右脚支撑
            L_foot(2,k_n_t_n+i) = footy(j);
            R_foot(2,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);           
        else
            R_foot(2,k_n_t_n+i) = footy(j);
            L_foot(2,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);
        end

        %%%%ZMP轨迹初始和末端时刻位置和速度
        zmp_condi = [a_n+footy(j);-a_n+footy(j)+Sy(j)*(-1)^(j+1);0;0];
        ttt_f = [(t_i)^3,(t_i)^2,(t_i)^1,1;(t_f)^3,(t_f)^2,t_f,1;3*(t_i)^2,2*(t_i),1,0;3*(t_f)^2,2*(t_f),1,0];    
        a_zmp = ttt_f\zmp_condi;
        tj = dt*i;
        zmp(2,k_n_t_n+i) = a_zmp(4) + a_zmp(3)*(tj) + a_zmp(2)*(tj)^2 + a_zmp(1)*(tj)^3;           
    end       
    
    %%%%%%%%%%======================================全局质心加速度最小=======================求解    
    a_i0 = a_i00; a_i1 = a_i10;  a_i2 = a_i20; 
    a_f0 = a_f00; a_f1 = a_f10;  a_f2 = a_f20; 
    xx_u_t0 = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_s_t0 = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_u_tn = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_tn = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);   
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    Ah = xx_u_t0 - xx_s_t0;
    Bh = xx_s_tn - xx_u_tn;
    % h_add1 = -(w^2)/4*[Ah*exp(-w*t_i);Bh*exp(w*t_f-w*t_fn)];
    h_add1 = 0;
    h_add2 = -w*[a_i2*(1-exp(-w*t_i));a_f2*(1-exp(w*t_f-w*t_fn))];
    h_add = h_add1 + h_add2 ;
    %%%%%%%%矩阵求解
    % Wpre1 = w^2/4*exp(-2*w*t_i);   
    % Wpost1 = w^2/4*exp(2*w*t_f-2*w*t_fn); 
    Wpre1 = 0;   
    Wpost1 = 0; 
    Wpre2 = w^3/8*(1-exp(-2*w*t_i));   
    Wpost2 = w^3/8*(1-exp(2*w*t_f-2*w*t_fn));   
    Wpre = Wpre1 + Wpre2 ;
    Wpost = Wpost1 + Wpost2;


    W = [Wpre,0;0,Wpost];
    G = [det_t^3/3,det_t^2/2;det_t^2/2,det_t]; 
    G_inv = inv(G);
    phi_u = [1/2;w/2];  phi_s = [1/2;-w/2];
    M = [1,det_t;0,1];
    H11 = -M*phi_s;  H1 = [H11,phi_u];
    H21 = -M*phi_u;  H2 = [H21,phi_s];

    %%%计算已知分量：
    %%%xx_s_t = a0 + (2*a2 - w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%xx_u_t = a0 + (2*a2 + w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%其中：xx_c_t = a0 + a1*t + a2*t^2 + (2*a2)/w^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    F = [xx_s_ti;xx_u_tf];
    A = W + H2'*G_inv*H2;
    h = (W - H2'*G_inv*H1)*F + h_add;
    Phi = A\h;

    x_u_ti = Phi(1);
    x_s_tf = Phi(2);


    %%%% 质心轨迹前半段；
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);    
    for i=1:n_t_i        
        com(2,k_n_t_n+i) = exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 + (2*a_i2)/w^2);
        comv(2,k_n_t_n+i) = w*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i1 + 2*a_i2*(dt*i));
        comacc(2,k_n_t_n+i) = w^2*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + 2*a_i2;
        zmp_eq(2,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2;
        comaccx(2) = comaccx(2)+ comacc(2,k_n_t_n+i)^2*dt;
        comxx(2) = comxx(2)+ com(2,k_n_t_n+i)^2;
        comvx(2) = comvx(2)+ comv(2,k_n_t_n+i)^2;
        E(2) = E(2) + abs(comv(2,k_n_t_n+i)* comacc(2,k_n_t_n+i)*dt);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i);
        
        r_a(2,k_n_t_n+i) = (R_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(2,k_n_t_n+i) = (L_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footy(j))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            r_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) +2*foota_i2)/2;
            l_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;  
        else     
            l_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            l_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) +2*foota_i2)/2;
            r_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;             
        end
        Er(2)=Er(2)+abs(r_v(2,k_n_t_n+i)* r_acc(2,k_n_t_n+i)*dt);
        El(2)=El(2)+abs(l_v(2,k_n_t_n+i)* l_acc(2,k_n_t_n+i)*dt);         
        
    end

    %%%% 质心轨迹后半段
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    for i=n_t_f+1:n_t_n
        com(2,k_n_t_n+i) = exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2+ (a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 + (2*a_f2)/w^2);
        comv(2,k_n_t_n+i) = -w*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + (a_f1 + 2*a_f2*(dt*i));
        comacc(2,k_n_t_n+i) = w^2*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + 2*a_f2;
        zmp_eq(2,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2;
        comaccx(2) = comaccx(2)+ comacc(2,k_n_t_n+i)^2*dt;
        comxx(2) = comxx(2)+ com(2,k_n_t_n+i)^2; 
        comvx(2) = comvx(2)+ comv(2,k_n_t_n+i)^2;
        E(2) = E(2) + abs(comv(2,k_n_t_n+i)* comacc(2,k_n_t_n+i)*dt);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i);
        
        r_a(2,k_n_t_n+i) = (R_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(2,k_n_t_n+i) = (L_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footy(j))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            r_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;  
        else     
            l_v(2,k_n_t_n+i) = ( comv(2,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
            l_acc(2,k_n_t_n+i) = ( comacc(2,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;             
        end
        Er(2)=Er(2)+abs(r_v(2,k_n_t_n+i)* r_acc(2,k_n_t_n+i)*dt);
        El(2)=El(2)+abs(l_v(2,k_n_t_n+i)* l_acc(2,k_n_t_n+i)*dt);           
    end

    %%%%求解中间阶段的质心最优加速度和对应的质心轨迹和ZMP轨迹
    d_ti_tf = H2*Phi+H1*F;
    B = [0;1];
    for i=n_t_i+1:n_t_f  
        E_A_ti = [1,0;t_f-i*dt,1];
        comacc(2,k_n_t_n+i) = B'* E_A_ti * G_inv * d_ti_tf;
        %%%这里使用v_t1 = v_t2+(t1-t2)*a_t2;
        comv(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)*dt + comv(2,k_n_t_n+i-1);
        com(2,k_n_t_n+i) = comv(2,k_n_t_n+i)*dt + com(2,k_n_t_n+i-1);
        zmp_eq(2,k_n_t_n+i) = com(2,k_n_t_n+i) - comacc(2,k_n_t_n+i)/(w^2);
        comaccx(2) = comaccx(2)+ comacc(2,k_n_t_n+i)^2*dt;
        comxx(2) = comxx(2)+ com(2,k_n_t_n+i)^2;
        comvx(2) = comvx(2)+ comv(2,k_n_t_n+i)^2;
        E(2) = E(2) + abs(comv(2,k_n_t_n+i)* comacc(2,k_n_t_n+i)*dt);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i);
        
        r_a(2,k_n_t_n+i) = (R_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(2,k_n_t_n+i) = (L_foot(2,k_n_t_n+i)+com(2,k_n_t_n+i)+footy(j))/2+Cxs;
        r_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
        l_v(2,k_n_t_n+i) = comv(2,k_n_t_n+i)/2;
        r_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;
        l_acc(2,k_n_t_n+i) = comacc(2,k_n_t_n+i)/2;  
        Er(2)=Er(2)+abs(r_v(2,k_n_t_n+i)* r_acc(2,k_n_t_n+i)*dt);
        El(2)=El(2)+abs(l_v(2,k_n_t_n+i)* l_acc(2,k_n_t_n+i)*dt);           
    end
    for i=1:n_t_n
        zmp_eq(2,k_n_t_n+i) = zmp_eq(2,k_n_t_n+i) +footy(j);
        com(2,k_n_t_n+i) = com(2,k_n_t_n+i) + footy(j);
    end

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================有双足相,直线式参考轨迹===========
    % =============================================
    %%%从中间到两边拓展:det_t 为过渡相（也就是双足相时间）
    det_t_n = round(det_t/dt);
    t_i = t_ref(j,1)-round(det_t_n/2)*dt; 
    n_t_i = round(t_i/dt);
    t_i = n_t_i*dt;
    t_f = t_i+det_t_n*dt;
    n_t_f = round(t_f/dt);    
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 10*(-1)^(j+1);
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sy(j)-a_n;    


    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = a_n;     x_if = 0;   x_velo_i = (x_if-x_ii)/(t_i);              %相对于当前坐标footy(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：（利用对称性），修正摆动腿足部轨迹为二次多项式
    if j==1
        xfoot_ii = Sy(1); 
    else
        xfoot_ii = (Sy(j-1)*((-1)^(j+1))+Sy(j)*((-1)^(j+1)))/2;
    end    
    xfoot_if = Sy(j)*((-1)^(j+1));   xfoot_velo_ii = 0;               %相对于当前坐标footy(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1);
    
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end


    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-foota_i1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ex)*(MM/(MB+ML));

    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_i00 = ax(4,j); a_i10 = ax(3,j);  a_i20 = ax(2,j); 
    for i=1:n_t_i
        zmp(3,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footy(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(3,k_n_t_n+i) = footy(j);
            R_foot(3,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);               
        else
            R_foot(3,k_n_t_n+i) = footy(j);
            L_foot(3,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);              
        end        
    end
    
    %%%%========================================== 质心轨迹后半段
    % =============================================
    x_fi = 0;  x_ff = -a_n;     x_velo_f = (x_ff-x_fi)/(t_fn-t_f);         %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi = Sy(j)*(-1)^j;    xfoot_ff =(Sy(j)*((-1)^(j))+Sy(j+1)*((-1)^(j)))/2;     xfoot_velo_f = 0;                  %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);              
    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数  
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ex)*(MM/(MB+ML))+Sy(j)*(-1)^(j+1);  %%%加上Sy(j)，转化为相对 坐标式footy(j) 
    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_f00 = ax(4,j); a_f10 = ax(3,j);  a_f20 = ax(2,j); 
    for i=n_t_f:n_t_n    
        zmp(3,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footy(j)+Sy(j)*(-1)^(j+1); 
        if rem(j,2)==0 %右脚支撑
            L_foot(3,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            R_foot(3,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);           
        else
            R_foot(3,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            L_foot(3,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);   
        end               
    end     

    %%%% zmp轨迹过渡中间相
    for i=n_t_i+1:n_t_f-1
        %%% 3次多项式求解参考ZMP和质心轨迹。
        if rem(j,2)==0 %右脚支撑
            L_foot(3,k_n_t_n+i) = footy(j);
            R_foot(3,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);           
        else
            R_foot(3,k_n_t_n+i) = footy(j);
            L_foot(3,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);
        end

        %%%%ZMP轨迹初始和末端时刻位置和速度
        zmp_condi = [a_n+footy(j);-a_n+footy(j)+Sy(j)*(-1)^(j+1);0;0];
        ttt_f = [(t_i)^3,(t_i)^2,(t_i)^1,1;(t_f)^3,(t_f)^2,t_f,1;3*(t_i)^2,2*(t_i),1,0;3*(t_f)^2,2*(t_f),1,0];    
        a_zmp = ttt_f\zmp_condi;
        tj = dt*i;
        zmp(3,k_n_t_n+i) = a_zmp(4) + a_zmp(3)*(tj) + a_zmp(2)*(tj)^2 + a_zmp(1)*(tj)^3;           
    end       
    
    %%%%%%%%%%======================================全局质心加速度最小=======================求解    
    a_i0 = a_i00; a_i1 = a_i10;  a_i2 = a_i20; 
    a_f0 = a_f00; a_f1 = a_f10;  a_f2 = a_f20; 
    xx_u_t0 = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_s_t0 = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_u_tn = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_tn = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);   
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    Ah = xx_u_t0 - xx_s_t0;
    Bh = xx_s_tn - xx_u_tn;
    % h_add1 = -(w^2)/4*[Ah*exp(-w*t_i);Bh*exp(w*t_f-w*t_fn)];
    h_add1 = 0;
    h_add2 = -w*[a_i2*(1-exp(-w*t_i));a_f2*(1-exp(w*t_f-w*t_fn))];
    h_add = h_add1 + h_add2 ;
    %%%%%%%%矩阵求解
    % Wpre1 = w^2/4*exp(-2*w*t_i);   
    % Wpost1 = w^2/4*exp(2*w*t_f-2*w*t_fn); 
    Wpre1 = 0;   
    Wpost1 = 0; 
    Wpre2 = w^3/8*(1-exp(-2*w*t_i));   
    Wpost2 = w^3/8*(1-exp(2*w*t_f-2*w*t_fn));   
    Wpre = Wpre1 + Wpre2 ;
    Wpost = Wpost1 + Wpost2;


    W = [Wpre,0;0,Wpost];
    G = [det_t^3/3,det_t^2/2;det_t^2/2,det_t]; 
    G_inv = inv(G);
    phi_u = [1/2;w/2];  phi_s = [1/2;-w/2];
    M = [1,det_t;0,1];
    H11 = -M*phi_s;  H1 = [H11,phi_u];
    H21 = -M*phi_u;  H2 = [H21,phi_s];

    %%%计算已知分量：
    %%%xx_s_t = a0 + (2*a2 - w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%xx_u_t = a0 + (2*a2 + w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%其中：xx_c_t = a0 + a1*t + a2*t^2 + (2*a2)/w^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    F = [xx_s_ti;xx_u_tf];
    A = W + H2'*G_inv*H2;
    h = (W - H2'*G_inv*H1)*F + h_add;
    Phi = A\h;

    x_u_ti = Phi(1);
    x_s_tf = Phi(2);


    %%%% 质心轨迹前半段；
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);    
    for i=1:n_t_i        
        com(3,k_n_t_n+i) = exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 + (2*a_i2)/w^2);
        comv(3,k_n_t_n+i) = w*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i1 + 2*a_i2*(dt*i));
        comacc(3,k_n_t_n+i) = w^2*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + 2*a_i2;
        zmp_eq(3,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2;
        comaccx(3) = comaccx(3)+ comacc(3,k_n_t_n+i)^2*dt;
        comxx(3) = comxx(3)+ com(3,k_n_t_n+i)^2;
        comvx(3) = comvx(3)+ comv(3,k_n_t_n+i)^2;
        E(3) = E(3) + abs(comv(3,k_n_t_n+i)* comacc(3,k_n_t_n+i)*dt);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i);
        
        r_a(3,k_n_t_n+i) = (R_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(3,k_n_t_n+i) = (L_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footy(j))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            r_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) +2*foota_i2)/2;
            l_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;  
        else     
            l_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            l_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) +2*foota_i2)/2;
            r_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;             
        end
        Er(3)=Er(3)+abs(r_v(3,k_n_t_n+i)* r_acc(3,k_n_t_n+i)*dt);
        El(3)=El(3)+abs(l_v(3,k_n_t_n+i)* l_acc(3,k_n_t_n+i)*dt);          
    end

    %%%% 质心轨迹后半段
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    for i=n_t_f+1:n_t_n
        com(3,k_n_t_n+i) = exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2+ (a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 + (2*a_f2)/w^2);
        comv(3,k_n_t_n+i) = -w*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + (a_f1 + 2*a_f2*(dt*i));
        comacc(3,k_n_t_n+i) = w^2*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + 2*a_f2;
        zmp_eq(3,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2;
        comaccx(3) = comaccx(3)+ comacc(3,k_n_t_n+i)^2*dt;
        comxx(3) = comxx(3)+ com(3,k_n_t_n+i)^2; 
        comvx(3) = comvx(3)+ comv(3,k_n_t_n+i)^2;
        E(3) = E(3) + abs(comv(3,k_n_t_n+i)* comacc(3,k_n_t_n+i)*dt);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i);
        
        r_a(3,k_n_t_n+i) = (R_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(3,k_n_t_n+i) = (L_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footy(j))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            r_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;  
        else     
            l_v(3,k_n_t_n+i) = ( comv(3,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
            l_acc(3,k_n_t_n+i) = ( comacc(3,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;             
        end
        Er(3)=Er(3)+abs(r_v(3,k_n_t_n+i)* r_acc(3,k_n_t_n+i)*dt);
        El(3)=El(3)+abs(l_v(3,k_n_t_n+i)* l_acc(3,k_n_t_n+i)*dt);           
    end

    %%%%求解中间阶段的质心最优加速度和对应的质心轨迹和ZMP轨迹
    d_ti_tf = H2*Phi+H1*F;
    B = [0;1];
    for i=n_t_i+1:n_t_f  
        E_A_ti = [1,0;t_f-i*dt,1];
        comacc(3,k_n_t_n+i) = B'* E_A_ti * G_inv * d_ti_tf;
        %%%这里使用v_t1 = v_t2+(t1-t2)*a_t2;
        comv(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)*dt + comv(3,k_n_t_n+i-1);
        com(3,k_n_t_n+i) = comv(3,k_n_t_n+i)*dt + com(3,k_n_t_n+i-1);
        zmp_eq(3,k_n_t_n+i) = com(3,k_n_t_n+i) - comacc(3,k_n_t_n+i)/(w^2);
        comaccx(3) = comaccx(3)+ comacc(3,k_n_t_n+i)^2*dt;
        comxx(3) = comxx(3)+ com(3,k_n_t_n+i)^2;
        comvx(3) = comvx(3)+ comv(3,k_n_t_n+i)^2;
        E(3) = E(3) + abs(comv(3,k_n_t_n+i)* comacc(3,k_n_t_n+i)*dt);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i);
        
        r_a(3,k_n_t_n+i) = (R_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(3,k_n_t_n+i) = (L_foot(3,k_n_t_n+i)+com(3,k_n_t_n+i)+footy(j))/2+Cxs;
        r_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
        l_v(3,k_n_t_n+i) = comv(3,k_n_t_n+i)/2;
        r_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;
        l_acc(3,k_n_t_n+i) = comacc(3,k_n_t_n+i)/2;  
        Er(3)=Er(3)+abs(r_v(3,k_n_t_n+i)* r_acc(3,k_n_t_n+i)*dt);
        El(3)=El(3)+abs(l_v(3,k_n_t_n+i)* l_acc(3,k_n_t_n+i)*dt);           
        
    end
    for i=1:n_t_n
        zmp_eq(3,k_n_t_n+i) = zmp_eq(3,k_n_t_n+i) +footy(j);
        com(3,k_n_t_n+i) = com(3,k_n_t_n+i) + footy(j);
    end
    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ================有双足相,抛物线式参考轨迹===========
    % =============================================
    %%%从中间到两边拓展:det_t 为过渡相（也就是双足相时间）
    det_t_n = round(det_t/dt);
    t_i = t_ref(j,1)-round(det_t_n/2)*dt; 
    n_t_i = round(t_i/dt);
    t_i = n_t_i*dt;
    t_f = t_i+det_t_n*dt;
    n_t_f = round(t_f/dt);    
    
    %%% a_n为调整因子，正号为正斜率。
    a_n = 10*(-1)^(j+1);
    %%%% 参考ZMP前半单足相终止点和后半双足相终止点
    xfoot(1) = a_n; xfoot(2) = Sy(j)-a_n;    


    %%%%========================================== 质心轨迹前半段
    % =============================================
    %%%参考ZMP轨迹：最大二阶（抛物线形）轨迹(注意这里直接写成相对于支撑点的坐标形式)
    x_ii = a_n;     x_if = 0;   x_velo_i = 0;                              %相对于当前坐标footy(j);        
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(0),1,0];
    a_ti = ttt_i\[x_ii;x_if;x_velo_i];
    a_i0 = a_ti(3); a_i1 = a_ti(2);  a_i2 = a_ti(1); 
    
    %%%摆动腿与固定腿足部轨迹：（利用对称性），修正摆动腿足部轨迹为二次多项式
    if j==1
        xfoot_ii = Sy(1); 
    else
        xfoot_ii = (Sy(j-1)*((-1)^(j+1))+Sy(j)*((-1)^(j+1)))/2;
    end    
    xfoot_if = Sy(j)*((-1)^(j+1));   xfoot_velo_ii = 0;               %相对于当前坐标footy(j);     
    ttt_i = [(0)^2,(0)^1,1;(t_i)^2,t_i,1;2*(t_i),1,0];
    foota_ti = ttt_i\[xfoot_ii;xfoot_if;xfoot_velo_ii];
    foota_i0 = foota_ti(3); foota_i1 = foota_ti(2);  foota_i2 = foota_ti(1);
    
    if rem(j,2)==0   %%%左脚支撑
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys=ss(2);  Czs = Zc/2;              %%%右脚摆动腿质心调整参数              
    else
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数             
    end


    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_i2-foota_i2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_i1-foota_i1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_i0-(foota_i0-Czs/g*2*foota_i2)*ML/(2*MM)-Ex)*(MM/(MB+ML));

    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_i00 = ax(4,j); a_i10 = ax(3,j);  a_i20 = ax(2,j); 
    for i=1:n_t_i
        zmp(4,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 +footy(j);
        if rem(j,2)==0   %%%左脚支撑
            L_foot(4,k_n_t_n+i) = footy(j);
            R_foot(4,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);               
        else
            R_foot(4,k_n_t_n+i) = footy(j);
            L_foot(4,k_n_t_n+i) = foota_i2*(dt*i)^2+foota_i1*(dt*i)^1+foota_i0+footy(j);              
        end        
    end
    
    %%%%========================================== 质心轨迹后半段
    % =============================================
    x_fi =0;  x_ff =  -a_n;     x_velo_f = 0;         %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_fn),1,0];
    a_fi = ttt_f\[x_fi;x_ff;x_velo_f];
    a_f0 = a_fi(3) ; a_f1 = a_fi(2);  a_f2 = a_fi(1);                
    
    xfoot_fi = Sy(j)*(-1)^j;    xfoot_ff =(Sy(j)*((-1)^(j))+Sy(j+1)*((-1)^(j)))/2;     xfoot_velo_f = 0;                  %相对于全局坐标式footy(j+1)=footy(j)+Sy(j);         
    ttt_f = [(t_f)^2,(t_f)^1,1;(t_fn)^2,t_fn,1;2*(t_f),1,0];
    foota_fi = ttt_f\[xfoot_fi;xfoot_ff;xfoot_velo_f];
    foota_f0 = foota_fi(3) ; foota_f1 = foota_fi(2);  foota_f2 = foota_fi(1);              
    if rem(j,2)==0 %右脚支撑
        Cxt = ss(1); Cyt=ss(2); Czt = Zc/2;               %%%右脚支撑腿质心调整参数
        Cxs = ss(1); Cys=-ss(2);  Czs = Zc/2;             %%%左脚摆动腿质心调整参数  
    else
        Cxt = ss(1); Cyt=-ss(2); Czt = Zc/2;              %%%左脚支撑腿质心调整参数
        Cxs = ss(1); Cys= ss(2);  Czs = Zc/2;             %%%右脚摆动腿质心调整参数     
    end

    Ex = (MB*Cxb+ML*Cxt+ML*Cxs)/MM;
    Ez = (MB*Czb+ML*Czt/2+ML*Czs/2)/MM;

    %%%等效Zc/g：
    w=sqrt(g*(MB+ML)/(Ez*MM));  
    %%%%% 求解出等效的zmp轨迹
    ax(1,j)=0;   ax(2,j)= (a_f2-foota_f2*ML/(2*MM))*(MM/(MB+ML));  
    ax(3,j)=(a_f1-foota_f1*ML/(2*MM))*(MM/(MB+ML)); ax(4,j)=(a_f0-(foota_f0-Czs/g*2*foota_f2)*ML/(2*MM)-Ex)*(MM/(MB+ML))+Sy(j)*(-1)^(j+1);  %%%加上Sy(j)，转化为相对 坐标式footy(j) 
    %%%%%%%%%%%%***************
    %%%%%%%%%%%% 把变量赋值给等效ZMP（替代原来的ZMP值进行计算）
    %%%%%%%%%%%%***************
    a_f00 = ax(4,j); a_f10 = ax(3,j);  a_f20 = ax(2,j); 
    for i=n_t_f:n_t_n    
        zmp(4,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 +footy(j)+Sy(j)*(-1)^(j+1); 
        if rem(j,2)==0 %右脚支撑
            L_foot(4,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            R_foot(4,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);           
        else
            R_foot(4,k_n_t_n+i) = foota_f2*(dt*i)^2+foota_f1*(dt*i)^1+foota_f0+footy(j)+Sy(j)*(-1)^(j+1);
            L_foot(4,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);   
        end               
    end     

    %%%% zmp轨迹过渡中间相
    for i=n_t_i+1:n_t_f-1
        %%% 3次多项式求解参考ZMP和质心轨迹。
        if rem(j,2)==0 %右脚支撑
            L_foot(4,k_n_t_n+i) = footy(j);
            R_foot(4,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);           
        else
            R_foot(4,k_n_t_n+i) = footy(j);
            L_foot(4,k_n_t_n+i) = footy(j)+Sy(j)*(-1)^(j+1);
        end

        %%%%ZMP轨迹初始和末端时刻位置和速度
        zmp_condi = [a_n+footy(j);-a_n+footy(j)+Sy(j)*(-1)^(j+1);0;0];
        ttt_f = [(t_i)^3,(t_i)^2,(t_i)^1,1;(t_f)^3,(t_f)^2,t_f,1;3*(t_i)^2,2*(t_i),1,0;3*(t_f)^2,2*(t_f),1,0];    
        a_zmp = ttt_f\zmp_condi;
        tj = dt*i;
        zmp(4,k_n_t_n+i) = a_zmp(4) + a_zmp(3)*(tj) + a_zmp(2)*(tj)^2 + a_zmp(1)*(tj)^3;           
    end       
    
    %%%%%%%%%%======================================全局质心加速度最小=======================求解    
    a_i0 = a_i00; a_i1 = a_i10;  a_i2 = a_i20; 
    a_f0 = a_f00; a_f1 = a_f10;  a_f2 = a_f20; 
    xx_u_t0 = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_s_t0 = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*0))/w^2 + a_i1*0 + a_i2*0^2;
    xx_u_tn = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_tn = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_fn))/w^2 + a_f1*t_fn + a_f2*(t_fn)^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);   
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    Ah = xx_u_t0 - xx_s_t0;
    Bh = xx_s_tn - xx_u_tn;
    % h_add1 = -(w^2)/4*[Ah*exp(-w*t_i);Bh*exp(w*t_f-w*t_fn)];
    h_add1 = 0;
    h_add2 = -w*[a_i2*(1-exp(-w*t_i));a_f2*(1-exp(w*t_f-w*t_fn))];
    h_add = h_add1 + h_add2 ;
    %%%%%%%%矩阵求解
    % Wpre1 = w^2/4*exp(-2*w*t_i);   
    % Wpost1 = w^2/4*exp(2*w*t_f-2*w*t_fn); 
    Wpre1 = 0;   
    Wpost1 = 0; 
    Wpre2 = w^3/8*(1-exp(-2*w*t_i));   
    Wpost2 = w^3/8*(1-exp(2*w*t_f-2*w*t_fn));   
    Wpre = Wpre1 + Wpre2 ;
    Wpost = Wpost1 + Wpost2;


    W = [Wpre,0;0,Wpost];
    G = [det_t^3/3,det_t^2/2;det_t^2/2,det_t]; 
    G_inv = inv(G);
    phi_u = [1/2;w/2];  phi_s = [1/2;-w/2];
    M = [1,det_t;0,1];
    H11 = -M*phi_s;  H1 = [H11,phi_u];
    H21 = -M*phi_u;  H2 = [H21,phi_s];

    %%%计算已知分量：
    %%%xx_s_t = a0 + (2*a2 - w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%xx_u_t = a0 + (2*a2 + w*(a1 + 2*a2*t))/w^2 + a1*t + a2*t^2
    %%%其中：xx_c_t = a0 + a1*t + a2*t^2 + (2*a2)/w^2;
    xx_s_ti = a_i0 + (2*a_i2 - w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);
    xx_u_tf = a_f0 + (2*a_f2 + w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);

    F = [xx_s_ti;xx_u_tf];
    A = W + H2'*G_inv*H2;
    h = (W - H2'*G_inv*H1)*F + h_add;
    Phi = A\h;

    x_u_ti = Phi(1);
    x_s_tf = Phi(2);


    %%%% 质心轨迹前半段；
    xx_u_ti = a_i0 + (2*a_i2 + w*(a_i1 + 2*a_i2*t_i))/w^2 + a_i1*t_i + a_i2*(t_i^2);    
    for i=1:n_t_i        
        com(4,k_n_t_n+i) = exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2 + (2*a_i2)/w^2);
        comv(4,k_n_t_n+i) = w*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + (a_i1 + 2*a_i2*(dt*i));
        comacc(4,k_n_t_n+i) = w^2*exp(-w*(t_i - dt*i))*(x_u_ti-xx_u_ti)/2 + 2*a_i2;
        zmp_eq(4,k_n_t_n+i) = a_i0 + a_i1*(dt*i) + a_i2*(dt*i)^2;
        comaccx(4) = comaccx(4)+ comacc(4,k_n_t_n+i)^2*dt;
        comxx(4) = comxx(4)+ com(4,k_n_t_n+i)^2;
        comvx(4) = comvx(4)+ comv(4,k_n_t_n+i)^2;
        E(4) = E(4) + abs(comv(4,k_n_t_n+i)* comacc(4,k_n_t_n+i)*dt);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i);
        
        r_a(4,k_n_t_n+i) = (R_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(4,k_n_t_n+i) = (L_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footy(j))/2+Cxs;
        if rem(j,2)==0   %%%左脚支撑
            r_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            l_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            r_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) +2*foota_i2)/2;
            l_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;  
        else     
            l_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) +2*foota_i2*(dt*i)+foota_i1 )/2;
            r_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            l_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) +2*foota_i2)/2;
            r_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;             
        end
        Er(4)=Er(4)+abs(r_v(4,k_n_t_n+i)* r_acc(4,k_n_t_n+i)*dt);
        El(4)=El(4)+abs(l_v(4,k_n_t_n+i)* l_acc(4,k_n_t_n+i)*dt);        
    end

    %%%% 质心轨迹后半段
    xx_s_tf = a_f0 + (2*a_f2 - w*(a_f1 + 2*a_f2*t_f))/w^2 + a_f1*t_f + a_f2*(t_f^2);
    for i=n_t_f+1:n_t_n
        com(4,k_n_t_n+i) = exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2+ (a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2 + (2*a_f2)/w^2);
        comv(4,k_n_t_n+i) = -w*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + (a_f1 + 2*a_f2*(dt*i));
        comacc(4,k_n_t_n+i) = w^2*exp(-w*(dt*i-t_f))*(x_s_tf-xx_s_tf)/2 + 2*a_f2;
        zmp_eq(4,k_n_t_n+i) = a_f0 + a_f1*(dt*i) + a_f2*(dt*i)^2;
        comaccx(4) = comaccx(4)+ comacc(4,k_n_t_n+i)^2*dt;
        comxx(4) = comxx(4)+ com(4,k_n_t_n+i)^2; 
        comvx(4) = comvx(4)+ comv(4,k_n_t_n+i)^2;
        E(4) = E(4) + abs(comv(4,k_n_t_n+i)* comacc(4,k_n_t_n+i)*dt);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i);
        
        r_a(4,k_n_t_n+i) = (R_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(4,k_n_t_n+i) = (L_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footy(j))/2+Cxs;
        if rem(j,2)==1   %%%左脚支撑
            r_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            l_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            r_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) + 2*foota_f2)/2;
            l_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;  
        else     
            l_v(4,k_n_t_n+i) = ( comv(4,k_n_t_n+i) + 2*foota_f2*(dt*i)+foota_f1 )/2;
            r_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
            l_acc(4,k_n_t_n+i) = ( comacc(4,k_n_t_n+i) + 2*foota_f2)/2;
            r_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;             
        end
        Er(4)=Er(4)+abs(r_v(4,k_n_t_n+i)* r_acc(4,k_n_t_n+i)*dt);
        El(4)=El(4)+abs(l_v(4,k_n_t_n+i)* l_acc(4,k_n_t_n+i)*dt);         
    end

    %%%%求解中间阶段的质心最优加速度和对应的质心轨迹和ZMP轨迹
    d_ti_tf = H2*Phi+H1*F;
    B = [0;1];
    for i=n_t_i+1:n_t_f  
        E_A_ti = [1,0;t_f-i*dt,1];
        comacc(4,k_n_t_n+i) = B'* E_A_ti * G_inv * d_ti_tf;
        %%%这里使用v_t1 = v_t2+(t1-t2)*a_t2;
        comv(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)*dt + comv(4,k_n_t_n+i-1);
        com(4,k_n_t_n+i) = comv(4,k_n_t_n+i)*dt + com(4,k_n_t_n+i-1);
        zmp_eq(4,k_n_t_n+i) = com(4,k_n_t_n+i) - comacc(4,k_n_t_n+i)/(w^2);
        comaccx(4) = comaccx(4)+ comacc(4,k_n_t_n+i)^2*dt;
        comxx(4) = comxx(4)+ com(4,k_n_t_n+i)^2;
        comvx(4) = comvx(4)+ comv(4,k_n_t_n+i)^2;
        E(4) = E(4) + abs(comv(4,k_n_t_n+i)* comacc(4,k_n_t_n+i)*dt);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i);
        
        r_a(4,k_n_t_n+i) = (R_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footy(j))/2+Cxs;
        l_a(4,k_n_t_n+i) = (L_foot(4,k_n_t_n+i)+com(4,k_n_t_n+i)+footy(j))/2+Cxs;
        r_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
        l_v(4,k_n_t_n+i) = comv(4,k_n_t_n+i)/2;
        r_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;
        l_acc(4,k_n_t_n+i) = comacc(4,k_n_t_n+i)/2;  
        Er(4)=Er(4)+abs(r_v(4,k_n_t_n+i)* r_acc(4,k_n_t_n+i)*dt);
        El(4)=El(4)+abs(l_v(4,k_n_t_n+i)* l_acc(4,k_n_t_n+i)*dt);          
    end
    for i=1:n_t_n
        zmp_eq(4,k_n_t_n+i) = zmp_eq(4,k_n_t_n+i) +footy(j);
        com(4,k_n_t_n+i) = com(4,k_n_t_n+i) + footy(j);
    end
        
    
    
    
    
n_t_h = round((n_t_i + n_t_f)/2)+k_n_t_n;
Exx(1) = Exx(1)+(comv(1,n_t_h))^2 - 1/2*(comv(1,k_n_t_n+1))^2 - 1/2*(comv(1,end))^2; 
Exx(2) = Exx(2)+(comv(2,n_t_h))^2 - 1/2*(comv(2,k_n_t_n+1))^2 - 1/2*(comv(2,end))^2; 
Exx(3) = Exx(3)+(comv(3,n_t_h))^2 - 1/2*(comv(3,k_n_t_n+1))^2 - 1/2*(comv(3,end))^2; 
Exx(4) = Exx(4)+(comv(4,n_t_h))^2 - 1/2*(comv(4,k_n_t_n+1))^2 - 1/2*(comv(4,end))^2;     
    
end

det_comy = 0;
com_i = com(2,1:(n_t_n-n_t_f));
for i =1:n_t_n-n_t_f
    com_i(i) = com(2,n_t_n-n_t_f+1-i);
end

com_f =  com(2,end-n_t_i+1:end);
[a,b] = size(com);
for i = 1:n_t_i
    com_f(i) = com(2,b+1-i);
end
com1 = com(2,:)+det_comy;
com_f = com_f+2*det_comy;
comy = [com_i,com1,com_f];
%%%%左移动一个半周期
comy = comy(1,2*n_t_n-det_t_n+1:end);
com_start = zeros(1,n_t_n);
comy = [com_start,comy];

zmp=zmp_eq;
det_zmpy = 0;
zmp_i = zmp(2,1:(n_t_n-n_t_f));
for i =1:n_t_n-n_t_f
    zmp_i(i) = zmp(2,n_t_n-n_t_f+1-i);
end
zmp_f =  zmp(2,end-n_t_i+1:end);
[a,b] = size(zmp);
for i = 1:n_t_i
    zmp_f(i) = zmp(2,b+1-i);
end
zmp1 = zmp(2,:)+det_zmpy;
zmp_f = zmp_f+2*det_zmpy;
zmpy = [zmp_i,zmp1,zmp_f];
%%%%左移动一个半周期
zmpy = zmpy(1,2*n_t_n-det_t_n+1:end);
zmp_start = zeros(1,n_t_n);
zmpy = [zmp_start,zmpy];


for j = 1:nxx
    Sy(j) = stepwidth;
    if j == 1
        footy(j) = 0-stepwidth/2;
    else
        footy(j) = footy(j-1)+(-1)^j*Sy(j);
    end
    
end
end
%%%%%足部轨迹；
function [Lpx,LRx,Lvx,Lwx,Rpx,RRx,Rvx,Rwx,Lpp_ground,Rpp_ground]=FootpR_polynomial_mod2(n,n_jieduan,dt,steplength,stepwidth,T_ref,alpha)
%%%%左右脚的位置和角度值
global  rfootxyza; global  lfootxyza;

%%%%左右脚的位置和角度曲线系数
global rfoot_ax;global rfoot_ay;global rfoot_az;global rfoot_aa
global lfoot_ax;global lfoot_ay;global lfoot_az;global lfoot_aa

%%%发生周期调整时的左右脚位置和角度曲线系数；
global rfoot_axx ;  global rfoot_ayy ;    global rfoot_azz ;   global rfoot_aaa;
global lfoot_axx ;  global lfoot_ayy ;    global lfoot_azz ;   global lfoot_aaa ;   

global flag_foot
global Kt; global T1;
global xxx
global foot_updown;
global xxxt;

x_stair=steplength*ones(1,n);  %%步长 
x_stair(1)=steplength/2;
% x_stair(5)=60;x_stair(6)=30;x_stair(8)=10;
h_stair=0*ones(1,n);  %台阶高度
w_stair=stepwidth*ones(1,n); %步宽   
% w_stair(5)=120;w_stair(7)=80;

%%% 落脚点
xfootx=zeros(n,1);
yfooty=zeros(n,1);
yfooty(1)=55;
%%%%
w_stair(1)=105;
zfootz=zeros(n,1);
for xx=2:n
    xfootx(xx)=xfootx(xx-1)+ x_stair(xx-1);
    yfooty(xx)=yfooty(xx-1)+ (-1)^(xx-1)*w_stair(xx-1);
    zfootz(xx)=zfootz(xx-1)+h_stair(xx-1);
end

%周期时间
Tdd = alpha*T_ref;
n_td = round(Tdd/dt);
Tdd = n_td*dt;
Ts=(T_ref-Tdd)*ones(n,1);    Td=Tdd*ones(n,1);

rfootxyza=zeros(6,1);
lfootxyza=zeros(6,1);
rfoot_ax=zeros(10,n);   %步长多项式系数
rfoot_ay=zeros(10,n);   %步宽多项式系数
rfoot_az=zeros(10,n);   %步高多项式系数
rfoot_aa=zeros(10,n);
lfoot_ax=zeros(10,n);   %步长多项式系数
lfoot_ay=zeros(10,n);   %步宽多项式系数
lfoot_az=zeros(10,n);   %步高多项式系数
lfoot_aa=zeros(10,n);

T1 = Ts+Td;     
Kt=zeros(1,n);
for k=2:n
    Kt(k)=Kt(k-1)+round(T1(k-1)/dt);
end
N_sum = Kt(end)+round(T1(end)/dt);

flag_foot=zeros(N_sum,2);    
%%%%当做是是否改变周期的标志位
flag = zeros(1,N_sum);

x_a=0;
y_a=20;

Tini =6;  T =T_ref;  
xxx = round((n_jieduan-1)*T/dt);
nx = N_sum-xxx;

nxx = nx+round(Tini/dt)-2*round(T/dt);
Lpx=zeros(3,1,nxx); LRx=zeros(3,3,nxx); Lvx=zeros(3,1,nxx); Lwx=zeros(3,1,nxx);
Rpx=zeros(3,1,nxx); RRx=zeros(3,3,nxx); Rvx=zeros(3,1,nxx); Rwx=zeros(3,1,nxx);
La =zeros(nxx,1);
Ra =zeros(nxx,1);

%%%%起步阶段加上后1s都已经确定了
k_start = round(Tini/dt);
t_start = dt:dt:k_start*dt;
[Lp,LR,Lv,Lw,Rp,RR,Rv,Rw]=FootpR_start(t_start,xfootx,yfooty,Tini,dt,Ts,Td);
Lpx(:,:,1:k_start)=Lp; LRx(:,:,1:k_start)=LR; Lvx(:,:,1:k_start)=Lv; Lwx(:,:,1:k_start)=Lw;
Rpx(:,:,1:k_start)=Rp; RRx(:,:,1:k_start)=RR; Rvx(:,:,1:k_start)=Rv; Rwx(:,:,1:k_start)=Rw;


%R的每一列是一个时间点的R，需要用reshape(R(:,i),[3,3])来复原
%后半足长度。
b1=48.3;
%前半足长度。
b2=102.3;
%步长，注意其意义。
th1=0*pi/180;
%离地时，足与地面的夹角。
th2=0*pi/180;
%-----------------------------------规划出变步长的足部轨迹-----------------------------
%根据 footx1的维数确定左右足部的摆动周数，奇数为左脚位置（第一步要手动规划），偶数为右脚位置。
t1=0.001;  t2=0.001; h=15;
%注意足部落脚的序列中第三个点是左脚的第一步落脚点，第四个点是右脚的第二步落脚点   
txx_ini = Kt(3)+1;

xxxn = k_start+1-txx_ini;
for k = txx_ini:nx 
    
    %%当前时刻对应的周期
    b = find(k>Kt);
    i_zhouqi = b(end); 
    %%%上一周期截止时刻
    I = Kt(i_zhouqi);
    %%当前时刻在当前周期的位置
    t_des = (k-I)*dt;  
    i=i_zhouqi;  

    flag_foot(k,1) = i_zhouqi;
    
    flag_foot1=flag_foot(:,1);
    xxx=find(flag_foot1==i_zhouqi);
    xxn = length(xxx);
    an =0;
    for j=1:xxn
        kxxx=xxx(j);
        an=an+flag_foot(kxxx,2);
    end
    
    if (flag(k)==0)&&(an==0)
        flag_foot(k,2) = 0;
        t_plan=[0,Td(i,:),Td(i,:)+t1,Td(i,:)+Ts(i,:)/2,T1(i)-t2,T1(i)];    
        if rem(i,2)==0  %右脚
             Rfoot_plan=[[xfootx(i-2,:);yfooty(i-2,:);zfootz(i-2,:)], [xfootx(i-2,:);yfooty(i-2,:);zfootz(i-2,:)], [xfootx(i-2,:)+b2*(1-cos(th1));yfooty(i-2,:);zfootz(i-2,:)+b2*sin(th1)],...
                 [xfootx(i-1,:);(yfooty(i,:)+yfooty(i-2,:))/2;(zfootz(i,:)+zfootz(i-2,:))/2+h], [xfootx(i,:)-b1*(1-cos(th2));yfooty(i,:);zfootz(i,:)+b1*sin(th2)],...
                 [xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
             Rangle_plan=[0,0,th1,th2/2,-th2,0];

             Lfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
             Langle_plan=[0,0,0,0,0,0];         
        else
             Lfoot_plan=[[xfootx(i-2,:);yfooty(i-2,:);zfootz(i-2,:)], [xfootx(i-2,:);yfooty(i-2,:);zfootz(i-2,:)], [xfootx(i-2,:)+b2*(1-cos(th1));yfooty(i-2,:);zfootz(i-2,:)+b2*sin(th1)],...
                 [xfootx(i-1,:);(yfooty(i,:)+yfooty(i-2,:))/2;(zfootz(i,:)+zfootz(i-2,:))/2+h], [xfootx(i,:)-b1*(1-cos(th2));yfooty(i,:);zfootz(i,:)+b1*sin(th2)],...
                 [xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
             Langle_plan=[0,0,th1,th2/2,-th2,0];

             Rfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
             Rangle_plan=[0,0,0,0,0,0];                
        end

        %%%收尾各两个加速度和中间6个值
        A = [0,0,0,0,0,0,0,0,1,0;
            0,0,0,0,0,0,0,2,0,0;
            t_plan(1)^9,t_plan(1)^8,t_plan(1)^7,t_plan(1)^6,t_plan(1)^5,t_plan(1)^4,t_plan(1)^3,t_plan(1)^2,t_plan(1)^1,1;
            t_plan(2)^9,t_plan(2)^8,t_plan(2)^7,t_plan(2)^6,t_plan(2)^5,t_plan(2)^4,t_plan(2)^3,t_plan(2)^2,t_plan(2)^1,1;
            t_plan(3)^9,t_plan(3)^8,t_plan(3)^7,t_plan(3)^6,t_plan(3)^5,t_plan(3)^4,t_plan(3)^3,t_plan(3)^2,t_plan(3)^1,1;
            t_plan(4)^9,t_plan(4)^8,t_plan(4)^7,t_plan(4)^6,t_plan(4)^5,t_plan(4)^4,t_plan(4)^3,t_plan(4)^2,t_plan(4)^1,1;
            t_plan(5)^9,t_plan(5)^8,t_plan(5)^7,t_plan(5)^6,t_plan(5)^5,t_plan(5)^4,t_plan(5)^3,t_plan(5)^2,t_plan(5)^1,1;
            t_plan(6)^9,t_plan(6)^8,t_plan(6)^7,t_plan(6)^6,t_plan(6)^5,t_plan(6)^4,t_plan(6)^3,t_plan(6)^2,t_plan(6)^1,1;
            9*t_plan(6)^8,8*t_plan(6)^7,7*t_plan(6)^6,6*t_plan(6)^5,5*t_plan(6)^4,4*t_plan(6)^3,3*t_plan(6)^2,2*t_plan(6),1,0;
            72*t_plan(6)^7,56*t_plan(6)^6,42*t_plan(6)^5,30*t_plan(6)^4,20*t_plan(6)^3,12*t_plan(6)^2,6*t_plan(6),2,0,0;
            ];

        Arfootx =[0,0,Rfoot_plan(1,:),0,0]'; Arfooty =[0,0,Rfoot_plan(2,:),0,0]';Arfootz =[0,0,Rfoot_plan(3,:),0,0]';Arfoota =[0,0,Rangle_plan(1,:),0,0]';
        Alfootx =[0,0,Lfoot_plan(1,:),0,0]'; Alfooty =[0,0,Lfoot_plan(2,:),0,0]';Alfootz =[0,0,Lfoot_plan(3,:),0,0]';Alfoota =[0,0,Langle_plan(1,:),0,0]';
        rfoot_ax(:,i) = A\Arfootx;  rfoot_ay(:,i) = A\Arfooty;    rfoot_az(:,i) = A\Arfootz;   rfoot_aa(:,i) = A\Arfoota;
        lfoot_ax(:,i) = A\Alfootx;  lfoot_ay(:,i) = A\Alfooty;    lfoot_az(:,i) = A\Alfootz;   lfoot_aa(:,i) = A\Alfoota;   

        rfootx1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_ax(:,i);
        rfooty1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_ay(:,i);
        rfootz1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_az(:,i);
        rfoota1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_aa(:,i);

        lfootx1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_ax(:,i);
        lfooty1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_ay(:,i);
        lfootz1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_az(:,i);
        lfoota1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_aa(:,i);        

        %%%%%
        Rp1 = [rfootx1;rfooty1;rfootz1];
        
        Rp = Rp1;
        Rv = zeros(3,1);Rw = zeros(3,1);

        Lp1 = [lfootx1;lfooty1;lfootz1];
        
        Lp = Lp1;
        Lv = zeros(3,1);Lw = zeros(3,1);   

        Rfootz = pchip(t_plan,Rfoot_plan(3,:));Lfootz = pchip(t_plan,Lfoot_plan(3,:));
        Raaa = pchip(t_plan,Rangle_plan); Laaa = pchip(t_plan,Langle_plan);
%         Rfootz = spline(t_plan,Rfoot_plan(3,:));Lfootz = spline(t_plan,Lfoot_plan(3,:));
%         Raaa = spline(t_plan,Rangle_plan); Laaa = spline(t_plan,Langle_plan);        
        Rp(3)=ppval(Rfootz,t_des); Lp(3)=ppval(Lfootz,t_des);
        rfoota1=ppval(Raaa,t_des); lfoota1=ppval(Laaa,t_des);
        
        RR = Eularangle(rfoota1);
        LR = Eularangle(lfoota1);
        
        rfootxyza(1:3,k)=Rp; rfootxyza(5,k)=rfoota1; lfootxyza(1:3,k)=Lp; lfootxyza(5,k)=lfoota1;
        La(k)=lfoota1; Ra(k)=rfoota1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% 进入调整模式
    else
        flag_foot(k,2) = 1;        
        if flag(k)==1
            xxxt=k;
        end       
        xx_t = (xxxt-I)*dt;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %%%%%调整时刻发生在双组相时刻
        if xx_t<=Td(i)
            if flag(k)==1
                %%%%调整周期和步长
                %%%注意这里周期
                Td(i_zhouqi)=t_des;
                Td_n = round(Td(i_zhouqi)/dt);
                Td(i_zhouqi)=Td_n*dt;
                Ts(i_zhouqi)=Ts(i_zhouqi)/2;
                Ts_n = round(Ts(i_zhouqi)/dt);
                Ts(i_zhouqi)=Ts_n*dt;                
                Td(i_zhouqi+1)=Td(i_zhouqi+1)/2;
                Td_n = round(Td(i_zhouqi+1)/dt);
                Td(i_zhouqi+1)=Td_n*dt;                

                x_stair(i_zhouqi-1) = x_stair(i_zhouqi-1)+x_a;
                w_stair(i_zhouqi-1) = w_stair(i_zhouqi-1)+y_a;
                T1 = Ts+Td;                 
                %重现计算修改后的各个周期对应的采样点数
                for xx=2:n
                    Kt(xx)=Kt(xx-1)+round(T1(xx-1)/dt);
                end
                %%%%%%重新计算落脚点
                for xx=2:n
                    xfootx(xx)=xfootx(xx-1)+ x_stair(xx-1);
                    yfooty(xx)=yfooty(xx-1)+ (-1)^(xx-1)*w_stair(xx-1);
                    zfootz(xx)=zfootz(xx-1)+h_stair(xx-1);
                end                
            end
            
            if t_des== Td(i)
                %%%%当前时刻仍然是双足相
                rfootx1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_ax(:,i);
                rfooty1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_ay(:,i);
                rfootz1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_az(:,i);
                rfoota1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_aa(:,i);

                lfootx1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_ax(:,i);
                lfooty1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_ay(:,i);
                lfootz1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_az(:,i);
                lfoota1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_aa(:,i);            
                %%%%%
                Rp1 = [rfootx1;rfooty1;rfootz1];
                RR = Eularangle(rfoota1);
                Rv = zeros(3,1);Rw = zeros(3,1);

                Lp1 = [lfootx1;lfooty1;lfootz1];
                LR = Eularangle(lfoota1); 
                Lv = zeros(3,1);Lw = zeros(3,1);    
                
                rfoota1=Ra(k-1);
                lfoota1=La(k-1);
                Rp = Rp1;
                Lp = Lp1;
                RR = Eularangle(rfoota1);
                LR = Eularangle(lfoota1);                
                Rp(3)= 2*rfootxyza(3,k-1)-rfootxyza(3,k-2);
                Lp(3)= 2*lfootxyza(3,k-1)-lfootxyza(3,k-2);
                rfootxyza(1:3,k)=Rp; rfootxyza(5,k)=rfoota1; lfootxyza(1:3,k)=Lp; lfootxyza(5,k)=lfoota1;
                La(k)=lfoota1; Ra(k)=rfoota1;
            else
                t_yu=t_des-Td(i);
                t_plan=[0,t1,Ts(i,:)/2,Ts(i,:)-t2,Ts(i,:)];

                %%%初始条件
                rax= rfoot_ax(1,i)*(Td(i))^9+rfoot_ax(2,i)*(Td(i))^8+rfoot_ax(3,i)*(Td(i))^7+rfoot_ax(4,i)*(Td(i))^6+rfoot_ax(5,i)*(Td(i))^5+rfoot_ax(6,i)*(Td(i))^4+rfoot_ax(7,i)*(Td(i))^3+rfoot_ax(8,i)*(Td(i))^2+rfoot_ax(9,i)*(Td(i))+rfoot_ax(10,i);                
                ravx= 9*rfoot_ax(1,i)*(Td(i))^8+8*rfoot_ax(2,i)*(Td(i))^7+7*rfoot_ax(3,i)*(Td(i))^6+6*rfoot_ax(4,i)*(Td(i))^5+5*rfoot_ax(5,i)*(Td(i))^4+4*rfoot_ax(6,i)*(Td(i))^3+3*rfoot_ax(7,i)*(Td(i))^2+2*rfoot_ax(8,i)*(Td(i))^1+rfoot_ax(9,i);
                raaccx=72*rfoot_ax(1,i)*(Td(i))^7+56*rfoot_ax(2,i)*(Td(i))^6+42*rfoot_ax(3,i)*(Td(i))^5+30*rfoot_ax(4,i)*(Td(i))^4+20*rfoot_ax(5,i)*(Td(i))^3+12*rfoot_ax(6,i)*(Td(i))^2+6*rfoot_ax(7,i)*(Td(i))^1+2*rfoot_ax(8,i);

                ray= rfoot_ay(1,i)*(Td(i))^9+rfoot_ay(2,i)*(Td(i))^8+rfoot_ay(3,i)*(Td(i))^7+rfoot_ay(4,i)*(Td(i))^6+rfoot_ay(5,i)*(Td(i))^5+rfoot_ay(6,i)*(Td(i))^4+rfoot_ay(7,i)*(Td(i))^3+rfoot_ay(8,i)*(Td(i))^2+rfoot_ay(9,i)*(Td(i))+rfoot_ay(10,i); 
                ravy= 9*rfoot_ay(1,i)*(Td(i))^8+8*rfoot_ay(2,i)*(Td(i))^7+7*rfoot_ay(3,i)*(Td(i))^6+6*rfoot_ay(4,i)*(Td(i))^5+5*rfoot_ay(5,i)*(Td(i))^4+4*rfoot_ay(6,i)*(Td(i))^3+3*rfoot_ay(7,i)*(Td(i))^2+2*rfoot_ay(8,i)*(Td(i))^1+rfoot_ay(9,i);
                raaccy=72*rfoot_ay(1,i)*(Td(i))^7+56*rfoot_ay(2,i)*(Td(i))^6+42*rfoot_ay(3,i)*(Td(i))^5+30*rfoot_ay(4,i)*(Td(i))^4+20*rfoot_ay(5,i)*(Td(i))^3+12*rfoot_ay(6,i)*(Td(i))^2+6*rfoot_ay(7,i)*(Td(i))^1+2*rfoot_ay(8,i);

                raz= rfoot_az(1,i)*(Td(i))^9+rfoot_az(2,i)*(Td(i))^8+rfoot_az(3,i)*(Td(i))^7+rfoot_az(4,i)*(Td(i))^6+rfoot_az(5,i)*(Td(i))^5+rfoot_az(6,i)*(Td(i))^4+rfoot_az(7,i)*(Td(i))^3+rfoot_az(8,i)*(Td(i))^2+rfoot_az(9,i)*(Td(i))+rfoot_az(10,i); 
                ravz= 9*rfoot_az(1,i)*(Td(i))^8+8*rfoot_az(2,i)*(Td(i))^7+7*rfoot_az(3,i)*(Td(i))^6+6*rfoot_az(4,i)*(Td(i))^5+5*rfoot_az(5,i)*(Td(i))^4+4*rfoot_az(6,i)*(Td(i))^3+3*rfoot_az(7,i)*(Td(i))^2+2*rfoot_az(8,i)*(Td(i))^1+rfoot_az(9,i);
                raaccz=72*rfoot_az(1,i)*(Td(i))^7+56*rfoot_az(2,i)*(Td(i))^6+42*rfoot_az(3,i)*(Td(i))^5+30*rfoot_az(4,i)*(Td(i))^4+20*rfoot_az(5,i)*(Td(i))^3+12*rfoot_az(6,i)*(Td(i))^2+6*rfoot_az(7,i)*(Td(i))^1+2*rfoot_az(8,i);

                raa= rfoot_aa(1,i)*(Td(i))^9+rfoot_aa(2,i)*(Td(i))^8+rfoot_aa(3,i)*(Td(i))^7+rfoot_aa(4,i)*(Td(i))^6+rfoot_aa(5,i)*(Td(i))^5+rfoot_aa(6,i)*(Td(i))^4+rfoot_aa(7,i)*(Td(i))^3+rfoot_aa(8,i)*(Td(i))^2+rfoot_aa(9,i)*(Td(i))+rfoot_aa(10,i);                     
                rava= 9*rfoot_aa(1,i)*(Td(i))^8+8*rfoot_aa(2,i)*(Td(i))^7+7*rfoot_aa(3,i)*(Td(i))^6+6*rfoot_aa(4,i)*(Td(i))^5+5*rfoot_aa(5,i)*(Td(i))^4+4*rfoot_aa(6,i)*(Td(i))^3+3*rfoot_aa(7,i)*(Td(i))^2+2*rfoot_aa(8,i)*(Td(i))^1+rfoot_aa(9,i);
                raacca=72*rfoot_aa(1,i)*(Td(i))^7+56*rfoot_aa(2,i)*(Td(i))^6+42*rfoot_aa(3,i)*(Td(i))^5+30*rfoot_aa(4,i)*(Td(i))^4+20*rfoot_aa(5,i)*(Td(i))^3+12*rfoot_aa(6,i)*(Td(i))^2+6*rfoot_aa(7,i)*(Td(i))^1+2*rfoot_aa(8,i);

                lax= lfoot_ax(1,i)*(Td(i))^9+lfoot_ax(2,i)*(Td(i))^8+lfoot_ax(3,i)*(Td(i))^7+lfoot_ax(4,i)*(Td(i))^6+lfoot_ax(5,i)*(Td(i))^5+lfoot_ax(6,i)*(Td(i))^4+lfoot_ax(7,i)*(Td(i))^3+lfoot_ax(8,i)*(Td(i))^2+lfoot_ax(9,i)*(Td(i))+lfoot_ax(10,i);                                    
                lavx= 9*lfoot_ax(1,i)*(Td(i))^8+8*lfoot_ax(2,i)*(Td(i))^7+7*lfoot_ax(3,i)*(Td(i))^6+6*lfoot_ax(4,i)*(Td(i))^5+5*lfoot_ax(5,i)*(Td(i))^4+4*lfoot_ax(6,i)*(Td(i))^3+3*lfoot_ax(7,i)*(Td(i))^2+2*lfoot_ax(8,i)*(Td(i))^1+lfoot_ax(9,i);
                laaccx=72*lfoot_ax(1,i)*(Td(i))^7+56*lfoot_ax(2,i)*(Td(i))^6+42*lfoot_ax(3,i)*(Td(i))^5+30*lfoot_ax(4,i)*(Td(i))^4+20*lfoot_ax(5,i)*(Td(i))^3+12*lfoot_ax(6,i)*(Td(i))^2+6*lfoot_ax(7,i)*(Td(i))^1+2*lfoot_ax(8,i);

                lay= lfoot_ay(1,i)*(Td(i))^9+lfoot_ay(2,i)*(Td(i))^8+lfoot_ay(3,i)*(Td(i))^7+lfoot_ay(4,i)*(Td(i))^6+lfoot_ay(5,i)*(Td(i))^5+lfoot_ay(6,i)*(Td(i))^4+lfoot_ay(7,i)*(Td(i))^3+lfoot_ay(8,i)*(Td(i))^2+lfoot_ay(9,i)*(Td(i))+lfoot_ay(10,i); 
                lavy= 9*lfoot_ay(1,i)*(Td(i))^8+8*lfoot_ay(2,i)*(Td(i))^7+7*lfoot_ay(3,i)*(Td(i))^6+6*lfoot_ay(4,i)*(Td(i))^5+5*lfoot_ay(5,i)*(Td(i))^4+4*lfoot_ay(6,i)*(Td(i))^3+3*lfoot_ay(7,i)*(Td(i))^2+2*lfoot_ay(8,i)*(Td(i))^1+lfoot_ay(9,i);
                laaccy=72*lfoot_ay(1,i)*(Td(i))^7+56*lfoot_ay(2,i)*(Td(i))^6+42*lfoot_ay(3,i)*(Td(i))^5+30*lfoot_ay(4,i)*(Td(i))^4+20*lfoot_ay(5,i)*(Td(i))^3+12*lfoot_ay(6,i)*(Td(i))^2+6*lfoot_ay(7,i)*(Td(i))^1+2*lfoot_ay(8,i);

                laz= lfoot_az(1,i)*(Td(i))^9+lfoot_az(2,i)*(Td(i))^8+lfoot_az(3,i)*(Td(i))^7+lfoot_az(4,i)*(Td(i))^6+lfoot_az(5,i)*(Td(i))^5+lfoot_az(6,i)*(Td(i))^4+lfoot_az(7,i)*(Td(i))^3+lfoot_az(8,i)*(Td(i))^2+lfoot_az(9,i)*(Td(i))+lfoot_az(10,i); 
                lavz= 9*lfoot_az(1,i)*(Td(i))^8+8*lfoot_az(2,i)*(Td(i))^7+7*lfoot_az(3,i)*(Td(i))^6+6*lfoot_az(4,i)*(Td(i))^5+5*lfoot_az(5,i)*(Td(i))^4+4*lfoot_az(6,i)*(Td(i))^3+3*lfoot_az(7,i)*(Td(i))^2+2*lfoot_az(8,i)*(Td(i))^1+lfoot_az(9,i);
                laaccz=72*lfoot_az(1,i)*(Td(i))^7+56*lfoot_az(2,i)*(Td(i))^6+42*lfoot_az(3,i)*(Td(i))^5+30*lfoot_az(4,i)*(Td(i))^4+20*lfoot_az(5,i)*(Td(i))^3+12*lfoot_az(6,i)*(Td(i))^2+6*lfoot_az(7,i)*(Td(i))^1+2*lfoot_az(8,i);

                laa= lfoot_aa(1,i)*(Td(i))^9+lfoot_aa(2,i)*(Td(i))^8+lfoot_aa(3,i)*(Td(i))^7+lfoot_aa(4,i)*(Td(i))^6+lfoot_aa(5,i)*(Td(i))^5+lfoot_aa(6,i)*(Td(i))^4+lfoot_aa(7,i)*(Td(i))^3+lfoot_aa(8,i)*(Td(i))^2+lfoot_aa(9,i)*(Td(i))+lfoot_aa(10,i);                     
                lava= 9*lfoot_aa(1,i)*(Td(i))^8+8*lfoot_aa(2,i)*(Td(i))^7+7*lfoot_aa(3,i)*(Td(i))^6+6*lfoot_aa(4,i)*(Td(i))^5+5*lfoot_aa(5,i)*(Td(i))^4+4*lfoot_aa(6,i)*(Td(i))^3+3*lfoot_aa(7,i)*(Td(i))^2+2*lfoot_aa(8,i)*(Td(i))^1+lfoot_aa(9,i);
                laacca=72*lfoot_aa(1,i)*(Td(i))^7+56*lfoot_aa(2,i)*(Td(i))^6+42*lfoot_aa(3,i)*(Td(i))^5+30*lfoot_aa(4,i)*(Td(i))^4+20*lfoot_aa(5,i)*(Td(i))^3+12*lfoot_aa(6,i)*(Td(i))^2+6*lfoot_aa(7,i)*(Td(i))^1+2*lfoot_aa(8,i);
                
                
                raz=2*rfootxyza(3,xxxt-1)-rfootxyza(3,xxxt-2);
                laz= 2*lfootxyza(3,xxxt-1)-lfootxyza(3,xxxt-2);
%                 ravz=0;raaccz=0;
%                 lavz=0;laaccz=0;
                raa=Ra(xxxt);laa=La(xxxt);
                %%%%注意这里左右脚位置的初始值应该是由系数进行计算的。不然会导致不连续
                if rem(i,2)==0  %右脚
                     Rfoot_plan=[[rax;ray;raz], [rax+b2*(1-cos(th1));ray;raz+b2*sin(th1)],...
                         [(xfootx(i,:)+rax)/2;(yfooty(i,:)+ray)/2;h], [xfootx(i,:)-b1*(1-cos(th2));yfooty(i,:);zfootz(i,:)+b1*sin(th2)],...
                         [xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
                     Rangle_plan=[raa+0,raa+th1,raa+th2/2,raa-th2,raa+0];

                     Lfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                         [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                         [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
                     Langle_plan=[laa+0,laa+0,laa+0,laa+0,laa+0];         
                else
                     Lfoot_plan=[[lax;lay;laz], [lax+b2*sin(th1);lay;laz+b2*sin(th1)],...
                         [(xfootx(i,:)+lax)/2;(yfooty(i,:)+lay)/2;h], [xfootx(i,:)-b1*(1-cos(th2));yfooty(i,:);zfootz(i,:)+b1*sin(th2)],...
                         [xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
                     Langle_plan=[laa+0,laa+th1,laa+th2/2,laa-th2,laa+0];

                     Rfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                         [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                         [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
                     Rangle_plan=[raa+0,raa+0,raa+0,raa+0,raa+0];                
                end

                %%%收尾各两个加速度和中间6个值
                A = [0,0,0,0,0,0,0,1,0;
                    0,0,0,0,0,0,2,0,0;
                    t_plan(1)^8,t_plan(1)^7,t_plan(1)^6,t_plan(1)^5,t_plan(1)^4,t_plan(1)^3,t_plan(1)^2,t_plan(1)^1,1;
                    t_plan(2)^8,t_plan(2)^7,t_plan(2)^6,t_plan(2)^5,t_plan(2)^4,t_plan(2)^3,t_plan(2)^2,t_plan(2)^1,1;
                    t_plan(3)^8,t_plan(3)^7,t_plan(3)^6,t_plan(3)^5,t_plan(3)^4,t_plan(3)^3,t_plan(3)^2,t_plan(3)^1,1;
                    t_plan(4)^8,t_plan(4)^7,t_plan(4)^6,t_plan(4)^5,t_plan(4)^4,t_plan(4)^3,t_plan(4)^2,t_plan(4)^1,1;
                    t_plan(5)^8,t_plan(5)^7,t_plan(5)^6,t_plan(5)^5,t_plan(5)^4,t_plan(5)^3,t_plan(5)^2,t_plan(5)^1,1;
                    8*t_plan(5)^7,7*t_plan(5)^6,6*t_plan(5)^5,5*t_plan(5)^4,4*t_plan(5)^3,3*t_plan(5)^2,2*t_plan(5),1,0;
                    56*t_plan(5)^6,42*t_plan(5)^5,30*t_plan(5)^4,20*t_plan(5)^3,12*t_plan(5)^2,6*t_plan(5),2,0,0;
                    ];                    

                Arfootx =[ravx,raaccx,Rfoot_plan(1,:),0,0]'; Arfooty =[ravy,raaccy,Rfoot_plan(2,:),0,0]';Arfootz =[ravz,raaccz,Rfoot_plan(3,:),0,0]';Arfoota =[rava,raacca,Rangle_plan(1,:),0,0]';
                Alfootx =[lavx,laaccx,Lfoot_plan(1,:),0,0]'; Alfooty =[lavy,laaccy,Lfoot_plan(2,:),0,0]';Alfootz =[lavz,laaccz,Lfoot_plan(3,:),0,0]';Alfoota =[lava,laacca,Langle_plan(1,:),0,0]';
                rfoot_axx = A\Arfootx;  rfoot_ayy = A\Arfooty;    rfoot_azz = A\Arfootz;   rfoot_aaa = A\Arfoota;
                lfoot_axx = A\Alfootx;  lfoot_ayy = A\Alfooty;    lfoot_azz = A\Alfootz;   lfoot_aaa = A\Alfoota;   

                rfootx1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_axx;
                rfooty1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_ayy;
                rfootz1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_azz;
                rfoota1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_aaa;
                lfootx1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_axx;
                lfooty1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_ayy;
                lfootz1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_azz;
                lfoota1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_aaa;        

                %%%%%
                Rp1 = [rfootx1;rfooty1;rfootz1];
                Rp = Rp1;
                Rv = zeros(3,1);Rw = zeros(3,1);

                Lp1 = [lfootx1;lfooty1;lfootz1];
                Lp = Lp1;
                Lv = zeros(3,1);Lw = zeros(3,1);   

                
                Rfootz = pchip(t_plan,Rfoot_plan(3,:));Lfootz = pchip(t_plan,Lfoot_plan(3,:));
                Raaa = pchip(t_plan,Rangle_plan); Laaa = pchip(t_plan,Langle_plan);
                Rp(3)=ppval(Rfootz,t_yu); Lp(3)=ppval(Lfootz,t_yu);
                rfoota1=ppval(Raaa,t_yu); lfoota1=ppval(Laaa,t_yu);    
                RR = Eularangle(rfoota1);
                LR = Eularangle(lfoota1);                
                rfootxyza(1:3,k)=Rp; rfootxyza(5,k)=rfoota1; lfootxyza(1:3,k)=Lp; lfootxyza(5,k)=lfoota1;  
                La(k)=lfoota1; Ra(k)=rfoota1;
            end

            %%%%%%%%%%%%%%%%%%%%%
            %%%%%%发生在单足项时刻
        else
            %%%计算当前时刻相对于单组相的时刻；   
            k_yu = k-I-round(Td(i_zhouqi)/dt);
            k_t = k_yu*dt;
                
            %%%%当前时刻之进行不来的调整
            if (flag(k)==1) %%%%单足项调整的当前时刻
                
                if k_t<Ts(i_zhouqi)/2
                    foot_updown=1;           %%%上升阶段
                else
                    foot_updown=0;           %%%下降阶段
                end 
                Ts(i_zhouqi)=(k_t+Ts(i_zhouqi))/2;
                Ts_n = round(Ts(i_zhouqi)/dt);
                Ts(i_zhouqi)=Ts_n*dt;                 
                Td(i_zhouqi+1)=Td(i_zhouqi+1)/2;
                Td_n = round(Td(i_zhouqi+1)/dt);
                Td(i_zhouqi+1)=Td_n*dt;                  
                x_stair(i_zhouqi-1) = x_stair(i_zhouqi-1)+x_a;
                w_stair(i_zhouqi-1) = w_stair(i_zhouqi-1)+y_a;             
                T1 = Ts+Td;                 
                %重现计算修改后的各个周期对应的采样点数
                for xx=2:n
                    Kt(xx)=Kt(xx-1)+round(T1(xx-1)/dt);
                end       
                %%%%%%重新计算落脚点
                for xx=2:n
                    xfootx(xx)=xfootx(xx-1)+ x_stair(xx-1);
                    yfooty(xx)=yfooty(xx-1)+ (-1)^(xx-1)*w_stair(xx-1);
                    zfootz(xx)=zfootz(xx-1)+h_stair(xx-1);
                end                   
                %%%%%初始时刻仍然使用以前的结论
                rfootx1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_ax(:,i);
                rfooty1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_ay(:,i);
                rfootz1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_az(:,i);
                rfoota1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*rfoot_aa(:,i);

                lfootx1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_ax(:,i);
                lfooty1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_ay(:,i);
                lfootz1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_az(:,i);
                lfoota1 = [t_des^9,t_des^8,t_des^7,t_des^6,t_des^5,t_des^4,t_des^3,t_des^2,t_des^1,1]*lfoot_aa(:,i);            
                %%%%%
                Rp1 = [rfootx1;rfooty1;rfootz1];
                Rp = Rp1;
                Rv = zeros(3,1);Rw = zeros(3,1);

                Lp1 = [lfootx1;lfooty1;lfootz1];
                Lp = Lp1;                
                Lv = zeros(3,1);Lw = zeros(3,1);    
                
                rfoota1=Ra(k-1);
                lfoota1=La(k-1);

                RR = Eularangle(rfoota1);
                LR = Eularangle(lfoota1);                  
                Rp(3)= 2*rfootxyza(3,k-1)-rfootxyza(3,k-2);
                Lp(3)= 2*lfootxyza(3,k-1)-lfootxyza(3,k-2);
                rfootxyza(1:3,k)=Rp; rfootxyza(5,k)=rfoota1; lfootxyza(1:3,k)=Lp; lfootxyza(5,k)=lfoota1;
                La(k)=lfoota1; Ra(k)=rfoota1;
                                
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%状态机切换：    
            else
                txxt = (xxxt-I)*dt;
                t_ini = (xxxt-I-round(Td(i)/dt))*dt;
                t_yu = (k-xxxt)*dt;
                %%%%%%这里暂且将初始时刻和整个单组项时刻不作区分
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                %%%%上升阶段
                   if foot_updown==1              
                       
                        t_plan=[0,t1,(Ts(i,:)-t_ini)/2,Ts(i,:)-t2-t_ini,Ts(i,:)-t_ini];

                        %%%初始条件
                        rax= rfoot_ax(1,i)*(txxt)^9+rfoot_ax(2,i)*(txxt)^8+rfoot_ax(3,i)*(txxt)^7+rfoot_ax(4,i)*(txxt)^6+rfoot_ax(5,i)*(txxt)^5+rfoot_ax(6,i)*(txxt)^4+rfoot_ax(7,i)*(txxt)^3+rfoot_ax(8,i)*(txxt)^2+rfoot_ax(9,i)*(txxt)+rfoot_ax(10,i);                
                        ravx= 9*rfoot_ax(1,i)*(txxt)^8+8*rfoot_ax(2,i)*(txxt)^7+7*rfoot_ax(3,i)*(txxt)^6+6*rfoot_ax(4,i)*(txxt)^5+5*rfoot_ax(5,i)*(txxt)^4+4*rfoot_ax(6,i)*(txxt)^3+3*rfoot_ax(7,i)*(txxt)^2+2*rfoot_ax(8,i)*(txxt)^1+rfoot_ax(9,i);
                        raaccx=72*rfoot_ax(1,i)*(txxt)^7+56*rfoot_ax(2,i)*(txxt)^6+42*rfoot_ax(3,i)*(txxt)^5+30*rfoot_ax(4,i)*(txxt)^4+20*rfoot_ax(5,i)*(txxt)^3+12*rfoot_ax(6,i)*(txxt)^2+6*rfoot_ax(7,i)*(txxt)^1+2*rfoot_ax(8,i);

                        ray= rfoot_ay(1,i)*(txxt)^9+rfoot_ay(2,i)*(txxt)^8+rfoot_ay(3,i)*(txxt)^7+rfoot_ay(4,i)*(txxt)^6+rfoot_ay(5,i)*(txxt)^5+rfoot_ay(6,i)*(txxt)^4+rfoot_ay(7,i)*(txxt)^3+rfoot_ay(8,i)*(txxt)^2+rfoot_ay(9,i)*(txxt)+rfoot_ay(10,i); 
                        ravy= 9*rfoot_ay(1,i)*(txxt)^8+8*rfoot_ay(2,i)*(txxt)^7+7*rfoot_ay(3,i)*(txxt)^6+6*rfoot_ay(4,i)*(txxt)^5+5*rfoot_ay(5,i)*(txxt)^4+4*rfoot_ay(6,i)*(txxt)^3+3*rfoot_ay(7,i)*(txxt)^2+2*rfoot_ay(8,i)*(txxt)^1+rfoot_ay(9,i);
                        raaccy=72*rfoot_ay(1,i)*(txxt)^7+56*rfoot_ay(2,i)*(txxt)^6+42*rfoot_ay(3,i)*(txxt)^5+30*rfoot_ay(4,i)*(txxt)^4+20*rfoot_ay(5,i)*(txxt)^3+12*rfoot_ay(6,i)*(txxt)^2+6*rfoot_ay(7,i)*(txxt)^1+2*rfoot_ay(8,i);

                        raz= rfoot_az(1,i)*(txxt)^9+rfoot_az(2,i)*(txxt)^8+rfoot_az(3,i)*(txxt)^7+rfoot_az(4,i)*(txxt)^6+rfoot_az(5,i)*(txxt)^5+rfoot_az(6,i)*(txxt)^4+rfoot_az(7,i)*(txxt)^3+rfoot_az(8,i)*(txxt)^2+rfoot_az(9,i)*(txxt)+rfoot_az(10,i); 
                        ravz= 9*rfoot_az(1,i)*(txxt)^8+8*rfoot_az(2,i)*(txxt)^7+7*rfoot_az(3,i)*(txxt)^6+6*rfoot_az(4,i)*(txxt)^5+5*rfoot_az(5,i)*(txxt)^4+4*rfoot_az(6,i)*(txxt)^3+3*rfoot_az(7,i)*(txxt)^2+2*rfoot_az(8,i)*(txxt)^1+rfoot_az(9,i);
                        raaccz=72*rfoot_az(1,i)*(txxt)^7+56*rfoot_az(2,i)*(txxt)^6+42*rfoot_az(3,i)*(txxt)^5+30*rfoot_az(4,i)*(txxt)^4+20*rfoot_az(5,i)*(txxt)^3+12*rfoot_az(6,i)*(txxt)^2+6*rfoot_az(7,i)*(txxt)^1+2*rfoot_az(8,i);

                        raa= rfoot_aa(1,i)*(txxt)^9+rfoot_aa(2,i)*(txxt)^8+rfoot_aa(3,i)*(txxt)^7+rfoot_aa(4,i)*(txxt)^6+rfoot_aa(5,i)*(txxt)^5+rfoot_aa(6,i)*(txxt)^4+rfoot_aa(7,i)*(txxt)^3+rfoot_aa(8,i)*(txxt)^2+rfoot_aa(9,i)*(txxt)+rfoot_aa(10,i);                     
                        rava= 9*rfoot_aa(1,i)*(txxt)^8+8*rfoot_aa(2,i)*(txxt)^7+7*rfoot_aa(3,i)*(txxt)^6+6*rfoot_aa(4,i)*(txxt)^5+5*rfoot_aa(5,i)*(txxt)^4+4*rfoot_aa(6,i)*(txxt)^3+3*rfoot_aa(7,i)*(txxt)^2+2*rfoot_aa(8,i)*(txxt)^1+rfoot_aa(9,i);
                        raacca=72*rfoot_aa(1,i)*(txxt)^7+56*rfoot_aa(2,i)*(txxt)^6+42*rfoot_aa(3,i)*(txxt)^5+30*rfoot_aa(4,i)*(txxt)^4+20*rfoot_aa(5,i)*(txxt)^3+12*rfoot_aa(6,i)*(txxt)^2+6*rfoot_aa(7,i)*(txxt)^1+2*rfoot_aa(8,i);

                        lax= lfoot_ax(1,i)*(txxt)^9+lfoot_ax(2,i)*(txxt)^8+lfoot_ax(3,i)*(txxt)^7+lfoot_ax(4,i)*(txxt)^6+lfoot_ax(5,i)*(txxt)^5+lfoot_ax(6,i)*(txxt)^4+lfoot_ax(7,i)*(txxt)^3+lfoot_ax(8,i)*(txxt)^2+lfoot_ax(9,i)*(txxt)+lfoot_ax(10,i);                                    
                        lavx= 9*lfoot_ax(1,i)*(txxt)^8+8*lfoot_ax(2,i)*(txxt)^7+7*lfoot_ax(3,i)*(txxt)^6+6*lfoot_ax(4,i)*(txxt)^5+5*lfoot_ax(5,i)*(txxt)^4+4*lfoot_ax(6,i)*(txxt)^3+3*lfoot_ax(7,i)*(txxt)^2+2*lfoot_ax(8,i)*(txxt)^1+lfoot_ax(9,i);
                        laaccx=72*lfoot_ax(1,i)*(txxt)^7+56*lfoot_ax(2,i)*(txxt)^6+42*lfoot_ax(3,i)*(txxt)^5+30*lfoot_ax(4,i)*(txxt)^4+20*lfoot_ax(5,i)*(txxt)^3+12*lfoot_ax(6,i)*(txxt)^2+6*lfoot_ax(7,i)*(txxt)^1+2*lfoot_ax(8,i);

                        lay= lfoot_ay(1,i)*(txxt)^9+lfoot_ay(2,i)*(txxt)^8+lfoot_ay(3,i)*(txxt)^7+lfoot_ay(4,i)*(txxt)^6+lfoot_ay(5,i)*(txxt)^5+lfoot_ay(6,i)*(txxt)^4+lfoot_ay(7,i)*(txxt)^3+lfoot_ay(8,i)*(txxt)^2+lfoot_ay(9,i)*(txxt)+lfoot_ay(10,i); 
                        lavy= 9*lfoot_ay(1,i)*(txxt)^8+8*lfoot_ay(2,i)*(txxt)^7+7*lfoot_ay(3,i)*(txxt)^6+6*lfoot_ay(4,i)*(txxt)^5+5*lfoot_ay(5,i)*(txxt)^4+4*lfoot_ay(6,i)*(txxt)^3+3*lfoot_ay(7,i)*(txxt)^2+2*lfoot_ay(8,i)*(txxt)^1+lfoot_ay(9,i);
                        laaccy=72*lfoot_ay(1,i)*(txxt)^7+56*lfoot_ay(2,i)*(txxt)^6+42*lfoot_ay(3,i)*(txxt)^5+30*lfoot_ay(4,i)*(txxt)^4+20*lfoot_ay(5,i)*(txxt)^3+12*lfoot_ay(6,i)*(txxt)^2+6*lfoot_ay(7,i)*(txxt)^1+2*lfoot_ay(8,i);

                        laz= lfoot_az(1,i)*(txxt)^9+lfoot_az(2,i)*(txxt)^8+lfoot_az(3,i)*(txxt)^7+lfoot_az(4,i)*(txxt)^6+lfoot_az(5,i)*(txxt)^5+lfoot_az(6,i)*(txxt)^4+lfoot_az(7,i)*(txxt)^3+lfoot_az(8,i)*(txxt)^2+lfoot_az(9,i)*(txxt)+lfoot_az(10,i); 
                        lavz= 9*lfoot_az(1,i)*(txxt)^8+8*lfoot_az(2,i)*(txxt)^7+7*lfoot_az(3,i)*(txxt)^6+6*lfoot_az(4,i)*(txxt)^5+5*lfoot_az(5,i)*(txxt)^4+4*lfoot_az(6,i)*(txxt)^3+3*lfoot_az(7,i)*(txxt)^2+2*lfoot_az(8,i)*(txxt)^1+lfoot_az(9,i);
                        laaccz=72*lfoot_az(1,i)*(txxt)^7+56*lfoot_az(2,i)*(txxt)^6+42*lfoot_az(3,i)*(txxt)^5+30*lfoot_az(4,i)*(txxt)^4+20*lfoot_az(5,i)*(txxt)^3+12*lfoot_az(6,i)*(txxt)^2+6*lfoot_az(7,i)*(txxt)^1+2*lfoot_az(8,i);

                        laa= lfoot_aa(1,i)*(txxt)^9+lfoot_aa(2,i)*(txxt)^8+lfoot_aa(3,i)*(txxt)^7+lfoot_aa(4,i)*(txxt)^6+lfoot_aa(5,i)*(txxt)^5+lfoot_aa(6,i)*(txxt)^4+lfoot_aa(7,i)*(txxt)^3+lfoot_aa(8,i)*(txxt)^2+lfoot_aa(9,i)*(txxt)+lfoot_aa(10,i);                     
                        lava= 9*lfoot_aa(1,i)*(txxt)^8+8*lfoot_aa(2,i)*(txxt)^7+7*lfoot_aa(3,i)*(txxt)^6+6*lfoot_aa(4,i)*(txxt)^5+5*lfoot_aa(5,i)*(txxt)^4+4*lfoot_aa(6,i)*(txxt)^3+3*lfoot_aa(7,i)*(txxt)^2+2*lfoot_aa(8,i)*(txxt)^1+lfoot_aa(9,i);
                        laacca=72*lfoot_aa(1,i)*(txxt)^7+56*lfoot_aa(2,i)*(txxt)^6+42*lfoot_aa(3,i)*(txxt)^5+30*lfoot_aa(4,i)*(txxt)^4+20*lfoot_aa(5,i)*(txxt)^3+12*lfoot_aa(6,i)*(txxt)^2+6*lfoot_aa(7,i)*(txxt)^1+2*lfoot_aa(8,i);

                        raz=2*rfootxyza(3,xxxt-1)-rfootxyza(3,xxxt-2);
                        laz= 2*lfootxyza(3,xxxt-1)-lfootxyza(3,xxxt-2);
        %                 ravz=0;raaccz=0;
        %                 lavz=0;laaccz=0;
                        ravz=0;raaccz=0;
                        lavz=0;laaccz=0;  
                        raa=Ra(xxxt);laa=La(xxxt);
                        %%%%注意这里左右脚位置的初始值应该是由系数进行计算的。不然会导致不连续
                        if rem(i,2)==0  %右脚
                             Rfoot_plan=[[rax;ray;raz], [(xfootx(i,:)+2*rax)/3;(yfooty(i,:)+2*ray)/3;(raz+h)/2],...
                                 [(xfootx(i,:)+xfootx(i-1,:))/2;(yfooty(i,:)+ray)/2;h], [xfootx(i,:)-b1*(1-cos(th2));yfooty(i,:);zfootz(i,:)+b1*sin(th2)],...
                                 [xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
                             Rangle_plan=[raa,th1/2,th2/4,-th2,0];

                             Lfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
                             Langle_plan=[laa+0,laa+0,laa+0,laa+0,laa+0];         
                        else
                             Lfoot_plan=[[lax;lay;laz], [(xfootx(i,:)+2*rax)/3;(yfooty(i,:)+ray)/2;(laz+h)/2],...
                                 [(xfootx(i,:)+xfootx(i-1,:))/2;(yfooty(i,:)+ray)/2;h], [xfootx(i,:)-b1*(1-cos(th2));yfooty(i,:);zfootz(i,:)+b1*sin(th2)],...
                                 [xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
                             Langle_plan=[raa,th1/2,th2/4,-th2,0];

                             Rfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],...
                                 [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
                             Rangle_plan=[raa+0,raa+0,raa+0,raa+0,raa+0];                
                        end
                        %%%收尾各两个加速度和中间6个值
                        A = [0,0,0,0,0,0,0,1,0;
                            0,0,0,0,0,0,2,0,0;
                            t_plan(1)^8,t_plan(1)^7,t_plan(1)^6,t_plan(1)^5,t_plan(1)^4,t_plan(1)^3,t_plan(1)^2,t_plan(1)^1,1;
                            t_plan(2)^8,t_plan(2)^7,t_plan(2)^6,t_plan(2)^5,t_plan(2)^4,t_plan(2)^3,t_plan(2)^2,t_plan(2)^1,1;
                            t_plan(3)^8,t_plan(3)^7,t_plan(3)^6,t_plan(3)^5,t_plan(3)^4,t_plan(3)^3,t_plan(3)^2,t_plan(3)^1,1;
                            t_plan(4)^8,t_plan(4)^7,t_plan(4)^6,t_plan(4)^5,t_plan(4)^4,t_plan(4)^3,t_plan(4)^2,t_plan(4)^1,1;
                            t_plan(5)^8,t_plan(5)^7,t_plan(5)^6,t_plan(5)^5,t_plan(5)^4,t_plan(5)^3,t_plan(5)^2,t_plan(5)^1,1;
                            8*t_plan(5)^7,7*t_plan(5)^6,6*t_plan(5)^5,5*t_plan(5)^4,4*t_plan(5)^3,3*t_plan(5)^2,2*t_plan(5),1,0;
                            56*t_plan(5)^6,42*t_plan(5)^5,30*t_plan(5)^4,20*t_plan(5)^3,12*t_plan(5)^2,6*t_plan(5),2,0,0;
                            ];                    

                        Arfootx =[ravx,raaccx,Rfoot_plan(1,:),0,0]'; Arfooty =[ravy,raaccy,Rfoot_plan(2,:),0,0]';Arfootz =[ravz,raaccz,Rfoot_plan(3,:),0,0]';Arfoota =[rava,raacca,Rangle_plan(1,:),0,0]';
                        Alfootx =[lavx,laaccx,Lfoot_plan(1,:),0,0]'; Alfooty =[lavy,laaccy,Lfoot_plan(2,:),0,0]';Alfootz =[lavz,laaccz,Lfoot_plan(3,:),0,0]';Alfoota =[lava,laacca,Langle_plan(1,:),0,0]';
                        rfoot_axx = A\Arfootx;  rfoot_ayy = A\Arfooty;    rfoot_azz = A\Arfootz;   rfoot_aaa = A\Arfoota;
                        lfoot_axx = A\Alfootx;  lfoot_ayy = A\Alfooty;    lfoot_azz = A\Alfootz;   lfoot_aaa = A\Alfoota;   

                        rfootx1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_axx;
                        rfooty1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_ayy;
                        rfootz1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_azz;
                        rfoota1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_aaa;
                        lfootx1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_axx;
                        lfooty1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_ayy;
                        lfootz1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_azz;
                        lfoota1 = [t_yu^8,t_yu^7,t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_aaa;        

                        %%%%%
                        Rp1 = [rfootx1;rfooty1;rfootz1];
                        Rp = Rp1;
                        Rv = zeros(3,1);Rw = zeros(3,1);

                        Lp1 = [lfootx1;lfooty1;lfootz1];
                        Lp = Lp1;
                        Lv = zeros(3,1);Lw = zeros(3,1);   
                        
                        Rfootz = pchip(t_plan,Rfoot_plan(3,:));Lfootz = pchip(t_plan,Lfoot_plan(3,:));
                        Raaa = pchip(t_plan,Rangle_plan); Laaa = pchip(t_plan,Langle_plan);
                        Rp(3)=ppval(Rfootz,t_yu); Lp(3)=ppval(Lfootz,t_yu);
                        rfoota1=ppval(Raaa,t_yu); lfoota1=ppval(Laaa,t_yu);  
       
                        RR = Eularangle(rfoota1);
                        LR = Eularangle(lfoota1);                        
                        rfootxyza(1:3,k)=Rp; rfootxyza(5,k)=rfoota1; lfootxyza(1:3,k)=Lp; lfootxyza(5,k)=lfoota1;   
                        La(k)=lfoota1; Ra(k)=rfoota1;
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                   %%%%%%%%%%%下降阶段
                   else                          
                        t_plan=[0,(Ts(i,:)-t_ini)/2,Ts(i,:)-t_ini];
                        %%%初始条件
                        rax= rfoot_ax(1,i)*(txxt)^9+rfoot_ax(2,i)*(txxt)^8+rfoot_ax(3,i)*(txxt)^7+rfoot_ax(4,i)*(txxt)^6+rfoot_ax(5,i)*(txxt)^5+rfoot_ax(6,i)*(txxt)^4+rfoot_ax(7,i)*(txxt)^3+rfoot_ax(8,i)*(txxt)^2+rfoot_ax(9,i)*(txxt)+rfoot_ax(10,i);                
                        ravx= 9*rfoot_ax(1,i)*(txxt)^8+8*rfoot_ax(2,i)*(txxt)^7+7*rfoot_ax(3,i)*(txxt)^6+6*rfoot_ax(4,i)*(txxt)^5+5*rfoot_ax(5,i)*(txxt)^4+4*rfoot_ax(6,i)*(txxt)^3+3*rfoot_ax(7,i)*(txxt)^2+2*rfoot_ax(8,i)*(txxt)^1+rfoot_ax(9,i);
                        raaccx=72*rfoot_ax(1,i)*(txxt)^7+56*rfoot_ax(2,i)*(txxt)^6+42*rfoot_ax(3,i)*(txxt)^5+30*rfoot_ax(4,i)*(txxt)^4+20*rfoot_ax(5,i)*(txxt)^3+12*rfoot_ax(6,i)*(txxt)^2+6*rfoot_ax(7,i)*(txxt)^1+2*rfoot_ax(8,i);

                        ray= rfoot_ay(1,i)*(txxt)^9+rfoot_ay(2,i)*(txxt)^8+rfoot_ay(3,i)*(txxt)^7+rfoot_ay(4,i)*(txxt)^6+rfoot_ay(5,i)*(txxt)^5+rfoot_ay(6,i)*(txxt)^4+rfoot_ay(7,i)*(txxt)^3+rfoot_ay(8,i)*(txxt)^2+rfoot_ay(9,i)*(txxt)+rfoot_ay(10,i); 
                        ravy= 9*rfoot_ay(1,i)*(txxt)^8+8*rfoot_ay(2,i)*(txxt)^7+7*rfoot_ay(3,i)*(txxt)^6+6*rfoot_ay(4,i)*(txxt)^5+5*rfoot_ay(5,i)*(txxt)^4+4*rfoot_ay(6,i)*(txxt)^3+3*rfoot_ay(7,i)*(txxt)^2+2*rfoot_ay(8,i)*(txxt)^1+rfoot_ay(9,i);
                        raaccy=72*rfoot_ay(1,i)*(txxt)^7+56*rfoot_ay(2,i)*(txxt)^6+42*rfoot_ay(3,i)*(txxt)^5+30*rfoot_ay(4,i)*(txxt)^4+20*rfoot_ay(5,i)*(txxt)^3+12*rfoot_ay(6,i)*(txxt)^2+6*rfoot_ay(7,i)*(txxt)^1+2*rfoot_ay(8,i);

                        raz= rfoot_az(1,i)*(txxt)^9+rfoot_az(2,i)*(txxt)^8+rfoot_az(3,i)*(txxt)^7+rfoot_az(4,i)*(txxt)^6+rfoot_az(5,i)*(txxt)^5+rfoot_az(6,i)*(txxt)^4+rfoot_az(7,i)*(txxt)^3+rfoot_az(8,i)*(txxt)^2+rfoot_az(9,i)*(txxt)+rfoot_az(10,i); 
                        ravz= 9*rfoot_az(1,i)*(txxt)^8+8*rfoot_az(2,i)*(txxt)^7+7*rfoot_az(3,i)*(txxt)^6+6*rfoot_az(4,i)*(txxt)^5+5*rfoot_az(5,i)*(txxt)^4+4*rfoot_az(6,i)*(txxt)^3+3*rfoot_az(7,i)*(txxt)^2+2*rfoot_az(8,i)*(txxt)^1+rfoot_az(9,i);
                        raaccz=72*rfoot_az(1,i)*(txxt)^7+56*rfoot_az(2,i)*(txxt)^6+42*rfoot_az(3,i)*(txxt)^5+30*rfoot_az(4,i)*(txxt)^4+20*rfoot_az(5,i)*(txxt)^3+12*rfoot_az(6,i)*(txxt)^2+6*rfoot_az(7,i)*(txxt)^1+2*rfoot_az(8,i);

                        raa= rfoot_aa(1,i)*(txxt)^9+rfoot_aa(2,i)*(txxt)^8+rfoot_aa(3,i)*(txxt)^7+rfoot_aa(4,i)*(txxt)^6+rfoot_aa(5,i)*(txxt)^5+rfoot_aa(6,i)*(txxt)^4+rfoot_aa(7,i)*(txxt)^3+rfoot_aa(8,i)*(txxt)^2+rfoot_aa(9,i)*(txxt)+rfoot_aa(10,i);                     
                        rava= 9*rfoot_aa(1,i)*(txxt)^8+8*rfoot_aa(2,i)*(txxt)^7+7*rfoot_aa(3,i)*(txxt)^6+6*rfoot_aa(4,i)*(txxt)^5+5*rfoot_aa(5,i)*(txxt)^4+4*rfoot_aa(6,i)*(txxt)^3+3*rfoot_aa(7,i)*(txxt)^2+2*rfoot_aa(8,i)*(txxt)^1+rfoot_aa(9,i);
                        raacca=72*rfoot_aa(1,i)*(txxt)^7+56*rfoot_aa(2,i)*(txxt)^6+42*rfoot_aa(3,i)*(txxt)^5+30*rfoot_aa(4,i)*(txxt)^4+20*rfoot_aa(5,i)*(txxt)^3+12*rfoot_aa(6,i)*(txxt)^2+6*rfoot_aa(7,i)*(txxt)^1+2*rfoot_aa(8,i);

                        lax= lfoot_ax(1,i)*(txxt)^9+lfoot_ax(2,i)*(txxt)^8+lfoot_ax(3,i)*(txxt)^7+lfoot_ax(4,i)*(txxt)^6+lfoot_ax(5,i)*(txxt)^5+lfoot_ax(6,i)*(txxt)^4+lfoot_ax(7,i)*(txxt)^3+lfoot_ax(8,i)*(txxt)^2+lfoot_ax(9,i)*(txxt)+lfoot_ax(10,i);                                    
                        lavx= 9*lfoot_ax(1,i)*(txxt)^8+8*lfoot_ax(2,i)*(txxt)^7+7*lfoot_ax(3,i)*(txxt)^6+6*lfoot_ax(4,i)*(txxt)^5+5*lfoot_ax(5,i)*(txxt)^4+4*lfoot_ax(6,i)*(txxt)^3+3*lfoot_ax(7,i)*(txxt)^2+2*lfoot_ax(8,i)*(txxt)^1+lfoot_ax(9,i);
                        laaccx=72*lfoot_ax(1,i)*(txxt)^7+56*lfoot_ax(2,i)*(txxt)^6+42*lfoot_ax(3,i)*(txxt)^5+30*lfoot_ax(4,i)*(txxt)^4+20*lfoot_ax(5,i)*(txxt)^3+12*lfoot_ax(6,i)*(txxt)^2+6*lfoot_ax(7,i)*(txxt)^1+2*lfoot_ax(8,i);

                        lay= lfoot_ay(1,i)*(txxt)^9+lfoot_ay(2,i)*(txxt)^8+lfoot_ay(3,i)*(txxt)^7+lfoot_ay(4,i)*(txxt)^6+lfoot_ay(5,i)*(txxt)^5+lfoot_ay(6,i)*(txxt)^4+lfoot_ay(7,i)*(txxt)^3+lfoot_ay(8,i)*(txxt)^2+lfoot_ay(9,i)*(txxt)+lfoot_ay(10,i); 
                        lavy= 9*lfoot_ay(1,i)*(txxt)^8+8*lfoot_ay(2,i)*(txxt)^7+7*lfoot_ay(3,i)*(txxt)^6+6*lfoot_ay(4,i)*(txxt)^5+5*lfoot_ay(5,i)*(txxt)^4+4*lfoot_ay(6,i)*(txxt)^3+3*lfoot_ay(7,i)*(txxt)^2+2*lfoot_ay(8,i)*(txxt)^1+lfoot_ay(9,i);
                        laaccy=72*lfoot_ay(1,i)*(txxt)^7+56*lfoot_ay(2,i)*(txxt)^6+42*lfoot_ay(3,i)*(txxt)^5+30*lfoot_ay(4,i)*(txxt)^4+20*lfoot_ay(5,i)*(txxt)^3+12*lfoot_ay(6,i)*(txxt)^2+6*lfoot_ay(7,i)*(txxt)^1+2*lfoot_ay(8,i);

                        laz= lfoot_az(1,i)*(txxt)^9+lfoot_az(2,i)*(txxt)^8+lfoot_az(3,i)*(txxt)^7+lfoot_az(4,i)*(txxt)^6+lfoot_az(5,i)*(txxt)^5+lfoot_az(6,i)*(txxt)^4+lfoot_az(7,i)*(txxt)^3+lfoot_az(8,i)*(txxt)^2+lfoot_az(9,i)*(txxt)+lfoot_az(10,i); 
                        lavz= 9*lfoot_az(1,i)*(txxt)^8+8*lfoot_az(2,i)*(txxt)^7+7*lfoot_az(3,i)*(txxt)^6+6*lfoot_az(4,i)*(txxt)^5+5*lfoot_az(5,i)*(txxt)^4+4*lfoot_az(6,i)*(txxt)^3+3*lfoot_az(7,i)*(txxt)^2+2*lfoot_az(8,i)*(txxt)^1+lfoot_az(9,i);
                        laaccz=72*lfoot_az(1,i)*(txxt)^7+56*lfoot_az(2,i)*(txxt)^6+42*lfoot_az(3,i)*(txxt)^5+30*lfoot_az(4,i)*(txxt)^4+20*lfoot_az(5,i)*(txxt)^3+12*lfoot_az(6,i)*(txxt)^2+6*lfoot_az(7,i)*(txxt)^1+2*lfoot_az(8,i);

                        laa= lfoot_aa(1,i)*(txxt)^9+lfoot_aa(2,i)*(txxt)^8+lfoot_aa(3,i)*(txxt)^7+lfoot_aa(4,i)*(txxt)^6+lfoot_aa(5,i)*(txxt)^5+lfoot_aa(6,i)*(txxt)^4+lfoot_aa(7,i)*(txxt)^3+lfoot_aa(8,i)*(txxt)^2+lfoot_aa(9,i)*(txxt)+lfoot_aa(10,i);                     
                        lava= 9*lfoot_aa(1,i)*(txxt)^8+8*lfoot_aa(2,i)*(txxt)^7+7*lfoot_aa(3,i)*(txxt)^6+6*lfoot_aa(4,i)*(txxt)^5+5*lfoot_aa(5,i)*(txxt)^4+4*lfoot_aa(6,i)*(txxt)^3+3*lfoot_aa(7,i)*(txxt)^2+2*lfoot_aa(8,i)*(txxt)^1+lfoot_aa(9,i);
                        laacca=72*lfoot_aa(1,i)*(txxt)^7+56*lfoot_aa(2,i)*(txxt)^6+42*lfoot_aa(3,i)*(txxt)^5+30*lfoot_aa(4,i)*(txxt)^4+20*lfoot_aa(5,i)*(txxt)^3+12*lfoot_aa(6,i)*(txxt)^2+6*lfoot_aa(7,i)*(txxt)^1+2*lfoot_aa(8,i);

                        raz=2*rfootxyza(3,xxxt-1)-rfootxyza(3,xxxt-2);
                        laz= 2*lfootxyza(3,xxxt-1)-lfootxyza(3,xxxt-2);
        %                 ravz=0;raaccz=0;
        %                 lavz=0;laaccz=0;
                        ravz=0;raaccz=0;
                        lavz=0;laaccz=0;                           
                        raa=Ra(xxxt);laa=La(xxxt);
                        %%%%注意这里左右脚位置的初始值应该是由系数进行计算的。不然会导致不连续
                        if rem(i,2)==0  %右脚
                             Rfoot_plan=[[rax;ray;raz],[(xfootx(i,:)+rax)/2;(yfooty(i,:)+ray)/2;(zfootz(i,:)+raz)/2], [xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
                             Rangle_plan=[raa+0,-abs(raa/2),0];

                             Lfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
                             Langle_plan=[laa+0,laa+0,laa+0];         
                        else
                             Lfoot_plan=[[lax;lay;laz],[(xfootx(i,:)+lax)/2;(yfooty(i,:)+lay)/2;(zfootz(i,:)+laz)/2],[xfootx(i,:);yfooty(i,:);zfootz(i,:)]];
                             Langle_plan=[laa+0,-abs(raa/2),0];

                             Rfoot_plan=[[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)], [xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)],[xfootx(i-1,:);yfooty(i-1,:);zfootz(i-1,:)]];
                             Rangle_plan=[raa+0,raa+0,raa+0];                
                        end
                        
                        %%%收尾各两个加速度和中间3个值
                        A = [0,0,0,0,0,1,0;
                            0,0,0,0,2,0,0;
                            t_plan(1)^6,t_plan(1)^5,t_plan(1)^4,t_plan(1)^3,t_plan(1)^2,t_plan(1)^1,1;
                            t_plan(2)^6,t_plan(2)^5,t_plan(2)^4,t_plan(2)^3,t_plan(2)^2,t_plan(2)^1,1;
                            t_plan(3)^6,t_plan(3)^5,t_plan(3)^4,t_plan(3)^3,t_plan(3)^2,t_plan(3)^1,1;
                            6*t_plan(3)^5,5*t_plan(3)^4,4*t_plan(3)^3,3*t_plan(3)^2,2*t_plan(3),1,0;
                            30*t_plan(3)^4,20*t_plan(3)^3,12*t_plan(3)^2,6*t_plan(3),2,0,0;
                            ];                       
                        Arfootx =[ravx,raaccx,Rfoot_plan(1,:),0,0]'; Arfooty =[ravy,raaccy,Rfoot_plan(2,:),0,0]';Arfootz =[ravz,raaccz,Rfoot_plan(3,:),0,0]';Arfoota =[rava,raacca,Rangle_plan(1,:),0,0]';
                        Alfootx =[lavx,laaccx,Lfoot_plan(1,:),0,0]'; Alfooty =[lavy,laaccy,Lfoot_plan(2,:),0,0]';Alfootz =[lavz,laaccz,Lfoot_plan(3,:),0,0]';Alfoota =[lava,laacca,Langle_plan(1,:),0,0]';
                        rfoot_axx = A\Arfootx;  rfoot_ayy = A\Arfooty;    rfoot_azz = A\Arfootz;   rfoot_aaa = A\Arfoota;
                        lfoot_axx = A\Alfootx;  lfoot_ayy = A\Alfooty;    lfoot_azz = A\Alfootz;   lfoot_aaa = A\Alfoota;   

                        rfootx1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_axx;
                        rfooty1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_ayy;
                        rfootz1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_azz;
                        rfoota1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*rfoot_aaa;
                        lfootx1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_axx;
                        lfooty1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_ayy;
                        lfootz1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_azz;
                        lfoota1 = [t_yu^6,t_yu^5,t_yu^4,t_yu^3,t_yu^2,t_yu^1,1]*lfoot_aaa;        

                        %%%%%
                        Rp1 = [rfootx1;rfooty1;rfootz1];
                        
                        Rp = Rp1;
                        Rv = zeros(3,1);Rw = zeros(3,1);

                        Lp1 = [lfootx1;lfooty1;lfootz1];
                        
                        Lp = Lp1;
                        Lv = zeros(3,1);Lw = zeros(3,1);   
                        Rfootz = pchip(t_plan,Rfoot_plan(3,:));Lfootz = pchip(t_plan,Lfoot_plan(3,:));
                        Raaa = pchip(t_plan,Rangle_plan); Laaa = pchip(t_plan,Langle_plan);
                        Rp(3)=ppval(Rfootz,t_yu); Lp(3)=ppval(Lfootz,t_yu);
                        rfoota1=ppval(Raaa,t_yu); lfoota1=ppval(Laaa,t_yu);  
                        RR = Eularangle(rfoota1);
                        LR = Eularangle(lfoota1);
                        
                        rfootxyza(1:3,k)=Rp; rfootxyza(5,k)=rfoota1; lfootxyza(1:3,k)=Lp; lfootxyza(5,k)=lfoota1; 
                        La(k)=lfoota1; Ra(k)=rfoota1;
                   end                
            end
        end
    end
    
    Lpx(:,:,k+xxxn)=Lp; LRx(:,:,k+xxxn)=LR; Lvx(:,:,k+xxxn)=Lv; Lwx(:,:,k+xxxn)=Lw; 
    Rpx(:,:,k+xxxn)=Rp; RRx(:,:,k+xxxn)=RR; Rvx(:,:,k+xxxn)=Rv; Rwx(:,:,k+xxxn)=Rw;
end
Lpp_ground = Lpx;
Rpp_ground = Rpx;
for k = txx_ini:nx
    rr=[0,0,45.19]';
    Rpx(:,:,k+xxxn)=Rpx(:,:,k+xxxn)+RRx(:,:,k+xxxn)*rr;
    Lpx(:,:,k+xxxn)=Lpx(:,:,k+xxxn)+LRx(:,:,k+xxxn)*rr;
end

% Lpp=reshape(Lpx,3,nxx);
% Rpp=reshape(Rpx,3,nxx);
% kkk = 1:1:nxx;
% figure (2)
% hold on
% plot(kkk,Lpp(1,:),'r')
% plot(kkk,Rpp(1,:));
% legend('lfootx','rfootx');
% 
% figure (3)
% hold on
% plot(kkk,Lpp(2,:),'r')
% plot(kkk,Rpp(2,:));
% legend('lfooty','rfooty');
% 
% figure (4)
% hold on
% plot(kkk,Lpp(3,:),'r')
% plot(kkk,Rpp(3,:));
% legend('lfootz','rfootz');
% 
% figure (5)
% hold on
% plot(kkk,La(:,1),'r')
% plot(kkk,Ra(:,1));
% legend('lfoot-angle','rfoot-angle');

end
function [Lp,LR,Lv,Lw,Rp,RR,Rv,Rw]=FootpR_start(t,xfootx,yfooty,T_ini,dt,SSP,DSP)
%-------------parameter
%后半足长度。
b1=48.3;
%前半足长度。
b2=102.3;
T=DSP(1,1)+SSP(1,1);
%步长，注意其意义。
th1=0*pi/180;
%离地时，足与地面的夹角。
th2=0*pi/180;
h0=20;
tN=length(t);

%起步阶段，右脚迈出第一步
RLorRini=-55;
RFirststep=(RLorRini+yfooty(2,:))/2;
Rstep=yfooty(2,:);
LLorRini=55;
LFirststep=(LLorRini+yfooty(3,:))/2;
Lstep=yfooty(3,:);       
Rt_plan=[dt,T_ini-3*T,T_ini-2*T-2*DSP(1,1),T_ini-2*T,T_ini-T,T_ini-2*DSP(1,1),T_ini-DSP(1,1)/2];
Rfoot_plan=[[0;RLorRini;0],[0;RLorRini;0],[0;RLorRini;0],[b2*(1-cos(th2));RLorRini;b2*sin(th2)],[xfootx(2,:)/2;RFirststep;h0],...
[xfootx(2,:)-b1*(1-cos(th1));Rstep;b1*sin(th1)],[xfootx(2,:);Rstep;0]];
Rangle_plan=[0,0,0,th2,th2/2,th1,0];
%在 T_ini~T_ini+T之内左脚迈出第一步
Lt_plan=[dt,T_ini,T_ini+DSP(2,1),T_ini+DSP(2,1)+0.2,T_ini+DSP(2,1)+SSP(2,1)/2,T_ini+T-0.05,T_ini+T];
Lfoot_plan=[[0;LLorRini;0],[0;LLorRini;0],[0;LLorRini;0],[b2*(1-cos(th2));LLorRini;b2*sin(th2)],[xfootx(2,:);LFirststep;h0],...
[xfootx(3,:)-b1*(1-cos(th1));Lstep;b1*sin(th1)],[xfootx(3,:);Lstep;0]];
Langle_plan=[0,0,0,th2,th2/2,-th1,0]; 

Rfoot_pp=pchip(Rt_plan,Rfoot_plan);%至此，pp是3x1的矩阵，三行分别为gait p的x，y，z
R_aa=pchip(Rt_plan,Rangle_plan);
Lfoot_pp=pchip(Lt_plan,Lfoot_plan);
L_aa=pchip(Lt_plan,Langle_plan);
Rfoot=ppval(Rfoot_pp,t);
Lfoot=ppval(Lfoot_pp,t);
%右脚的参数 
Rq2w=ppval(R_aa,t);%采样角度
Rq2w=reshape(Rq2w,[1,1,tN]);
RR=Rodrigues([0 1 0],Rq2w);
Rpbm=reshape( Rfoot,[3,1,tN]);
Rp=Rpbm-sum(bsxfun(@times,RR,repmat([0 0 -45.19],[1,1,tN])),2);
Rv=Diff(t,Rp);
Rw=bsxfun(@times,repmat([0 1 0]',[1,1,tN]),Diff(t,Rq2w));%foot的w曲线,若有需要可取出fnder，这里fw是一个矩阵，列数为length(t)
Rw=reshape(Rw,[3,1,tN]); 
%左脚参数 
Lq2w=ppval(L_aa,t);%采样角度
Lq2w=reshape(Lq2w,[1,1,tN]);
LR=Rodrigues([0 1 0],Lq2w);
Lpbm=reshape( Lfoot,[3,1,tN]);
Lp=Lpbm-sum(bsxfun(@times,LR,repmat([0 0 -45.19],[1,1,tN])),2);
Lv=Diff(t,Lp);
Lw=bsxfun(@times,repmat([0 1 0]',[1,1,tN]),Diff(t,Lq2w));%foot的w曲线,若有需要可取出fnder，这里fw是一个矩阵，列数为length(t)
Lw=reshape(Lw,[3,1,tN]); 
function dy=Diff(x,y)
    [a,~,b]=size(y);
    y=reshape(y,[a,b]);
    cs=csapi(x,y);
    pp1=fnder(cs);
    dy=fnval(pp1,x);
    dy=reshape(dy,[a,1,b]);
end
end
function [ e_out ] = Rodrigues( a,q )
%把单位矢量旋转轴a变为帽子矩阵（斜对称矩阵）
%e_out是三维数组3x3xN，a是矢量，3x1；q是三维矢量，1x1xN
a_skew=zeros(3,3);
a_skew(1,2)=-a(3);
a_skew(1,3)=a(2);
a_skew(2,3)=-a(1);
a_skew=-a_skew'+a_skew;
%按旋转矩阵计算出Rodrigues式
%测试expm效率，隐去expm，自编expm
%e_out=expm(a_skew*q);%考虑用expm几
RN=length(q);%要求q是一个1x1xN的三维向量
%e_out=eye(3)+a_skew*sin(q)+a_skew^2*(1-cos(q));
%-------------向量化,本意是实现上式
a_skew1=repmat(a_skew,[1,1,RN]);
a_skew2=repmat(a_skew^2,[1,1,RN]);
%sq=reshape(sin(q)',[1,1,RN]);Initial把q弄成三维的向量了。
sq=sin(q);
%cq_1=reshape(1-cos(q)',[1,1,RN]);
cq_1=1-cos(q);
factor2=bsxfun(@times,a_skew1,sq);%sq行向量
factor3=bsxfun(@times,a_skew2,cq_1);%cq_1行向量
e_out=repmat(eye(3),[1,1,RN])+factor2+factor3;
end  
%%%%多次求解旋转矩阵
function [Rrpy]=Eularangle(angle)
%%%输入角度为q_ref=3*N矩阵；输出为3*3*N
q_ref=[0,angle,0];

Rrpy = [cos(q_ref(3)),  -sin(q_ref(3)),0;
        sin(q_ref(3)), cos(q_ref(3)), 0;
        0,               0,               1]*...
       [cos(q_ref(2)), 0,  sin(q_ref(2));
        0,               1,  0; 
       -sin(q_ref(2)), 0,  cos(q_ref(2))]*...
       [1,               0,  0;
        0, cos(q_ref(1)),  -sin(q_ref(1));
        0, sin(q_ref(1)),  cos(q_ref(1))];    

end

%%%%正逆运动学解算
function Forward_Kinematics(j,tN)
global Link;
%j=2，遍历赋值,global的值不会被释放？需要在command中给个global Link
%运动学正解需要躯干的初始化。
if j==0 
    return;
end
if j~=1
    i=Link(j).motherID;%下面的本意如注释，sum等是为了把b变为三维数组并与R相乘,注意b的转置问题。输入的第二项一定是行向量
    %因为，@times是对应元素相乘，多行对应单行。这样，得到的各个杆的R都是3x3xN个，p为3x1xN个
    %Link(i).p+Link(i).R*Link(j).b;
    Link(j).p=Link(i).p+sum(bsxfun(@times,Link(i).R,repmat(Link(j).b',[1,1,tN])),2);
    %Link(j).R=Link(i).R*Rodrigues(Link(j).a, Link(j).q);
    %注意这里的p和R已经是j相对于世界坐标系的位置和姿态，因为母杆坐标的p和R
    %是对世界坐标的，这相当于0Tj=0Ti * iTj
    %下面的语句（到end结束）是实现两个3x3xN的数组在第三维的相乘
    temR=repmat(Link(i).R,[],3);
    rot=Rodrigues(Link(j).a,Link(j).q);
    rot2=permute(rot,[2 1 3]);
    [m,~,~]=size(rot2);
    idx=1:m;
    idx=idx(ones(1,3),:);
    rot3=rot2(idx(:),:,:);
    tem=sum(bsxfun(@times,temR,rot3),2);
    Link(j).R=reshape(tem,[3,3,tN]);
end
Forward_Kinematics(Link(j).brotherID,tN);
Forward_Kinematics(Link(j).childID,tN);
end
function [J]=InverseKinematics(TargetID, posRef)
%TargetID和posRef分别是目标杆的ID和其参考位姿（包括p和R）
global Link;
%给各个关节角度初值,如果都是0则jacobian是奇异状态。
ray=Target_ray(TargetID);
N=16;%迭代次数上限
[~,~,tN]=size(Link(1).p);
delta_q=zeros(5,1,tN);
for n=1:N
    J=CalJacobian(ray);
    err=CalErr(TargetID,posRef);
%     -----------------------------------需要考虑
      if norm(err(:))<1E-3  %判断矩阵
          break;
      end
      %delta_q=0.5*J\err(:);%deltaq是一列向量，每N个是一个关节的所有时刻，长度是5N
      for i1=1:length(J(1,1,:))%取第三维的长度
          delta_q(:,:,i1)=J(:,:,i1)\err(:,:,i1)*0.5;
      end
      %delta_q=reshape(delta_q,[],tN);
      for nn=2:length(ray)
          j=ray(nn);
          Link(j).q=Link(j).q+delta_q(nn-1,:,:);
      end
       Forward_Kinematics(2,tN);
end
%----------Traget_ray（）给出从躯干到目标连杆的ID序列
    function ray=Target_ray(TargetID)
     i=Link(TargetID).motherID;
     if i==1
         ray=[1;TargetID];
     else
         ray=[Target_ray(i);TargetID];% ray = 1 2 3
     end
    end
%------CalJacobian()给出从躯干到目标连杆这N个杆之间的雅格比矩阵(6*N)
    function J=CalJacobian(ray)
        pn=Link(ray(end)).p;
        %------------向量化---------------------
        jN=length(ray);
        j2=zeros(6,1,tN);j3=zeros(6,1,tN);j4=zeros(6,1,tN);j5=zeros(6,1,tN);
        j6=zeros(6,1,tN);%eval出现的变量都要再eval之前存在。
        %   有关优化，1:）J是否应该欲加载一下2）j的最后一列，不用算pn-pn,直接写个0
        for i=2:jN   
            a=sum(bsxfun(@times,Link(ray(i)).R,repmat(Link(ray(i)).a',[1,1,tN])),2);%a要为行向量
            eval(['j' int2str(i) '=[cross(a,pn-Link(ray(i)).p);a];']);               %cross可以处理三维的。;也可以，看来还是给数组用的
        end
        J=[j2,j3,j4,j5,j6];%J  6x5xN,N表示时刻,6x5表示了5个关节角度。
    end
%------------计算位姿误差的函数，p和R(w)的误差
    function err=CalErr(TargetID,PosRef)
        dp=PosRef.p-Link(TargetID).p;
        %向量化下面一句
        %dw=Link(TargetID).R*R2W(Link(TargetID).R'*PosRef.R);
        forR=array3D_multi(permute(Link(TargetID).R,[2 1 3]),PosRef.R);
        tem1=permute(R2W(forR),[2 1 3]);
        dw=sum(bsxfun(@times,Link(TargetID).R,tem1),2);
        err=[dp;dw];%6x1xn
    end
%------------姿态到角速度的转换，服务于姿态误差函数。
    function w=R2W(dR)%判断两个矩阵是否相等，要设置阈值，不要用==
        %-------------向量化---------------
        tol=1e-004;
        E=repmat(eye(3),[1,1,tN]);
        temr=max(max(dR(:,:,:)-E(:,:,:)));
        ray0=find(temr(:)<tol);
        w(:,:,ray0)=zeros(3,1,length(ray0));%初始化并在1：tN中找出E，记录角标到ray0
        temray=1:tN;temray(ray0)=[];ray1=temray;%把1：tN的角标中非单位阵的角标放到ray1
        tem1=dR(1,1,ray1)+dR(2,2,ray1)+dR(3,3,ray1);
        thta=acos((tem1-1)/2);%1x1xN
        tem2=bsxfun(@rdivide,thta,2*sin(thta));
        tem3=[dR(3,2,ray1)-dR(2,3,ray1);dR(1,3,ray1)-dR(3,1,ray1);dR(2,1,ray1)-dR(1,2,ray1)];
        w(:,:,ray1)=bsxfun(@times,tem2,tem3);
%         if norm(dR-eye(3))<1e-004;
%             w=[0 0 0]';
%         else
%             thta=acos((sum(diag(dR))-1)/2);
%             w=thta/(2*sin(thta))*[dR(3,2)-dR(2,3);dR(1,3)-dR(3,1);dR(2,1)-dR(1,2)];
%         end
    end
    function R=array3D_multi(R1,R2)%两个3x3xN的矩阵在第三维相乘
        tR1=repmat(R1,[],3);
        tR2=permute(R2,[2 1 3]);
        [m,~,tN]=size(tR2);
        idx=1:m;
        idx=idx(ones(1,3),:);
        ttR2=tR2(idx(:),:,:);
        tem=sum(bsxfun(@times,tR1,ttR2),2);
        R=reshape(tem,[3,3,tN]);
    end
end
function [bp,bodyR,bodyv,bodyw]=trans_bp(bodyp,t,N)
bp=bodyp;
bodyR=reshape(repmat(eye(3),[1,N]),[3,3,N]);
bodyv=Diff(t,bodyp);
bodyw=kron([0 0 0]',ones(1,length(t)));
bodyw=reshape(bodyw,[3,1,N]);
function dy=Diff(x,y)
        if length(size(y))==3
            [a,~,b]=size(y);
            y=reshape(y,[a,b]);
            cs=csapi(x,y);
            pp1=fnder(cs);
            dy=fnval(pp1,x);
            dy=reshape(dy,[a,1,b]);
        else 
            cs=csapi(x,y);
            pp1=fnder(cs);
            dy=fnval(pp1,x);
        end
    end
end
function com=calCom(~)
%该程序返回总体质心，格式为3*N
global Link;
    function mc=calmc(j)
        if j==0
            mc=0;
        else
            %mc为各个杆在世界坐标系的质心位置，局部坐标1为腰部的一杆，其他坐标系包含下面的杆。
            %向量化下面一句
            %mc=Link(j).m*(Link(j).p+Link(j).R*Link(j).c);
            Link(j).mc=Link(j).p+sum(bsxfun(@times,Link(j).R,Link(j).c),2);
            mc=Link(j).m*Link(j).mc;
            mc=mc+calmc(Link(j).brotherID)+calmc(Link(j).childID);
        end
    end
    function M=TotalMass(j)
        if j==0
            M=0;
        else
            m=Link(j).m;
            M=m+TotalMass(Link(j).brotherID)+TotalMass(Link(j).childID);
        end
    end
Mc=calmc(1);%1->j
com=Mc./TotalMass(1);
com=reshape(com,3,[]);
end

%%%%多刚体条件下ZMP轨迹计算
function [ ZMPall ZMPl1 ] = cal_linkZMP(c,Rflink,Lflink,t)
%该程序返回总体ZMP和躯干ZMP，格式均为2*N，第一行为X坐标。
global Link N;
global N1;
N =N1;
pz=0;
g=9800;%这里单位是mm,kg,s,则g应取为mm/s^2或mN/kg
[P L]=calPL(Rflink,Lflink);
M=TotalMass(1);
%------------------------可考虑用Polyfit处理总体的P和L得到趋势项
dP=Diff(t,P);
dL=Diff(t,L);
allpx = (M*g*c(1,:)+pz*dP(1,:)-dL(2,:))./(M*g+dP(3,:));
allpy = (M*g*c(2,:)+pz*dP(2,:)+dL(1,:))./(M*g+dP(3,:));
ZMPall=[allpx;allpy];
%-----------------------总体ZMP求解结束-----------------------
%-----------------------躯干ZMP-------------------------------
dPl1=Diff(t,reshape(Link(1).P,3,[]));
dLl1=Diff(t,reshape(Link(1).L,3,[]));
m=Link(1).m;
cc=reshape(Link(1).mc,3,[]);
zmpl1x = (m*g*cc(1,:)+pz*dPl1(1,:)-dLl1(2,:))./(m*g+dPl1(3,:));
zmpl1y = (m*g*cc(2,:)+pz*dPl1(2,:)+dLl1(1,:))./(m*g+dPl1(3,:));
ZMPl1=[zmpl1x;zmpl1y];
%------------------------躯干ZMP结束--------------------------
    function [P L]=calPL(Rflink,Lflink)%计算各杆的PL，并返回总的PL(3*N的矩阵)
        %给出左右脚为目标连杆的链的dq，这样才能调用速度函数（计算所有杆的速度）
        Caldq(6,Rflink);%右TargetID_R
        Caldq(11,Lflink);%左TargetID_L
        ForwardVelocity(2);%有了各个杆的速度，可以计算L和P
        P=reshape(calP(1),3,[]);
        L=reshape(calL(1),3,[]);
        function ForwardVelocity(j)
        if j==0
            return;
        end
        if j~=1
            i=Link(j).motherID;
            %Link(j).v=Link(i).v+cross(Link(i).w,Link(i).R*Link(j).b);
            %Link(j).w=Link(i).w+Link(i).R*Link(j).a*Link(j).dq;%计算dq
            %向量化
            tempRb=sum(bsxfun(@times,Link(i).R,repmat(Link(j).b',[1,1,N])),2);
            Link(j).v=Link(i).v+cross(Link(i).w,tempRb);
            tempRa=sum(bsxfun(@times,Link(i).R,repmat(Link(j).a',[1,1,N])),2);
            Link(j).w=Link(i).w+bsxfun(@times,tempRa,Link(j).dq);
        end
        ForwardVelocity(Link(j).brotherID);
        ForwardVelocity(Link(j).childID);
        end
        function P=calP(j)
        if j==0;
            P=0;
        else
            %tempc=Link(j).R*Link(j).c;
            %P=Link(j).m*(Link(j).v+cross(Link(j).w,tempc));
            %P=P+calP(Link(j).brotherID)+calP(Link(j).childID);
            tempRc=sum(bsxfun(@times,Link(j).R,Link(j).c),2);
            Link(j).P=Link(j).m*(Link(j).v+cross(tempRc,Link(j).w));
            P=Link(j).P;
            P=P+calP(Link(j).brotherID)+calP(Link(j).childID);
        end
        end
        function L=calL(j)
        if j==0;
            L=0;
        else
            %tempc=Link(j).R*Link(j).c;
            %p=Link(j).m*(Link(j).v+cross(Link(j).w,tempc));
            %L=cross(Link(j).c,p)+Link(j).R*Link(j).I*Link(j).R'*Link(j).w;
            %L=L+calL(Link(j).brotherID)+calL(Link(j).childID);
            tempRI=array3D_multi(Link(j).R,Link(j).I);
            tempRIRt=array3D_multi(tempRI,permute(Link(j).R,[2 1 3]));
            Link(j).L=cross(permute(Link(j).c,[2 1 3]),Link(j).P)+sum(bsxfun(@times,tempRIRt,Link(j).w),2);
            L=Link(j).L;
            L=L+calL(Link(j).brotherID)+calL(Link(j).childID);       
        end
        end
        function Caldq(TargetID,flink)%这里的flink只提供参考的速度和角速度;已向量化
        vb=Link(1).v;
        wb=Link(1).w;
        vt=flink.v;
        wt=flink.w;
        vd=vt-vb-cross(wb,(flink.p-Link(1).p));
        wd=wt-wb;
        %JACOBIAN
        ray=Target_ray(TargetID);
        J=CalJacobian(ray);
        vw=[vd;wd];
        dq=zeros(5,1,length(vw));
        for izm=1:length(J(1,1,:))%取第三维的长度
          dq(:,:,izm)=J(:,:,izm)\vw(:,:,izm);
        end
        %赋值于各个杆
        for qi=2:length(ray)
            j=ray(qi);
            Link(j).dq=dq(qi-1,:,:);%J 2 3 4 5 6
        end
        end
        function ray=Target_ray(TargetID)
     i=Link(TargetID).motherID;
     if i==1
         ray=[1;TargetID];
     else
         ray=[Target_ray(i);TargetID];% ray = 1 2 3
     end
        end
        function J=CalJacobian(ray)
            pn=Link(ray(end)).p;
            %------------向量化---------------------
            jN=length(ray);
            j2=zeros(6,1,N);j3=zeros(6,1,N);j4=zeros(6,1,N);j5=zeros(6,1,N);
            j6=zeros(6,1,N);%eval出现的变量都要再eval之前存在。
        %   有关优化，1:）J是否应该欲加载一下2）j的最后一列，不用算pn-pn,直接写个0
            for i=2:jN   
                a=sum(bsxfun(@times,Link(ray(i)).R,repmat(Link(ray(i)).a',[1,1,N])),2);%a要为行向量
                eval(['j' int2str(i) '=[cross(a,pn-Link(ray(i)).p);a];']);               %cross可以处理三维的。;也可以，看来还是给数组用的
            end
            J=[j2,j3,j4,j5,j6];%J  6x5xN,N表示时刻,6x5表示了5个关节角度。
        end
        function R=array3D_multi(R1,R2)
        tR1=repmat(R1,[],3);
        tR2=permute(R2,[2 1 3]);
        [m,~,tN]=size(tR2);
        idx=1:m;
        idx=idx(ones(1,3),:);
        ttR2=tR2(idx(:),:,:);
        tem=sum(bsxfun(@times,tR1,ttR2),2);
        R=reshape(tem,[3,3,tN]);
end
    end 
    function M=TotalMass(j)
        if j==0
            M=0;
        else
            m=Link(j).m;
            M=m+TotalMass(Link(j).brotherID)+TotalMass(Link(j).childID);
        end
    end
    function dy=Diff(x,y)
        if length(size(y))==3
            [a,~,b]=size(y);
            y=reshape(y,[a,b]);
            cs=csapi(x,y);
            pp1=fnder(cs);
            dy=fnval(pp1,x);
            dy=reshape(dy,[a,1,b]);
        else 
            cs=csapi(x,y);
            pp1=fnder(cs);
            dy=fnval(pp1,x);
        end
    end
end

%%%%%计算反馈量
function [ del_com,del_p ] = cal_del_x( zmpreal,t,zmpdesign)
del_p=zmpdesign-zmpreal;
%--------------------------A
Zc=310;
g=9800;%单位啊，shit！这里是mm/s^2
dt=t(2)-t(1);
a=-Zc/(g*dt^2);
b=2*Zc/(g*dt^2)+1;
c=-Zc/(g*dt^2);
N=length(t);
diag_1=ones(N,1)*b;
diag_2=ones(N-1,1)*c;
diag_3=ones(N-1,1)*a;
diag_1(1)=a+b;diag_1(end)=b+c;
A=diag(diag_2,1)+diag(diag_1)+diag(diag_3,-1);
%-----------------------A-end
del_comx=A\del_p(1,:)';%结果为列向量 
del_comy=A\del_p(2,:)';
del_comz=zeros(length(t),1);
del_com=[del_comx,del_comy,del_comz]';%第一行x，第二行y
end