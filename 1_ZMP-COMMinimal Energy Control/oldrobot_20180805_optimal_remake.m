function [yx,y,bodyp2,totZMP,ZMP_body,zmpdesign,Rp,Lp] = calcu20180805_remake
palse=1;
stnum =0;
num =20;
stepwidth=250;
stepx=100;
global Link N1;
%global comx;global bodyp;

%间隔取为两秒，if else使得低电平输入时，不用求解，从而提高效率
dt=0.02;
[t_goal,comx,comy,comz,Rp,RR,Rv,Rw,Lp,LR,Lv,Lw]=plan_com(dt,num,stepwidth,stepx);
bodyp1=[comx;comy;comz];
t=t_goal;
N1=length(t);
if palse==1
    Initial(N1);
    for pii=1:2 %这个for循环用来修正误差ZMP，循环次数为修正次数
        %-----这里调换交换左右脚的位置。--------
    bodyp2=bodyp1';
    bodyp=reshape(bodyp1,[3,1,N1]);
    [bp,bR,bv,bw]=trans_bp(bodyp,t,N1);%trans输出的bodyp依然是三维数组
    Link(1).p=bp;Link(1).R=bR;Link(1).v=bv;Link(1).w=bw;
    R_RefpR.v=Rv;R_RefpR.w=Rw;R_RefpR.p=Rp;R_RefpR.R=RR;
    L_RefpR.v=Lv;L_RefpR.w=Lw;L_RefpR.p=Lp;L_RefpR.R=LR;
    Forward_Kinematics(2,N1);
    InverseKinematics(6,R_RefpR);
    InverseKinematics(11,L_RefpR);
    com=calCom();%计算总体质心
    [totZMP,ZMP_body]=cal_linkZMP(com,R_RefpR,L_RefpR,t);
    [del_com,~,zmpdesign]=cal_del_x(totZMP,t,comx,comy);

    bodyp=bodyp+reshape(del_com,3,1,[]);%delc-com原为二维矩阵，转换为三维数组
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
clear global bodyp;
y=reshape(y,[],10);
yx=[t',y];

figure (11)
hold on;
plot(zmpdesign(1,:),zmpdesign(2,:),'--');
plot(totZMP(1,:),totZMP(2,:),'r');
legend('zmpdesign','totZMP');
figure (12)
plot(t,y);
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
b_Init={[-4 -128 -64],[0 0 -78],[0 0 -230],[0 -2 -230],[0 0 -78],...
    [-4 128 -64],[0 0 -78],[0 0 -230],[0 -2 -230],[0 0 -78]};
c_Init={[0,0,-0.79],[1.62,0.01,-6.13],[0,-0.85,-58.72],...
        [0,-0.4,-114.99],[1.62,0.01,-71.87],[9.9,0,-67.4],...
      [1.62,-0.01,-6.13],[0,0.85,-58.72],[0,0.4,-114.99],...
      [1.62,-0.01,-71.87],[9.9,0,-67.4]};
m_Init=[8.00,1.600,1.783,2.591,1.600,1.121,1.600,1.783,2.591,1.600,1.121];
dq_Init=zeros(1,11);%link(1)的dq没有意义。

I1=[147497.715,0.00,0;0.00,16281.200,0.00;0.00,0.00,163528.242];
I2=[3230.591,-0.121,70.560;-0.121,1547.051,46.679;70.560,46.679,3029.561];
I3=[13916.918,0.00,0.00;0.00,12449.332,971.630;0.00,971.630,2283.117];
I4=[41751.921,0.00,0.00;0.00,41282.239,117.776;0.00,117.776,1316.250];
I5=[11439.657,-0.121,-272.800;-0.121,9756.118,-48.452;-272.800,-48.452,3029.561];
I6=[5749.713,0.00,-811.110;0.00,7164.997,0.00;-811.110,0.00,1876.293];
I7=[3230.591,0.121,70.560;0.121,1547.051,-46.679;70.560,-46.679,3029.561];
I8=[13916.918,0.00,0.00;0.00,12449.332,-971.630;0.00,-971.630,2283.117];
I9=[41751.921,0.00,0.00;0.00,41282.239,-117.776;0.00,-117.776,1316.250];
I10=[11439.657,0.121,-272.800;0.121,9756.118,48.452;-272.800,48.452,3029.561];
I11=[5749.713,0.00,-811.110;0.00,7164.997,0.00;-811.110,0.00,1876.293];
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


function [t_goal,comx,comy,comz,Rp,RR,Rv,Rw,Lp,LR,Lv,Lw]=plan_com(dt,num,stepwidth,stepx)   %%质心轨迹规划
%===========关键参数初始化======
Tsup=0.8;
Tdou=1-Tsup;
T=Tsup+Tdou;
g=9800;
Zc=736;
Tc=sqrt(Zc/g);
%计算capture-point与足部位置调整的位置
t_adjust=Tdou/2;              %一个步行周期内调整步长的时间，分段处理，        ---------------可以考虑去掉
t_des=3.7;                    %计算当前时刻的步行周期。
theta=0;                      %计算capture-point时当前足部偏转角
%---静止状态到稳定状态的过渡
T_ini=4;                      %第一步迈右脚所用时间
t0=dt:dt:T_ini;               %第一步时间采样序列
b1=100;                       %后（前）半足长度
COMy0=130;%133.5              %机器人自然站姿半宽度
shor=20;
%==比较有意思的是这里需要首先确定稳定状态下倒立摆的质心和轨迹=======
[comx,comy,footx1,footy1,t1,SSP,DSP]=position_capture_region_Inout(stepwidth,stepx,t_adjust,t_des,theta,Tc,Tsup,num,dt);
%注意这里现有双足相，再有单足相！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
%加入一段，选取comy第一个极小值对应的下标，主要K中可能出现联号的情况。这里由于只用到K2，暂时不需要管（以后改进）。
%右脚迈出第一步：进入第二步的双足相起始阶段
%也就是1s,根据规划，第1s的坐标（71.65,48.81）
N1=fix(T/dt);
tn=length(t1);
                                         
 Dest1=[dt 3*dt T_ini/8 T_ini/4 T_ini/2 3*T_ini/4 T_ini-T/2 T_ini-Tdou T_ini-dt T_ini];       %初始阶段中只完成右脚向前跨一步，并且质心也移动到右脚上
 x1=[0 -0.1 -20 -5 b1/8 b1/4 b1/2 b1-shor comx(N1-1,:) comx(N1,:)];
 y1=[1 2 1*stepwidth/4 stepwidth/2 COMy0 COMy0  (COMy0+stepwidth/2)/2 3*stepwidth/8 comy(N1-1,:) comy(N1,:)]; 

pp=pchip(Dest1,x1);
comx0=ppval(pp,t0);
comy0=ppval(pp2,t0);
comx1=comx((N1+1):tn,:);
comy1=comy((N1+1):tn,:);
pp2=pchip(Dest1,y1);
comx=[comx0,comx1'];
comy=[comy0,comy1'];

Destx=[0 0.01 0.02 0.03 T_ini/2-2*dt T_ini/2-dt T_ini/2 3*T_ini/4 T_ini];
z1=[780 779.95 779.85 779.7 736.41 736.1 736 736 736];
pp3=pchip(Destx,z1);
comz0=ppval(pp3,t0);
comz1=736*ones(1,tn-N1);
comz=[comz0,comz1];
N_goal=length(comx);
t_goal=dt:dt:(dt*N_goal);
%=========足部轨迹=============
%直接利用以前的函数，修改第一步的参数
[Lp,LR,Lv,Lw,Rp,RR,Rv,Rw]=FootpR(t_goal,footx1,footy1,T_ini,dt,SSP,DSP);

end
%%%生成ZMP，COM和足部轨迹


function [x,y,footx,footy,t,SSP,DSP]=position_capture_region_Inout(stepwidth,stepx,t_adjust,t_des,theta,Tc,Tsup,num,dt)
%%%%%
[x,y,Vx,Vy,footx,footy,px,py,t,SSP,DSP]=position_ZMP_COM_threemassmodel_footplace(stepwidth,stepx,t_adjust,Tc,Tsup,num,dt);

end
%计算ZMP-质心位移和速度--


function  [x,y,Vx,Vy,footx,footy,px,py,t,SSP,DSP]=position_ZMP_COM_threemassmodel_footplace(stepwidth,stepx,t_adjust,Tc,Tsup,num,dt)   
%%================三质点计算质心位置========================
%足部位置就是理想的位置，不需要跟踪、t_adjust表示开始调整的时刻、步态序列初始化

N=num;%------周期变成序列，成为可调周期，便于写成统一形式
SSP=Tsup*ones(N,1);
T=1;
TN=N*T;
t=0.02:dt:TN;  %%步行num个周期的采样点数之和
[~,XN]=size(t);
DSP=T*ones(N,1)-SSP;
T_int=0;      %时刻的初始点
LZMPy0=130;   %静止状态的左、右脚步宽
length=75;    %足部前半足长度
width=55;     %足部半宽度
Sx=stepx*ones(N,1);  
Sy=stepwidth*ones(N,1);  %步行参数序列,
%---第一步迈右脚，步宽序列---
Sy(1,:)=LZMPy0+stepwidth/2;


%%%%%%%%%%%%%%%%变量初始化%%%%%%%%%%%%%%%%%
%每个步行单元质心初始位置和终止位置、速度，初始化，
%注意这里位置和速度都是相对于局部坐标系，速度在局部坐标系和全局坐标系中是相同的。
x_init=zeros(N,1);
y_init=zeros(N,1);
Vx_init=zeros(N,1);
Vy_init=zeros(N,1);

x_end=zeros(N,1);
y_end=zeros(N,1);
Vx_end=zeros(N,1);
Vy_end=zeros(N,1);
%------规划的落脚点位置，也就是各个时期局部坐标系的初始位置（相对于全局坐标系）-----
footx=zeros(N,1);
footy=zeros(N,1);

%-------各个时刻的质心位置和速度（相对于全局坐标系） 先定义数组，后放入数据
x=zeros(round((N)*(T)/dt),1);
y=zeros(round((N)*(T)/dt),1);
Vx=zeros(round((N)*(T)/dt),1);
Vy=zeros(round((N)*(T)/dt),1);
%状态向量初始状态
Rx=zeros(round((N)*(T)/dt),2);
Ry=zeros(round((N)*(T)/dt),2);
%双足相期间五次多项式系数初始化
 ax=zeros(N,6);
 ay=zeros(N,6);

%-----各个时刻的ZMP位置（相对于全局坐标系）====
px=zeros(round((N)*(T)/dt),1);
py=zeros(round((N)*(T)/dt),1);
%%%%%%%%%%%%%%几个重要的矩阵，参考论文。
%%%%%%%%%%%%%matlab中双曲函数都已经有了啊。。。。
%计算ZMP位置和质心位置。
%这里要改变对周期的认识，认为第一周期就是从双足相开始的，认为双足相是机器人开始单周期的一个准备阶段。
 for n_int=1:1:N  %%这里的N指一共有N个周期，
      %第一步先有双足相,比较特殊，目前专门写出来------
        if n_int==1    %第一步，认为双足相是每一步的开始，即初始时刻已经完成双足相，并处于双足相初始时刻
            Kn=round(T_int/dt);  
            %局部坐标系------
            footx(n_int)=0;
            footy(n_int)=LZMPy0;  %认为第一步
%----------------------------------------------
            K=round(DSP(n_int)/dt);       %双足相时间序列
            
            %定义单足相的始、末质心位置和速度
            %！！！！！！！不同原参考文献中的质心初末速度会随着步长和单足相变化而变化的情况，这里直接令某一单双足向的初末质心速度相等，且确保零输入相应，；
            %其余的问题交给双足相来处理。双足相的边界条件与要由本周期目标单足相的速度，位置来定。%（是本方法的与原参考文献的区别） 
            %需要针对不同步长继续测试
            %%比有意思的是，这里应该是先有单足相结束时刻的位置和速度，然后根据零输入条件写出初始时刻的位置。
            x_end(n_int)=(Sx(n_int)-10*DSP(n_int)*length)/2;               %跟本周期的步长有关（？？需要更改末状态位置参数）双足相结束时，质心位置应该靠近下一步支撑足，zmp应在下一支撑足边缘
            y_end(n_int)=(-1)^(n_int)*(Sy(n_int)-10*DSP(n_int)*width)/2;   %跟本周期的步宽有关（？？需要更改末状态位置参数）            
            Vx_end(n_int)= Tc*coth(Tc*SSP(n_int)/2)*x_end(n_int);          %LIPM参见《仿人机器人》
            Vy_end(n_int)= Tc*tanh(Tc*SSP(n_int)/2)*y_end(n_int);
            %单足相初始时刻-----------
            x_init(n_int)=-x_end(n_int);              %反对称
            y_init(n_int)=y_end(n_int);               %正对称
            Vx_init(n_int)= Vx_end(n_int);
            Vy_init(n_int)=-Vy_end(n_int);                    
            %=======单足相,最小能量控制，求解析解==
            
            Kssp=round(SSP(n_int)/dt);       %单足相时间序列
            D_T=calcu_D_T(Tc,SSP(n_int));
            E_AT=calcu_E_AT(Tc,SSP(n_int));
            %参见最小能量控制论文，状态转移矩阵以及格拉姆矩阵  所谓状态转移矩阵，是LIPM末端位置、速度计算矩阵，《仿人机器人》
            for i=1:Kssp
            px(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([x_end(n_int);Vx_end(n_int)]-E_AT*[x_init(n_int);Vx_init(n_int)])+footx(n_int);   %单足相ZMP最小能量理论计算
              %注意转移到全局坐标系上
            py(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([y_end(n_int);Vy_end(n_int)]-E_AT*[y_init(n_int);Vy_init(n_int)])+footy(n_int);   %注意转移到全局坐标系上
            end
                    
            %---五次多项式插值的边界条件,由于要用到单足相起始时刻的加速度，这里在程序上先计算的是单足相的位置。
            %求五次多项式。。
            %关于ZMP的五次多项式。。
            px_init=0;   vpx_init=0; apx_init=0;  px_end=px(Kn+K+1,1);   vpx_end=(px(Kn+K+2,1)-px(Kn+K+1,1))/dt; apx_end=0; 
            py_init=0;   vpy_init=0; apy_init=0;  py_end=py(Kn+K+1,1);   vpy_end=(py(Kn+K+2,1)-py(Kn+K+1,1))/dt; apy_end=0; 
            
            ax(n_int,:)=[1,0,0,0,0,0; 1,DSP(n_int),(DSP(n_int))^(2),(DSP(n_int))^(3),(DSP(n_int))^(4),(DSP(n_int))^(5);...
                          0,1,0,0,0,0; 0,1,2*DSP(n_int),3*(DSP(n_int))^(2),4*(DSP(n_int))^(3),5*(DSP(n_int))^(4);...
                          0,0,2,0,0,0; 0,0,2,6*DSP(n_int),12*(DSP(n_int))^(2),20*(DSP(n_int))^(3)  ]\[px_init;px_end;vpx_init;vpx_end;apx_init;apx_end];
            a0x=ax(n_int,1);  a1x=ax(n_int,2);a2x=ax(n_int,3);a3x=ax(n_int,4);a4x=ax(n_int,5);a5x=ax(n_int,6);
            
            ay(n_int,:)=[1,0,0,0,0,0; 1,DSP(n_int),(DSP(n_int))^(2),(DSP(n_int))^(3),(DSP(n_int))^(4),(DSP(n_int))^(5);...
                          0,1,0,0,0,0; 0,1,2*DSP(n_int),3*(DSP(n_int))^(2),4*(DSP(n_int))^(3),5*(DSP(n_int))^(4);...
                          0,0,2,0,0,0; 0,0,2,6*DSP(n_int),12*(DSP(n_int))^(2),20*(DSP(n_int))^(3)  ]\[py_init;py_end;vpy_init;vpy_end;apy_init;apy_end];
             a0y=ay(n_int,1);  a1y=ay(n_int,2);a2y=ay(n_int,3);a3y=ay(n_int,4);a4y=ay(n_int,5);a5y=ay(n_int,6);         
            for i=1:K
           %多项式，求解
            px(Kn+i,1)= a0x+a1x*(i*dt)+a2x*(i*dt)^2+a3x*(i*dt)^3+a4x*(i*dt)^4+a5x*(i*dt)^5;   %注意已经是全局坐标系
            py(Kn+i,1)= a0y+a1y*(i*dt)+a2y*(i*dt)^2+a3y*(i*dt)^3+a4y*(i*dt)^4+a5y*(i*dt)^5;   %注意已经是全局坐标系   双足相ZMP五次多项式插值
            end      
         %最后要加一句话，把末状态重新赋值为本周期的末状态，为下一周期的双足相的起始状态。
            x_end(n_int)=(Sx(n_int)-10*DSP(n_int)*length)/2;               
            y_end(n_int)=(-1)^(n_int)*(Sy(n_int)-10*DSP(n_int)*width)/2;
            Vx_end(n_int)= Tc*coth(Tc*SSP(n_int)/2)*x_end(n_int);
            Vy_end(n_int)= Tc*tanh(Tc*SSP(n_int)/2)*y_end(n_int);              
        else   %第一个周期后的周期
            %%%%周期性处理双足相%%%%%%%%%%
            T_int=T_int+T;
            Kn=round(T_int/dt);
            %---三次多项式插值的边界条件
            %局部坐标系的位置------
            footx(n_int)=footx(n_int-1)+Sx(n_int-1);
            footy(n_int)=footy(n_int-1)+(-1)^(n_int-1)*Sy(n_int-1);
            K=round(DSP(n_int)/dt);   
            %%%%%%%%%！！！！！！！！！！！！！！！！！！！！！！%%%%%%%%
            %为满足单足相零输入相应条件，应先重点规划单足相，
            %定义单足相的始、末质心位置和速度：先定义单足相末端位置和速度
            x_end(n_int)=(Sx(n_int)-10*DSP(n_int)*length)/2;                %跟本周期的步长有关
            y_end(n_int)=(-1)^(n_int)*(Sy(n_int)-10*DSP(n_int)*width)/2;    %跟本周期的步宽有关            
            Vx_end(n_int)= Tc*coth(Tc*SSP(n_int)/2)*x_end(n_int);
            Vy_end(n_int)= Tc*tanh(Tc*SSP(n_int)/2)*y_end(n_int);  
            %单足相初始时刻-----------
            x_init(n_int)=-x_end(n_int);              %反对称
            y_init(n_int)=y_end(n_int);               %正对称
            Vx_init(n_int)= Vx_end(n_int);
            Vy_init(n_int)=-Vy_end(n_int);           
            %=======单足相,最小能量控制，求解析解==
            Kssp=round(SSP(n_int)/dt);
            D_T=calcu_D_T(Tc,SSP(n_int));
            E_AT=calcu_E_AT(Tc,SSP(n_int));            
            %第三周期受外力,持续时间为t_adjust
            
            if n_int==3
                if t_adjust<=DSP(n_int)
                    for i=1:Kssp
                    px(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([x_end(n_int);Vx_end(n_int)]-E_AT*[x_init(n_int);Vx_init(n_int)])+footx(n_int);   %注意转移到全局坐标系上
                    py(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([y_end(n_int);Vy_end(n_int)]-E_AT*[y_init(n_int);Vy_init(n_int)])+footy(n_int);   %注意转移到全局坐标系上
                    det_time=i*dt;   %%%系统时间
                    S_t=calcu_S_T(Tc,det_time);
                    E_At=calcu_E_AT(Tc,i*dt);
                    Rx(Kn+K+i,:)=E_At*[x_init(n_int);Vx_init(n_int)]+S_t*D_T*([x_end(n_int);Vx_end(n_int)]-E_AT*[x_init(n_int);Vx_init(n_int)])/2;
                    Ry(Kn+K+i,:)=E_At*[y_init(n_int);Vy_init(n_int)]+S_t*D_T*([y_end(n_int);Vy_end(n_int)]-E_AT*[y_init(n_int);Vy_init(n_int)])/2;
                    x(Kn+K+i,1)=Rx(Kn+K+i,1) +footx(n_int); %注意转移到全局坐标系上
                    y(Kn+K+i,1)=Ry(Kn+K+i,1) +footy(n_int);  %注意转移到全局坐标系
                    Vx(Kn+K+i,1)=Rx(Kn+K+i,2);
                    Vy(Kn+K+i,1)=Ry(Kn+K+i,2);
                    end
                else
                    tN=round((t_adjust-DSP(n_int))/dt);
                    T_ad=T-t_adjust;
                    %前一段时间不进行调整
                    for i=1:tN
                    px(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([x_end(n_int);Vx_end(n_int)]-E_AT*[x_init(n_int);Vx_init(n_int)])+footx(n_int);   %注意转移到全局坐标系上
                    py(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([y_end(n_int);Vy_end(n_int)]-E_AT*[y_init(n_int);Vy_init(n_int)])+footy(n_int);   %注意转移到全局坐标系上
                    det_time=i*dt;
                    S_t=calcu_S_T(Tc,det_time);
                    E_At=calcu_E_AT(Tc,i*dt);
                    Rx(Kn+K+i,:)=E_At*[x_init(n_int);Vx_init(n_int)]+S_t*D_T*([x_end(n_int);Vx_end(n_int)]-E_AT*[x_init(n_int);Vx_init(n_int)])/2;
                    Ry(Kn+K+i,:)=E_At*[y_init(n_int);Vy_init(n_int)]+S_t*D_T*([y_end(n_int);Vy_end(n_int)]-E_AT*[y_init(n_int);Vy_init(n_int)])/2;
                    x(Kn+K+i,1)=Rx(Kn+K+i,1) +footx(n_int); %注意转移到全局坐标系上
                    y(Kn+K+i,1)=Ry(Kn+K+i,1) +footy(n_int);  %注意转移到全局坐标系
                    Vx(Kn+K+i,1)=Rx(Kn+K+i,2);
                    Vy(Kn+K+i,1)=Ry(Kn+K+i,2);
                    end     
                    D_T=calcu_D_T(Tc,T_ad);
                    E_AT=calcu_E_AT(Tc,T_ad);
                    %步长和步宽调整量
%                     Sx(n_int)=Sx(n_int)+50;
%                     Sy(n_int)=Sy(n_int)-50;
                    x_end(n_int)=(Sx(n_int)-10*DSP(n_int)*length)/2;                %跟本周期的步长有关
                    y_end(n_int)=(-1)^(n_int)*(Sy(n_int)-10*DSP(n_int)*width)/2;    %跟本周期的步宽有关            
                    Vx_end(n_int)= Tc*coth(Tc*SSP(n_int)/2)*x_end(n_int);
                    Vy_end(n_int)= Tc*tanh(Tc*SSP(n_int)/2)*y_end(n_int);                     
                    %后一段时间开始调整落脚点
                    for i=(tN+1):Kssp
                    det_time=(i-tN)*dt;
                    px(Kn+K+i,1)=[sinh(Tc*det_time),cosh(Tc*det_time)]*D_T*([x_end(n_int);Vx_end(n_int)]-E_AT*[x_init(n_int);Vx(Kn+K+tN,1)])+footx(n_int);   %注意转移到全局坐标系上
                    py(Kn+K+i,1)=[sinh(Tc*det_time),cosh(Tc*det_time)]*D_T*([y_end(n_int);Vy_end(n_int)]-E_AT*[y_init(n_int);Vy(Kn+K+tN,1)])+footy(n_int);   %注意转移到全局坐标系上
                    end                                                                                    
                end
            else   
                for i=1:Kssp
                px(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([x_end(n_int);Vx_end(n_int)]-E_AT*[x_init(n_int);Vx_init(n_int)])+footx(n_int);   %注意转移到全局坐标系上
                py(Kn+K+i,1)=[sinh(Tc*i*dt),cosh(Tc*i*dt)]*D_T*([y_end(n_int);Vy_end(n_int)]-E_AT*[y_init(n_int);Vy_init(n_int)])+footy(n_int);   %注意转移到全局坐标系上
                end  
            end
            %---双足相始、末ZMP位置和速度, 
            %求五次多项式。。
            %关于ZMP的五次多项式。。
            px_init=px(Kn,1);   vpx_init=(px(Kn,1)-px(Kn-1,1))/dt; apx_init=0;  px_end=px(Kn+K+1,1);   vpx_end=(px(Kn+K+2,1)-px(Kn+K+1,1))/dt; apx_end=0; 
            py_init=py(Kn,1);   vpy_init=(py(Kn,1)-py(Kn-1,1))/dt; apy_init=0;  py_end=py(Kn+K+1,1);   vpy_end=(py(Kn+K+2,1)-py(Kn+K+1,1))/dt; apy_end=0; 
            
            ax(n_int,:)=[1,0,0,0,0,0; 1,DSP(n_int),(DSP(n_int))^(2),(DSP(n_int))^(3),(DSP(n_int))^(4),(DSP(n_int))^(5);...
                          0,1,0,0,0,0; 0,1,2*DSP(n_int),3*(DSP(n_int))^(2),4*(DSP(n_int))^(3),5*(DSP(n_int))^(4);...
                          0,0,2,0,0,0; 0,0,2,6*DSP(n_int),12*(DSP(n_int))^(2),20*(DSP(n_int))^(3)  ]\[px_init;px_end;vpx_init;vpx_end;apx_init;apx_end];
            a0x=ax(n_int,1);  a1x=ax(n_int,2);a2x=ax(n_int,3);a3x=ax(n_int,4);a4x=ax(n_int,5);a5x=ax(n_int,6);
            
            ay(n_int,:)=[1,0,0,0,0,0; 1,DSP(n_int),(DSP(n_int))^(2),(DSP(n_int))^(3),(DSP(n_int))^(4),(DSP(n_int))^(5);...
                          0,1,0,0,0,0; 0,1,2*DSP(n_int),3*(DSP(n_int))^(2),4*(DSP(n_int))^(3),5*(DSP(n_int))^(4);...
                          0,0,2,0,0,0; 0,0,2,6*DSP(n_int),12*(DSP(n_int))^(2),20*(DSP(n_int))^(3)  ]\[py_init;py_end;vpy_init;vpy_end;apy_init;apy_end];
             a0y=ay(n_int,1);  a1y=ay(n_int,2);a2y=ay(n_int,3);a3y=ay(n_int,4);a4y=ay(n_int,5);a5y=ay(n_int,6);         
            for i=1:K
           %多项式，求解
            px(Kn+i,1)= a0x+a1x*(i*dt)+a2x*(i*dt)^2+a3x*(i*dt)^3+a4x*(i*dt)^4+a5x*(i*dt)^5;   %注意已经是全局坐标系
            py(Kn+i,1)= a0y+a1y*(i*dt)+a2y*(i*dt)^2+a3y*(i*dt)^3+a4y*(i*dt)^4+a5y*(i*dt)^5;   %注意已经是全局坐标系
            end                  
        end
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%最后的大招，对整个时间历程的ZMP进行求解。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pz=736*ones(XN,1);
[~, ~, ~, v,comx,comy]=ZMP2COM(t,px,py,pz);
for i=1:K
    x(Kn+i,1)=comx(i,:);
    y(Kn+i,1)=comy(i,:);
end    
x=comx;
y=comy;
Vx=v(1,:)';
Vy=v(2,:)';

%------------画图---    
    figure(1)
    hold on;
    plot(footx(1:N,:),footy(1:N,:),'r');
    plot(px(1:XN,:),py(1:XN,:),'--');
    plot(x(1:XN,:),y(1:XN,:),'g');
    legend('理想的落地点','ZMP位移','质心轨迹');
    xlabel('x(mm)');ylabel('y(mm)');
   
    figure(2)
    subplot(2,1,1);
    plot(t,Vx(1:XN,:));
    xlabel('t(s)');ylabel('Vx(mm/s)');
    subplot(2,1,2);
    plot(t,Vy(1:XN,:));
    xlabel('t(s)');ylabel('Vy(mm/s)');   
    
    figure(3)
    subplot(2,2,1);
    plot(t,x(1:XN,:));
    xlabel('t(s)');ylabel('x(mm/s)');
    subplot(2,2,2);
    plot(t,y(1:XN,:));
    xlabel('t(s)');ylabel('y(mm/s)'); 
    subplot(2,2,3);
    plot(t,px(1:XN,:));
    xlabel('t(s)');ylabel('px(mm/s)');
    subplot(2,2,4);
    plot(t,py(1:XN,:));
    xlabel('t(s)');ylabel('py(mm/s)');    

end

function E_AT=calcu_E_AT(Tc,Tsup)
   E_AT=[cosh(Tc*Tsup),   (sinh(Tc*Tsup))/Tc;...
        Tc*sinh(Tc*Tsup), cosh(Tc*Tsup)];
end

function D_T=calcu_D_T(Tc,Tsup)
D_T=2/(sinh((Tc*Tsup)^2)-(Tc*Tsup)^2)*[Tc*Tsup*cosh(Tc*Tsup)+sinh(Tc*Tsup),  -Tsup*sinh(Tc*Tsup);-Tc*Tsup*sinh(Tc*Tsup),(Tc*Tsup*cosh(Tc*Tsup)-sinh(Tc*Tsup))/Tc];
end

function S_T=calcu_S_T(Tc,Tsup)
S_T=[sinh(Tc*Tsup)-Tc*Tsup*cosh(Tc*Tsup),  -Tc*Tsup*sinh(Tc*Tsup);...
     -(Tc)^2*Tsup*sinh(Tc*Tsup),           -(Tc*Tsup*cosh(Tc*Tsup)+sinh(Tc*Tsup))/Tc];
end

%%
 function [comx,comy]=three_mass_ZMP_COM(px,py,num,mB,mL,T,dt)  %%num-步行周期数，这里计算的质心位移不包括迈右脚的的数据
%%%基于三质点，根据规划的ZMP以及足部轨迹计算腰部位置,周期循环，应考虑下蹲的情况
N=length(px);
M=mB+2*mL;
M_c=mB+mL;
w=sprt((mB+mL)*g/(M*Ez));
Nt=T/dt;
x_B=zeros(N);
y_B=zeros(N); %定义躯干位置序列
 px_ini=zeros(num);
 py_ini=zeros(num);
 px_end=zeros(num);
 py_end=zeros(num);
for i=1:1:num
    px_ini(i)=px(Nt*(i-1)+1);
    py_ini=py(Nt*(i-1)+1);
    px_end=px(Nt*i);
    py_end=py(Nt*i);
end

[Ax,Bx,Ay,By]=three_mass_model_parameter_resolution(px_ini,py_ini,px_end,py_end,T,px,py);%%三质点待定系数

 for i=1:1:num   %%分周期进行计算
    for k=1:1:Nt
    tr=dt*i;
    Ft=diff(add,tr);%%此处有function
    xB((i-1)*Nt+k)=Ax(i)*exp(-w*tr)+Bx(i)*exp(w*tr)+M*(Ft-Ex)/M_c;   %%躯干位置计算
    yB((i-1)*Nt+k)=Ay(i)*exp(-w*tr)+By(i)*exp(w*tr)+M*(Ft-Ex)/M_c;
    end
 end
comx=x_B;
comy=y_B;
 end

%%
function [Ax,Bx,Ay,By]=three_mass_model_parameter_resolution(px,py,T,dt,num,M,Xs,Ys)  
%================根据边界条件确定参数值===================%
%%LIPM求单足相始末质心位置，x=
Ax=zeros(1,1,num);
Ay=zeros(1,1,num);
Bx=zeros(1,1,num);
By=zeros(1,1,num);
Nt=T/dt;
% x_B=zeros(N);
% y_B=zeros(N); %定义躯干位置序列
% px_ini=zeros(1,1,num);
% py_ini=zeros(1,1,num);
% px_end=zeros(1,1,num);
% py_end=zeros(1,1,num);
% for i=1:1:num
%     px_ini(i)=px(Nt*(i-1)+1);
%     py_ini(i)=py(Nt*(i-1)+1);
%     px_end(i)=px(Nt*i);
%     py_end(i)=py(Nt*i);
%  end
[dfx_ini,dfx_end,dfy_ini,dfy_end]=resolution_ZMP(px,py,Xs,Ys,mL,mB,t,M);
[Ex,Ey,~]=parameter_cal(mB,mL,Zc,stepx,stepy);
    for j=1:1:num

    xB_end=(Sx(j)-10*DSP*length)/2;  %%??
    xB_ini=-xB_end;
    Ax(j)=(xB_ini*exp(w*T)-xB_end+M*(dfx_end(j)-Ex)/mc-exp(w*T)*M*(dfx_ini(j)-Ex)/mc)/(exp(w*T)-exp(-w*T));
    Ay(j)=(yB_ini*exp(w*T)-yB_end+M*(dfy_end(j)-Ey)/mc-exp(w*T)*M*(dfy_ini(j)-Ey)/mc)/(exp(w*T)-exp(-w*T));
    Bx(j)=(xB_end-xB_ini*exp(-w*T)+M*(dfx_end(j)-Ex)/mc-exp(-w*T)*M*(dfx_ini(j)-Ex)/mc)/(exp(w*T)-exp(-w*T));
    By(j)=(yB_end-yB_ini*exp(-w*T)+M*(dfy_end(j)-Ey)/mc-exp(-w*T)*M*(dfy_ini(j)-Ey)/mc)/(exp(w*T)-exp(-w*T));
    end
end
%%  -----------------------------------------------
function [dfx_ini,dfx_end,dfy_ini,dfy_end]=resolution_ZMP(px,py,Xs,Ys,mL,M,t)  %%ZMP——COM计算各周期初始位置信息，计算微分方程参数Ax、Ay、Bx、By
dXs=Diff(t,Xs);   %%注意足部位置的维数
ddXs=Diff(t,dXs);
dYs=Diff(t,Ys);
ddYs=Diff(t,dYs);
Fx=px-mL*(Xs-Zs*ddXs/g)/(2*M);
Fy=py-mL*(Ys-Zs*ddYs/g)/(2*M);
dFx=Diff(t,Fx);
dFy=Diff(t,Fy);

dfx_ini=zeros(num,1);
dfy_ini=zeros(num,1);
dfx_end=zeros(num,1);
dfy_end=zeros(num,1);

 for i=1:1:num
     dfx_ini(i)=dFx(Nt*(i-1)+1);
     dfy_ini(i)=dFy(Nt*(i-1)+1);
     dfx_end(i)=dFx(Nt*i);
     dfy_end(i)=dFy(Nt*i);
 end

 
end
%%
function [Ex,Ey,Ez]=parameter_cal(mB,mL,Zc,stepx,stepwidth)  %%简化参数计算
[Zs,Zw,Cxt,Cxs,Ysp,Ysw]=position_threemass_leg_com(Zc,stepx,stepwidth);

Ex=(mB*CxB+mL*Cxs+mL*Cxw)/(mB+mL*2);
Ey=(mB*CyB+mL*Cys+mL*Cyw)/(mB+mL*2);
Ez=(mB*Zc+mL*Zs/2+mL*Zw/2)/(mB+mL*2);
        function [Zs,Zw,Cxs,Cxw,Ysp,Ysw]=position_threemass_leg_com(Zc,stepx,stepwidth,mB,mL)  %%计算初始下蹲姿态支撑腿/Zs,摆动腿/Zw质心高度 可以在仿真环境直接测量
        %通过正、逆运动学计算
        end
    
%%下蹲姿态后质心位置、姿态
[t_goal,comx,comy,comz,Rp,RR,Rv,Rw,Lp,LR,Lv,Lw]=plan_com(dt,num,stepwidth,stepx);
bodyp1=[comx;comy;comz];
t=t_goal;
N1=length(t);
if palse==1
    Initial(N1);
    for pii=1:2 %这个for循环用来修正误差ZMP，循环次数为修正次数
        %-----这里调换交换左右脚的位置。--------
    bodyp2=bodyp1';
    bodyp=reshape(bodyp1,[3,1,N1]);
    [bp,bR,bv,bw]=trans_bp(bodyp,t,N1);%trans输出的bodyp依然是三维数组
    Link(1).p=bp;Link(1).R=bR;Link(1).v=bv;Link(1).w=bw;
    R_RefpR.v=Rv;R_RefpR.w=Rw;R_RefpR.p=Rp;R_RefpR.R=RR;
    L_RefpR.v=Lv;L_RefpR.w=Lw;L_RefpR.p=Lp;L_RefpR.R=LR;
    Forward_Kinematics(2,N1);
    InverseKinematics(6,R_RefpR);
    InverseKinematics(11,L_RefpR);
    com=calCom();%计算总体质心
    [totZMP,ZMP_body]=cal_linkZMP(com,R_RefpR,L_RefpR,t);
    [del_com,~,zmpdesign]=cal_del_x(totZMP,t,comx,comy);

    bodyp=bodyp+reshape(del_com,3,1,[]);%delc-com原为二维矩阵，转换为三维数组
    end 
    c=zeros(N1,10);
    for j=2:11
        c(:,j-1)=Link(j).q(:);%Nx10
    end
    y=c;
end
end
%%

function [bodyp1,bodyp,bodyv,v,comx,comy]=ZMP2COM(rt,px,py,pz)
%scatter:离散 
%parameters
Zc=736;%下面还一个!!!!cal_del_x
%tde=6;%初始下蹲姿势的时间
g=9800;%单位是mm/s^2
dt=rt(2)-rt(1);
% if stnum~=0
%     dohead=3;
% else dohead=0;
% end
%t=rt(1)-dohead:dt:rt(end)+0.5;%调试发现，当输入不同的ZMP段时，由于离散算法的原因
t=rt(1):dt:rt(end);
%重合段相同的时间点的重心位移不同，比如,0:2和2:4在2s时的重心有差别，这里的方法是在
%两端各加一段时间，再取中间的时间。
a=-Zc/(g*dt^2);
b=2*Zc/(g*dt^2)+1;
c=-Zc/(g*dt^2);
N=length(t);
diag_1=ones(N,1)*b;
diag_2=ones(N-1,1)*c;
diag_3=ones(N-1,1)*a;
diag_1(1)=a+b;diag_1(end)=b+c;
A=diag(diag_2,1)+diag(diag_1)+diag(diag_3,-1);
%[px py com_z]=Design_ZMP(t);
vx=Diff(t,px);
px(1)=px(1)+a*vx(1)*dt;
px(end)=px(end)-a*vx(end)*dt;
com_x=A\px;
vy=Diff(t,py);
py(1)=py(1)+a*vy(1)*dt;
py(end)=py(end)-a*vy(end)*dt;
com_y=A\py;
vz=Diff(t,pz);
pz(1)=pz(1)+a*vz(1)*dt;
pz(end)=pz(end)-a*vz(end)*dt;
com_z=A\pz;
bodyp1=[com_x com_y com_z]';
%向量化,change to 3_D array
bodyp=reshape(bodyp1,[3,1,N]);
%repmat可以复制矩阵repmat(eye(3),[1,N])
bodyv=Diff(t,bodyp);%Diff把各个时刻的p集中在一起，最后第三维存的是求导后的时间。
picst=find(abs(t-rt(1))<0.0001);
picen=find(abs(t-rt(end))<0.0001);
comx=com_x(picst:picen,:);
comy=com_y(picst:picen,:);
bodyp1=bodyp1(:,picst:picen);
bodyp=bodyp(:,:,picst:picen);
bodyv=bodyv(:,:,picst:picen);
v=reshape(bodyv,[3,N]);
end   %%基于拉格朗日中值定理求躯干位置

%%
function [Lp,LR,Lv,Lw,Rp,RR,Rv,Rw]=FootpR(t,footx1,footy1,T_ini,dt,SSP,DSP)
%-------------parameter
%注意，输入的t的间隔的达到0.005,0.001基本没有误差。
%本函数返回的fw是角速度矢量，fv是速度矢量，p的每一列是一个时间点的p，
%R的每一列是一个时间点的R，需要用reshape(R(:,i),[3,3])来复原
%后半足长度。
b1=100;
%前半足长度。
b2=100;
T=DSP(1,1)+SSP(1,1);
%步长，注意其意义。
th1=2*pi/180; 
%离地时，足与地面的夹角。
th2=2*pi/180;
%-----------------------------------规划出变步长的足部轨迹-----------------------------
%根据 footx1的维数确定左右足部的摆动周数，周期规划的（第一步）奇数为左脚位置（实际运动第一步迈右脚要手动规划），偶数为右脚位置。
tN=length(t);

K=length(footx1);    
    %起步阶段，右脚迈出第一步
    RLorRini=-130;   %%笛卡尔坐标系，向左为x轴正方向，向右为负方向
    RFirststep=(RLorRini+footy1(2,:))/2;
    Rstep=footy1(2,:);
    LLorRini=130;
    LFirststep=(LLorRini+footy1(3,:))/2;
    Lstep=footy1(3,:);       
    Rt_plan=[dt,1.2,T_ini-2*T,T_ini-T,T_ini-2*DSP(1,1),T_ini-DSP(1,1)/2];  
    Rfoot_plan=[[0;RLorRini;0],[0;RLorRini;0],[b2*(1-cos(th2));RLorRini;b2*sin(th2)],[footx1(2,:)/2;RFirststep;15],...
        [footx1(2,:)-b1*(1-cos(th1));Rstep;b1*sin(th1)],[footx1(2,:);Rstep;0]];
    Rangle_plan=[0,0,th2,th2/2,th1,0];
    %在 T_ini~T_ini+T之内左脚迈出第一步
    Lt_plan=[dt,T_ini,T_ini+DSP(2,1),T_ini+DSP(2,1)+0.2,T_ini+DSP(2,1)+SSP(2,1)/2,T_ini+T-0.02,T_ini+T];
    Lfoot_plan=[[0;LLorRini;0],[0;LLorRini;0],[0;LLorRini;0],[b2*(1-cos(th2));LLorRini;b2*sin(th2)],[footx1(2,:);LFirststep;15],...
        [footx1(3,:)-b1*(1-cos(th1));Lstep;b1*sin(th1)],[footx1(3,:);Lstep;0]];
    Langle_plan=[0,0,0,th2,th2/2,-th1,0]; 
   %注意position_ZMP_COM_energymin_footplace.m足部落脚的序列中第三个点是左脚的第一步落脚点，第四个点是右脚的第二步落脚点
 for i=4:K
     if rem(i,2)==0  %右脚
         Rt_plan=[Rt_plan,(T_ini+fix((2*i-5)/2)*T),(T_ini+fix((2*i-5)/2)*T+DSP((i-1),:)),(T_ini+fix((2*i-5)/2)*T+DSP((i-1),:)+0.02),...
             (T_ini+fix((2*i-5)/2)*T+DSP((i-1),:)+SSP((i-1),:)/2),(T_ini+fix((2*i-3)/2)*T-0.02),(T_ini+fix((2*i-3)/2)*T),(T_ini+fix((2*i-3)/2)*T)+DSP((i),:)];
         Rfoot_plan=[Rfoot_plan,[footx1(i-2,:);footy1(i-2,:);0],[footx1(i-2,:);footy1(i-2,:);0],...
             [footx1(i-2,:)+b2*(1-cos(th2));footy1(i-2,:);b2*sin(th2)],[footx1(i-1,:);(footy1(i,:)+footy1(i-2,:))/2;40],[footx1(i,:)-b1*(1-cos(th1));footy1(i,:);b1*sin(th1)],[footx1(i,:);footy1(i,:);0],[footx1(i,:);footy1(i,:);0]];
         Rangle_plan=[Rangle_plan,0,0,th2,th2/2,-th1,0,0];
     else            %左脚
         Lt_plan=[Lt_plan,(T_ini+fix((2*i-5)/2)*T),(T_ini+fix((2*i-5)/2)*T+DSP((i-1),:)),(T_ini+fix((2*i-5)/2)*T+DSP((i-1),:)+0.02),...
             (T_ini+fix((2*i-5)/2)*T+DSP((i-1),:)+SSP((i-1),:)/2),(T_ini+fix((2*i-3)/2)*T-0.02),(T_ini+fix((2*i-3)/2)*T),(T_ini+fix((2*i-3)/2)*T)+DSP(i,:)];
         Lfoot_plan=[Lfoot_plan,[footx1(i-2,:);footy1(i-2,:);0],[footx1(i-2,:);footy1(i-2,:);0],...
             [footx1(i-2,:)+b2*(1-cos(th2));footy1(i-2,:);b2*sin(th2)],[footx1(i-1,:);(footy1(i,:)+footy1(i-2,:))/2;40],[footx1(i,:)-b1*(1-cos(th1));footy1(i,:);b1*sin(th1)],[footx1(i,:);footy1(i,:);0],[footx1(i,:);footy1(i,:);0]];
         Langle_plan=[Langle_plan,0,0,th2,th2/2,-th1,0,0];         
     end
 end
 Rfoot_pp=pchip(Rt_plan,Rfoot_plan);%至此，pp是3x1的矩阵，三行分别为gait p的x，y，z
 R_aa=pchip(Rt_plan,Rangle_plan);
 Lfoot_pp=pchip(Lt_plan,Lfoot_plan);
 L_aa=pchip(Lt_plan,Langle_plan);
 Rfoot=ppval(Rfoot_pp,t);
 Lfoot=ppval(Lfoot_pp,t);

 %右脚的参数 
Rq2w=ppval(R_aa,t);%采样角度
Rq2w=reshape(Rq2w,[1,1,tN]);
RR=Rodrigues([0 1 0],Rq2w); %%旋转矩阵
Rpbm=reshape( Rfoot,[3,1,tN]);
Rp=Rpbm-sum(bsxfun(@times,RR,repmat([0 0 -122],[1,1,tN])),2);
Rv=Diff(t,Rp);
Rw=bsxfun(@times,repmat([0 1 0]',[1,1,tN]),Diff(t,Rq2w));%foot的w曲线,若有需要可取出fnder，这里fw是一个矩阵，列数为length(t)
Rw=reshape(Rw,[3,1,tN]); 

%左脚参数 
Lq2w=ppval(L_aa,t);%采样角度
Lq2w=reshape(Lq2w,[1,1,tN]);
LR=Rodrigues([0 1 0],Lq2w);
Lpbm=reshape( Lfoot,[3,1,tN]);
Lp=Lpbm-sum(bsxfun(@times,LR,repmat([0 0 -122],[1,1,tN])),2);
Lv=Diff(t,Lp);
Lw=bsxfun(@times,repmat([0 1 0]',[1,1,tN]),Diff(t,Lq2w));%foot的w曲线,若有需要可取出fnder，这里fw是一个矩阵，列数为length(t)
Lw=reshape(Lw,[3,1,tN]); 
end

%%
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
%%
function InverseKinematics(TargetID, posRef)
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
    %--------------向量化下面程序--------
%     if norm(err)<1E-5
%         %mark=1;
%         Forward_Kinematics(2);
%     end
%         break;
%     %矩阵的秩，需调整delta_q=0.5*J\err;
%     delta_q=0.5*pinv(J)*err;
%     for nn=2:length(ray)
%         j=ray(nn);
%         Link(j).q=Link(j).q+delta_q(nn-1);
%     end
%     Forward_Kinematics(2);
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
%%
function [ del_com,del_p,zmpdesign] = cal_del_x(zmpreal,t,comx,comy)
Zc=736;
g=9800;%单位啊，这里是mm/s^2
comx_v=Diff(t,comx); %%求一阶微分
comy_v=Diff(t,comy);
comx_acc=Diff(t,comx_v);
comy_acc=Diff(t,comy_v);
zmpdesignx=comx-Zc/g*comx_acc;     %%改为多刚体的zmp计算方法
zmpdesigny=comy-Zc/g*comy_acc;     %%改为多刚体的zmp计算方法
zmpdesign=[zmpdesignx;zmpdesigny];%%第一行x，第二行y
del_p=zmpdesign-zmpreal;
%--------------------------A
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
%%
function [bp,bodyR,bodyv,bodyw]=trans_bp(bodyp,t,N)  %%求躯干旋转矢量（不变）、速度、角速度（无）
bp=bodyp;
bodyR=reshape(repmat(eye(3),[1,N]),[3,3,N]);
bodyv=Diff(t,bodyp);
bodyw=kron([0 0 0]',ones(1,length(t)));
bodyw=reshape(bodyw,[3,1,N]);
end
%%
function com=calCom(~)%%计算总体质心
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
%%
function [ ZMPall ZMPl1 ] = cal_linkZMP(c,Rflink,Lflink,t)
%该程序返回总体ZMP和躯干ZMP，格式均为2*N，第一行为X坐标。
global Link N;
pz=0;
g=9800;%这里单位是mm,kg,s,则g应取为mm/s^2或mN/kg
[P L]=calPL(Rflink,Lflink);
M=TotalMass(1);
% % % % % %------------------------可考虑用Polyfit处理总体的P和L得到趋势项
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
        function J=CalJacobian(ray)         %雅可比矩阵求解
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
    function M=TotalMass(j)         %总体质量计算
        if j==0
            M=0;
        else
            m=Link(j).m;
            M=m+TotalMass(Link(j).brotherID)+TotalMass(Link(j).childID);
        end
    end
end
%%
function [ e_out ]=Rodrigues(a,q)       %把单位矢量旋转轴a变为帽子矩阵（斜对称矩阵）,再求旋转矩阵
                            
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

%%
function dy=Diff(x,y)
        if length(size(y))==3
            [a,~,b]=size(y);
            y=reshape(y,[a,b]);
            cs=csapi(x,y);
            pp1=fnder(cs);
            dy=fnval(pp1,x);
            dy=reshape(dy,[a,1,b]);
        else 
            cs=csapi(x,y);          %三次样条函数
            pp1=fnder(cs);          %求样条函数的微分
            dy=fnval(pp1,x);        %计算在给定点处的样条函数值
        end
end