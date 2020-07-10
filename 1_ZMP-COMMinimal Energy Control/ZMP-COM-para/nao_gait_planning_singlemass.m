function [yx,yyxx,zmp,totZMP,bodyp1,com,t] = Gait_Planning_LIPM_Singlemass()
% ,Rp,Lp,RR,bodypx,bodypy,comx,comy,comz
clear global;
palse=1;
num=10; T_r=1;  Ts=0.8;  dt=0.02; stepx=40; Zc=320; %%依次为循环周期、步行周期时间、单足相时间、采样间隔、步长、倒立摆质心高度
stepwidth=100;
h_m=20; %%摆动足部最大抬脚高度

global Link N1;
global bodyp;

[zmp,comx,comy,comz,~,Sx]=COM_ZMP_planning(T_r,Ts,dt,stepx,stepwidth,num,Zc);
[Lp,LR,Lv,Lw,Rp,RR,Rv,Rw,t_goal]=Selfplanning_Foot_position(T_r,Ts,dt,stepx,stepwidth,num,h_m,Sx);
bodyp1=[comx;comy;comz];

t=t_goal;
N1=length(t);
t1=t';
zmp_ref = zmp(1:2,:);
if palse==1
    Initial(N1);
    for pii=1:1 %这个for循环用来修正误差ZMP，循环次数为修正次数
        if isempty(bodyp)
            bodyp=reshape(bodyp1,[3,1,N1]);
        end
        [bp,bR,bv,bw]=trans_bp(bodyp,t,N1);        
        Link(1).p=bp;Link(1).R=bR;Link(1).v=bv;Link(1).w=bw;
        R_RefpR.v=Rv;R_RefpR.w=Rw;R_RefpR.p=Rp;R_RefpR.R=RR;
        L_RefpR.v=Lv;L_RefpR.w=Lw;L_RefpR.p=Lp;L_RefpR.R=LR;
        
        Forward_Kinematics(2,N1);
        InverseKinematics(6,R_RefpR);
        InverseKinematics(11,L_RefpR);
        com=calCom();%计算总体质心
        [totZMP,~]=cal_linkZMP(com,R_RefpR,L_RefpR,t);
        [del_com]=cal_del_x(totZMP,t,zmp_ref,Zc);
        del_com2=del_com/2;
        bodyp=bodyp+reshape(del_com2,3,1,[]);%delc-com原为二维矩阵，转换为三维数组
        
        [JR]=InverseKinematics(6,R_RefpR);
        [JL]=InverseKinematics(11,L_RefpR);
        com=calCom();%计算总体质心
        [totZMP,bodyZMP]=cal_linkZMP(com,R_RefpR,L_RefpR,t);
        
        %%%-------------反馈修正ZMP反馈量-------
        
        %         n_ini = round(T_ref/dt);
        %         [zmp_x,zpm_y]=size(bodyZMP);
        %         zmp_ref = zeros(zmp_x,zpm_y);
        %         zmp_ref(1,1:n_t_ini) = bodyZMP(1,1:n_t_ini);
        %         zmp_ref(1,n_t_ini+1:end) = zmp(1,n_ini+1:end);
        %         zmp_ref(2,1:n_t_ini) = bodyZMP(2,1:n_t_ini);
        %         zmp_ref(2,n_t_ini+1:end) = zmp(2,n_ini+1:end);
        [del_com,del_p]=cal_del_x(bodyZMP,t,zmp_ref,Zc);
        bodyp=bodyp+reshape(del_com,3,1,[]);
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

com=com';
zmp=zmp';

clear global bodyp;
clear global Link;
y=reshape(y,[],10);
yx=[t',y];
% yx=data_Interpolation(yx);
% yx=data_Interpolation2(yx);
yyxx = yx(:,2:end)';

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
    function [y]=data_Interpolation(yx)%%初始下蹲插值
        t1=yx(:,1);
        te=t1(end);
        % l=length(t1);
        dt=t1(2)-t1(1);
        N=fix(4/dt);
        yxx=yx(:,2:11);
        % yi=zeros(N,10);
        for i=1:10
            a=yxx(1,i);
            delt=a/N;
            a=a-delt;
            yi(:,i)=[0:delt:a]';
        end
        y1=[yi;yxx];
        tend=te+4;
        t_goal=[dt:dt:tend]';
        y=[t_goal,y1];
    end
    function [y]=data_Interpolation2(yx)%%%结束站立插值
        t1=yx(:,1);
        te=t1(end);
        n=length(t1);
        dt=t1(2)-t1(1);
        N=fix(4/dt);
        yxx=yx(:,2:11);
        % yi=zeros(N,10);
        for i=1:10
            a=yxx(n,i);
            delt=a/N;
            a=a-delt;
            yi(:,i)=[a:-delt:0]';
        end
        y1=[yxx;yi];
        tend=te+4;
        t_goal=[dt:dt:tend]';
        y=[t_goal,y1];
    end
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
I2=[2.76,-0.02,-4.11;-0.02,98003,0.00;-4.11,0.00,8.81];
I3=[1637.48,-0.92,85.88;-0.92,1592.21,-39.18;85.88,-39.18,303.98];
I4=[1182.83,-0.90,28.00,;-0.90,1128.28,-38.48;28.00,-38.48,191.45];
I5=[38.51,-0.06,3.87;-0.06,74.31,-0.00;3.87,-0.00,54.91];
I6=[269.30,5.88,139.13;5.88,643.47,-18.85;139.13,-18.85,525.03];
I7=[2.76,-0.02,-4.08;-0.02,98003,-0.00;-4.08,-0.00,8.81];
I8=[1637.20,-0.92,85.31;-0.92,1591.07,-38.36;85.31,-38.36,303.74];
I9=[1182.08,0.63,36.50,;0.63,1128.65,39.50;36.50,39.50,193.22];
I10=[38.51,-0.03,3.86;-0.03,74.27,-0.02;3.86,-0.02,54.87];
I11=[269.44,-5.70,139.38;-5.70,644.43,18.74;139.38,18.74,525.76];
I_Init={I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11};
v_Init=zeros(11,1);
w_Init=zeros(11,1);
q=[0 0 0 0.15 0 0 0 0 0.15 0 0];%第一个是没有用的
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
    Link(j).q=repmat(q(j),[1,1,tN]);
    if j>1
        Link(j).a=a_Init{1,j-1}';
        Link(j).b=b_Init{1,j-1}';
    end
end
end

function [zmp,comx,comy,comz,tN,Sx]=COM_ZMP_planning(T_r,Ts,dt,stepx,stepwidth,num,Zc)

[zmpx,~,~,tN,Sx]=Self_parabol_planning_ZMP_x(T_r,Ts,dt,stepx,stepwidth,num);
[zmpy,~,~,~]   = Self_parabol_planning_ZMP_y(T_r,Ts,dt,stepx,stepwidth,num);

[comx,comy,comz]=LIPM_cal_ZMP_COM(zmpx,zmpy,tN,dt,Zc);
zmp=[zmpx';zmpy'];

end

function [zmpx,Lfoot_x,Rfoot_x,t,Sx]=Self_parabol_planning_ZMP_x(T_r,Ts,dt,stepx,stepwidth,num)
%%%%--这里只用到了zmpx，单足向ZMP轨迹为抛物线形，X和Y轴上投影单独求解，后面同理-----%%%%%%%%
% num=10; T_r=4;  Ts=3.2;  dt=0.01; stepx=80; stepwidth=256; ml=0.67; mb=8.69;
Td=T_r-Ts;    DSP=Td;
TN=num*T_r;  t=dt:dt:TN;

[~,XN]=size(t);

t_s = Ts;   t_d = Td;   t_f = T_r;  %t_s、t_d、t_f为单足相、双足相、整周期时间
n_t_s = round(t_s/dt);
n_t_d = round(t_d/dt);
n_t_f = round(t_f/dt);

Rfoot_x=zeros(XN,1);
Lfoot_x=zeros(XN,1);
zmpx=zeros(XN,1);         %规划坐标系X轴方向ZMP轨迹,参数调用前先初始化，

an=1; %+1前向，-1反向


Wx=10; Wy=20;  %%Wx、Wy为单足向zmp轨迹在X和Y方向的范围

dl=zeros(10,1);
Sx=stepx*ones(num,1);  Sy=stepwidth*ones(num,1);
% Sx(3,1)=stepx+100;
% Sx(5,1)=stepx-40;
% dl(3)=(40)/120;
% dl(5)=(-40)/120;
Sx(num,1)=0; %%最后一步两足并行

%%=============落脚点控制，第一个落脚点为右脚(foot_x1,foot_y1)(0,-128),第一步迈左脚,落脚点为(foot_x2,foot_y2)-(stepx,128)============%%%%
footx = zeros(1,num);    footy = zeros(1,num);
for j=1:num
    if j == 1
        footx(1,1) = 0;
        footy(1,1) = -stepwidth/2;
        
    else
        footx(1,j) = footx(1,j-1) + Sx(j-1);   footy(1,j) = footy(1,j-1) + Sy(j)*((-1)^(j));
    end
end

%%==========足部规划，认为第一个落脚点为右脚foot_p(:,1)第一步迈左脚,落脚点为foot_p(:,2)=============%%
n_t_n=0;
for nxx=1:1:num
    if nxx==1
        %%%%================ZMP规划X坐标三次多项式插值===============%%%%%
        px_ii=0;   px_id=(an)*Wx;   vpx_ii=0;     vpx_id=0;
        ttt_i = [(0)^3,(0)^2,(0)^1,1;(t_s)^3,(t_s)^2,t_s,1;3*(0)^2,2*(0),1,0;3*(t_s)^2,2*(t_s),1,0];
        ay_i=ttt_i\[px_ii;px_id;vpx_ii;vpx_id];
        pax0=ay_i(4);pax1=ay_i(3);pax2=ay_i(2);pax3=ay_i(1);
        %%注意，ZMP规划Z坐标为0，Y坐标为四次多项式插值，保证单足相轨迹为抛物线
        
        %%%===========================zmp求解========================%%%%
        for i=1:1:n_t_s
            ti=dt*i;
            zmpx(n_t_n+n_t_d+i,1)=pax3*ti^3+pax2*ti^2+pax1*ti+pax0+footx(nxx);
        end
    else
        %%%%================ZMP规划X坐标三次多项式插值===============%%%%%
        px_ii=-(an)*Wx-dl(nxx)*Wx;   px_id=(an)*Wx+dl(nxx)*Wx;   vpx_ii=0;     vpx_id=0;
        ttt_i = [(0)^3,(0)^2,(0)^1,1;(t_s)^3,(t_s)^2,t_s,1;3*(0)^2,2*(0),1,0;3*(t_s)^2,2*(t_s),1,0];
        ay_i=ttt_i\[px_ii;px_id;vpx_ii;vpx_id];
        pax0=ay_i(4);pax1=ay_i(3);pax2=ay_i(2);pax3=ay_i(1);
        for i=1:n_t_s
            ti=dt*i;
            zmpx(n_t_n+n_t_d+i,1)=pax3*ti^3+pax2*ti^2+pax1*ti+pax0+footx(nxx);
        end
    end
    n_t_n=n_t_n+n_t_f;
end

k_n=0;
for i=1:1:num
    
    %%%===================================双足相ZMP轨迹插值==============================%%%%
    %%单足相ZMP轨迹有规划，根据双足相始末位置、速度，进行五次多项式插值%%
    %%%这里五次多项式插值形式：
    %%%y=a5*t^5+a4*t^4+a3*t^3+a2*t^2+a1*t^1+a0*t^0,六个未知系数，通过初始和结束时刻位置、速度、加速度计算系数，为了方便计算，系数的求解在后面写成了矩阵形式，ax为6X6的系数矩阵
    if i==1
        px_ini=0;    px_end=zmpx(k_n+n_t_d+1,1);     vpx_end=(zmpx(k_n+n_t_d+2,1)-zmpx(k_n+n_t_d+1,1))/dt;  vpx_ini=vpx_end;   apx_ini=0;  apx_end=0;
        ax(i,:)=[1,0,0,0,0,0; 1,DSP,(DSP)^(2),(DSP)^(3),(DSP)^(4),(DSP)^(5);...
            0,1,0,0,0,0; 0,1,2*DSP,3*(DSP)^(2),4*(DSP)^(3),5*(DSP)^(4);...
            0,0,2,0,0,0; 0,0,2,6*DSP,12*(DSP)^(2),20*(DSP)^(3)  ]\[px_ini;px_end;vpx_ini;vpx_end;apx_ini;apx_end];
        
        a0x=ax(i,1);  a1x=ax(i,2);a2x=ax(i,3);a3x=ax(i,4);a4x=ax(i,5);a5x=ax(i,6);
        
        for k=1:n_t_d
            zmpx(k_n+k,1)= a0x+a1x*(k*dt)+a2x*(k*dt)^2+a3x*(k*dt)^3+a4x*(k*dt)^4+a5x*(k*dt)^5;   %注意已经是全局坐标系
        end
    else
        px_init=zmpx(k_n,1);   px_end=zmpx(k_n+n_t_d+1,1);   vpx_init=(zmpx(k_n,1)-zmpx(k_n-1,1))/dt; vpx_end=vpx_init; apx_init=0;   apx_end=0;
        ax(i,:)=[1,0,0,0,0,0; 1,DSP,(DSP)^(2),(DSP)^(3),(DSP)^(4),(DSP)^(5);...
            0,1,0,0,0,0; 0,1,2*DSP,3*(DSP)^(2),4*(DSP)^(3),5*(DSP)^(4);...
            0,0,2,0,0,0; 0,0,2,6*DSP,12*(DSP)^(2),20*(DSP)^(3)  ]\[px_init;px_end;vpx_init;vpx_end;apx_init;apx_end];
        
        a0x=ax(i,1);  a1x=ax(i,2);a2x=ax(i,3);a3x=ax(i,4);a4x=ax(i,5);a5x=ax(i,6);
        
        for k=1:n_t_d
            zmpx(k_n+k,1)= a0x+a1x*(k*dt)+a2x*(k*dt)^2+a3x*(k*dt)^3+a4x*(k*dt)^4+a5x*(k*dt)^5;   %注意已经是全局坐标系
        end
    end
    %%%===================================双足相ZMP轨迹插值==============================%%%%
    k_n = k_n + n_t_f;
end

% figure(2);
% plot(zmpx,'--r','LineWidth',2);
% hold on;
% plot(zmpx_eq,'--','LineWidth',2);
% plot(bodypx,'LineWidth',2);
% legend('zmpx','zmpx_e_q','com');
% title('ZMP_X坐标');
end  %%生成zmp、足部位置X方向的位移
function [zmpy,Lfoot_y,Rfoot_y,t,Sy]=Self_parabol_planning_ZMP_y(T_r,Ts,dt,stepx,stepwidth,num)%%生成zmp、足部位置Y方向的位移
% num=10; T_r=4;  Ts=3.2;  dt=0.01; stepx=80; stepwidth=256; ml=0.67; mb=8.69;
Td=T_r-Ts;    DSP=Td;
TN=num*T_r;  t=dt:dt:TN;

[~,XN]=size(t);

t_s = Ts;   t_d = Td;   t_f = T_r;  %t_s、t_d、t_f为单足相、双足相、整周期时间
n_t_s = round(t_s/dt);
n_t_d = round(t_d/dt);
n_t_f = round(t_f/dt);

Wx=20; Wy=5;

zmpy=zeros(XN,1);

Sx=stepx*ones(num,1);   Sy=stepwidth*ones(num,1);
footx = zeros(1,num);    footy = zeros(1,num);
for j=1:num
    if j == 1
        footx(1,1) = 0;
        footy(1,1) = -stepwidth/2;
    else
        footx(1,j) = footx(1,j-1) + Sx(j-1);   footy(1,j) = footy(1,j-1) + Sy(j)*((-1)^(j));
    end
end

%%==========足部规划，认为第一个落脚点为右脚foot_p(:,1)第一步迈左脚,落脚点为foot_p(:,2)=============%%

Lfoot_y=stepwidth/2*ones(XN,1);
Rfoot_y=-stepwidth/2*ones(XN,1);

n_t_n=0;
for nxx=1:1:num
    py_ii=0;   py_id=0;  py_if=Wy*(-1)^(nxx);    vpy_if=0;
    ttt_i=[(0)^3,(0)^2,0,1;
        (t_s/2)^3,(t_s/2)^2,t_s/2,1;
        (t_s)^3,(t_s)^2,t_s,1;
        3*(t_s/2)^2,2*(t_s/2),1,0];
    ay_i=ttt_i\[py_ii;py_if;py_id;vpy_if];
    pay0=ay_i(4);pay1=ay_i(3);pay2=ay_i(2);pay3=ay_i(1);  %%注意，ZMP规划Z坐标为0，Y坐标为三次多项式插值，保证单足相轨迹为抛物线
    
    for i=1:1:n_t_s
        ti=dt*i;
        zmpy(n_t_n+n_t_d+i,1)=pay3*ti^3+pay2*ti^2+pay1*ti+pay0+footy(nxx);
    end
    n_t_n=n_t_n+n_t_f;
end

%%%=========ZMP规划=================%%
k_n=0;
ay=zeros(num,6);
for i=1:1:num
    
    %%%===================================双足相ZMP轨迹插值==============================%%%%
    %%单足相ZMP轨迹有规划，根据双足相始末位置、速度，进行五次多项式插值%%
    if i==1
        dsp = DSP+0.02;
        py_ini=0;   py_end=zmpy(n_t_d+1,1);  vpy_end=(zmpy(n_t_d+2,1)-zmpy(n_t_d+1,1))/dt;  vpy_ini=0;   apy_ini=0;  apy_end=0; %py_end=footy(i)+width*(-1)^(i+1)+delta*(-1)^;
        ay(i,:)=[1,0,0,0,0,0; 1,dsp,(dsp)^(2),(dsp)^(3),(dsp)^(4),(dsp)^(5);...
            0,1,0,0,0,0; 0,1,2*dsp,3*(dsp)^(2),4*(dsp)^(3),5*(dsp)^(4);...
            0,0,2,0,0,0; 0,0,2,6*dsp,12*(dsp)^(2),20*(dsp)^(3)  ]\[py_ini;py_end;vpy_ini;vpy_end;apy_ini;apy_end];
        
        a0y=ay(i,1);  a1y=ay(i,2);a2y=ay(i,3);a3y=ay(i,4);a4y=ay(i,5);a5y=ay(i,6);
        
        for k=1:n_t_d
            zmpy(k_n+k,1)= a0y+a1y*(k*dt)+a2y*(k*dt)^2+a3y*(k*dt)^3+a4y*(k*dt)^4+a5y*(k*dt)^5;   %注意已经是全局坐标系
        end
    else
        %footy(i-1)+width*(-1)^(i)+delta*(-1)^(i+1);footy(i)+width*(-1)^(i+1)+delta*(-1)^(i)
        py_ini=zmpy(k_n,1);  py_end=zmpy(k_n+n_t_d+1);  vpy_init=(zmpy(k_n,1)-zmpy(k_n-1,1))/dt; vpy_end=vpy_init; apy_ini=0; apy_end=0;
        ay(i,:)=[1,0,0,0,0,0; 1,dsp,(dsp)^(2),(dsp)^(3),(dsp)^(4),(dsp)^(5);...
            0,1,0,0,0,0; 0,1,2*dsp,3*(dsp)^(2),4*(dsp)^(3),5*(dsp)^(4);...
            0,0,2,0,0,0; 0,0,2,6*dsp,12*(dsp)^(2),20*(dsp)^(3)  ]\[py_ini;py_end;vpy_ini;vpy_end;apy_ini;apy_end];
        
        a0y=ay(i,1);  a1y=ay(i,2);a2y=ay(i,3);a3y=ay(i,4);a4y=ay(i,5);a5y=ay(i,6);
        for k=1:n_t_d
            zmpy(k_n+k,1)= a0y+a1y*(k*dt)+a2y*(k*dt)^2+a3y*(k*dt)^3+a4y*(k*dt)^4+a5y*(k*dt)^5;   %注意已经是全局坐标系
        end
    end
    k_n = k_n + n_t_f;
end

% figure(2);
% hold on;
% plot(t,zmpy,'r--','LineWidth',2);
% plot(t,zmpy_eq);
% legend('zmpydesign','zmpyequivalent');

%-----各个时刻的ZMP位置（相对于全局坐标系）====


end
function [comx,comy,comz]=LIPM_cal_ZMP_COM(px,py,tN,dt,Zc)  %%三维线性倒立摆模型，ZMP公式的简化解法，参见《仿人机器人》

g=9800;%单位是mm/s^2
t=tN;
%重合段相同的时间点的重心位移不同，比如,0:2和2:4在2s时的重心有差别，这里的方法是在
%两端各加一段时间，再取中间的时间。
N=length(t);
a=-Zc/(g*dt^2);
b=2*Zc/(g*dt^2)+1;
c=-Zc/(g*dt^2);
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
% vz=Diff(t,pz);
% pz(1)=pz(1)+a*vz(1)*dt;
% pz(end)=pz(end)-a*vz(end)*dt;

com_z=Zc*ones(N,1);
comx=com_x';comy=com_y';comz=com_z';
bodyp1=[com_x com_y com_z]';

bodyp=reshape(bodyp1,[3,1,N]);

bodyv=Diff(t,bodyp);%Diff把各个时刻的p集中在一起，最后第三维存的是求导后的时间。


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

function [Lp,LR,Lv,Lw,Rp,RR,Rv,Rw,t_goal]=Selfplanning_Foot_position(T_r,Ts,dt,stepx,stepwidth,num,h_m,Sx)         %num,T_r,Ts,dt,stepx,stepwidth
% num=10; T_r=4;  Ts=3.2;  dt=0.005; stepx=80; stepwidth=128; h_m=30;

tN=num*T_r;
N=fix(tN/dt);

%Nf=N+n_t_i;
% t_goal=dt:dt:Nf*dt;
t_goal=dt:dt:tN;
Nf=N;
t_s = Ts;
t_d = T_r-Ts;
t_f = T_r;  %t_s、t_d、t_f为单足相、双足相、整周期时间
n_t_s = round(t_s/dt);
n_t_d = round(t_d/dt);
n_t_f = round(t_f/dt);


Rfoot_x=zeros(1,N); Rfoot_y=zeros(1,N); Rfoot_h=zeros(1,N);
Lfoot_x=zeros(1,N); Lfoot_y=zeros(1,N); Lfoot_h=zeros(1,N);

Rfoot_q=zeros(1,Nf);  Lfoot_q=zeros(1,Nf);  %Lfoot_p=zeros(3,Nf);Rfoot_p=zeros(3,Nf);

footx = zeros(1,num);    footy = zeros(1,num);    Sy=stepwidth*ones(num,1); %Sx=stepx*ones(num,1);
for j=1:num
    if j == 1
        footx(1,1) = 0;
        footy(1,1) = -stepwidth/2;
    else
        footx(1,j) = footx(1,j-1) + Sx(j-1);   footy(1,j) = footy(1,j-1) + Sy(j)*((-1)^(j));
    end
end


%%%===============足部高度四次多项式规划================%%%
footh_i=0; footh_h=h_m; footh_d=0; footh_iv=0; footh_dv=0; footh_hv=0;
ttt_i = [0^4,(0)^3,(0)^2,0,1;(t_s/2)^4,(t_s/2)^3,(t_s/2)^2,t_s/2,1;(t_s)^4,(t_s)^3,(t_s)^2,t_s,1;4*(0)^3,3*(0)^2,2*(0),1,0;4*(t_s)^3,3*(t_s)^2,(2*t_s),1,0;4*(t_s/2)^3,3*(t_s/2)^2,2*(t_s/2),1,0];
ay_i=ttt_i\[footh_i;footh_h;footh_d;footh_iv;footh_dv;footh_hv];
ah0=ay_i(5);ah1=ay_i(4);ah2=ay_i(3);ah3=ay_i(2); ah4=ay_i(1);
n_t_n=0;
for nxx=1:1:num
    if nxx==1
        ttt_i = [(0)^3,(0)^2,(0)^1,1;(t_s)^3,(t_s)^2,t_s,1;(0)^2,0,1,0;3*(t_s)^2,2*(t_s),1,0];
        footx_i=0;  footx_d=Sx(nxx);  footx_iv=0;  footx_dv=0;
        ay_i=ttt_i\[footx_i;footx_d;footx_iv;footx_dv];
        fx0=ay_i(4);fx1=ay_i(3);fx2=ay_i(2); fx3=ay_i(1);
        %%足部X坐标三次多项式插值，footx=a3*t^3+a2*t^2+a1*t+a0;  %%注意，足部Y坐标不变，高度为四次多项式插值
        for i=1:n_t_s
            ti=i*dt;
            Rfoot_x(1,n_t_n+n_t_d+i) = footx(nxx);
            Rfoot_y(1,n_t_n+n_t_d+i)=-stepwidth/2;
            Rfoot_h(1,n_t_n+n_t_d+i)=0;
            
            Lfoot_x(1,n_t_n+n_t_d+i) = fx3*ti^3+fx2*ti^2+fx1*ti+fx0+footx(nxx);
            Lfoot_y(1,n_t_n+n_t_d+i)=stepwidth/2;
            Lfoot_h(1,n_t_n+n_t_d+i)=ah4*ti^4+ah3*ti^3+ah2*ti^2+ah1*ti+ah0;
        end
        
        for i=1:n_t_d
            Lfoot_x(1,n_t_n+i) = Lfoot_x(1,n_t_n+n_t_d+1);     Lfoot_y(1,n_t_n+i)= stepwidth/2;    Lfoot_h(1,n_t_n+i)=0;
            Rfoot_x(1,n_t_n+i) = Rfoot_x(1,n_t_n+n_t_d+1);     Rfoot_y(1,n_t_n+i)=-stepwidth/2;    Rfoot_h(1,n_t_n+i)=0;
        end
    elseif nxx==3
        footx_i=-Sx(nxx-1);  footx_d=Sx(nxx);  footx_iv=0;  footx_dv=0;
        ay_i=ttt_i\[footx_i;footx_d;footx_iv;footx_dv];
        fx0=ay_i(4);fx1=ay_i(3);fx2=ay_i(2); fx3=ay_i(1);
        %%足部X坐标三次多项式插值，footx=a3*t^3+a2*t^2+a1*t+a0; %%注意，足部Y坐标不变，高度为四次多项式插值
        for i=1:n_t_s
            ti=i*dt;
            Rfoot_x(1,n_t_n+n_t_d+i) = footx(nxx);
            Rfoot_y(1,n_t_n+n_t_d+i)=-stepwidth/2;
            Rfoot_h(1,n_t_n+n_t_d+i)=0;
            Lfoot_x(1,n_t_n+n_t_d+i) = fx3*ti^3+fx2*ti^2+fx1*ti+fx0+footx(nxx);
            Lfoot_y(1,n_t_n+n_t_d+i)=stepwidth/2;
            Lfoot_h(1,n_t_n+n_t_d+i)=ah4*ti^4+ah3*ti^3+ah2*ti^2+ah1*ti+ah0;
        end
        for i=1:n_t_d
            Lfoot_x(1,n_t_n+i) = Lfoot_x(1,n_t_n+n_t_d+1);     Lfoot_y(1,n_t_n+i)= stepwidth/2;    Lfoot_h(1,n_t_n+i)=0;
            Rfoot_x(1,n_t_n+i) = Rfoot_x(1,n_t_n+n_t_d+1);     Rfoot_y(1,n_t_n+i)=-stepwidth/2;    Rfoot_h(1,n_t_n+i)=0;
        end
    else
        %ttt_i = [(0)^3,(0)^2,(0)^1,1;(t_s)^3,(t_s)^2,t_s,1;(0)^2,0,1,0;3*(t_s)^2,2*(t_s),1,0];
        footx_i=-Sx(nxx-1);  footx_d=Sx(nxx);  footx_iv=0;  footx_dv=0;
        ay_i=ttt_i\[footx_i;footx_d;footx_iv;footx_dv];
        fx0=ay_i(4);fx1=ay_i(3);fx2=ay_i(2); fx3=ay_i(1);
        %%足部X坐标三次多项式插值，footx=a3*t^3+a2*t^2+a1*t+a0; %%注意，足部Y坐标不变，高度为四次多项式插值
        for i=1:n_t_s
            ti=i*dt;
            if rem(nxx,2)==0   %%%左脚支撑
                Rfoot_x(1,n_t_n+n_t_d+i) = fx3*ti^3+fx2*ti^2+fx1*ti+fx0+footx(nxx);
                Rfoot_y(1,n_t_n+n_t_d+i)=-stepwidth/2;
                Rfoot_h(1,n_t_n+n_t_d+i)=ah4*ti^4+ah3*ti^3+ah2*ti^2+ah1*ti+ah0;
                Lfoot_x(1,n_t_n+n_t_d+i) = footx(nxx);
                Lfoot_y(1,n_t_n+n_t_d+i)=stepwidth/2;
                Lfoot_h(1,n_t_n+n_t_d+i)=0;
            else
                Rfoot_x(1,n_t_n+n_t_d+i) = footx(nxx);
                Rfoot_y(1,n_t_n+n_t_d+i)=-stepwidth/2;
                Rfoot_h(1,n_t_n+n_t_d+i)=0;
                Lfoot_x(1,n_t_n+n_t_d+i) = fx3*ti^3+fx2*ti^2+fx1*ti+fx0+footx(nxx);
                Lfoot_y(1,n_t_n+n_t_d+i)=stepwidth/2;
                Lfoot_h(1,n_t_n+n_t_d+i)=ah4*ti^4+ah3*ti^3+ah2*ti^2+ah1*ti+ah0;
            end
        end
        %%%双足相足部位置
        for i=1:n_t_d
            Lfoot_x(1,n_t_n+i) = Lfoot_x(1,n_t_n+n_t_d+1);     Lfoot_y(1,n_t_n+i)= stepwidth/2;    Lfoot_h(1,n_t_n+i)=0;
            Rfoot_x(1,n_t_n+i) = Rfoot_x(1,n_t_n+n_t_d+1);     Rfoot_y(1,n_t_n+i)=-stepwidth/2;    Rfoot_h(1,n_t_n+i)=0;
        end
    end
    n_t_n=n_t_n+n_t_f;
end

Rfoot_p=[Rfoot_x;Rfoot_y;Rfoot_h];
Lfoot_p=[Lfoot_x;Lfoot_y;Lfoot_h];



t=t_goal;
% N=Nf;
Rfoot_q=reshape(Rfoot_q,[1,1,N]);
RR=Rodrigues([0 1 0],Rfoot_q);
Rpbm=reshape( Rfoot_p,[3,1,N]);
Rp=Rpbm-sum(bsxfun(@times,RR,repmat([0 0 -96],[1,1,N])),2);
Rv=Diff(t_goal,Rp);
Rw=bsxfun(@times,repmat([0 1 0]',[1,1,N]),Diff(t_goal,Rfoot_q));%foot的w曲线,若有需要可取出fnder，这里fw是一个矩阵，列数为length(t)
Rw=reshape(Rw,[3,1,N]);

Lfoot_q=reshape(Lfoot_q,[1,1,N]);
LR=Rodrigues([0 1 0],Lfoot_q);
Lpbm=reshape( Lfoot_p,[3,1,N]);
Lp=Lpbm-sum(bsxfun(@times,LR,repmat([0 0 -96],[1,1,N])),2);
Lv=Diff(t_goal,Lp);
Lw=bsxfun(@times,repmat([0 1 0]',[1,1,N]),Diff(t,Lfoot_q));%foot的w曲线,若有需要可取出fnder，这里fw是一个矩阵，列数为length(t)
Lw=reshape(Lw,[3,1,N]);

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
end

function [ e_out ] = Rodrigues( a,q )       %把单位矢量旋转轴a变为帽子矩阵（斜对称矩阵）

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
function Forward_Kinematics(j,tN)  %%正运动学  可模块化使用
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
    %     Link(j).mp=Link(i).p+sum(bsxfun(@times,Link(i).R,repmat(Link(j).c,[1,1,tN])),2);
    [~,~,N]=size(Link(i).R);
    a=Link(i).p;
    aR=Link(i).R;
    
    for n=1:N
        b(:,n)=a(:,n)+aR(:,:,n)*Link(j).b;
    end
    %Link(j).R=Link(i).R*Rodrigues(Link(j).a, Link(j).q);
    %注意这里的p和R已经是j相对于世界坐标系的位置和姿态，因为母杆坐标的p和R
    %是对世界坐标的，这相当于0Tj=0Ti * iTj
    %下面的语句（到end结束）是实现两个3x3xN的数组在第三维的相乘
    %     temR=repmat(Link(i).R,[],3);
    temR=repmat(Link(i).R,3,1);
    rot=Rodrigues(Link(j).a,Link(j).q);
    rot2=permute(rot,[2 1 3]);
    [m,~,~]=size(rot2);
    idx=1:m;
    idx=idx(ones(1,3),:);
    rot3=rot2(idx(:),:,:);
    tem=sum(bsxfun(@times,temR,rot3),2);
    Link(j).R=reshape(tem,[3,3,tN]);
    %Link(j).v=Link(i).v+sum(bsxfun(@times,Link(i).R,repmat(Link(j).b',[1,1,tN])),2);
end
Forward_Kinematics(Link(j).brotherID,tN);
Forward_Kinematics(Link(j).childID,tN);
end
function [J]=InverseKinematics(TargetID, posRef)   %%逆运动学，可模块化使用
%TargetID和posRef分别是目标杆的ID和其参考位姿（包括p和R）
global Link;
%给各个关节角度初值,如果都是0则jacobian是奇异状态。
ray=Target_ray(TargetID);
N=20;%迭代次数上限
[~,~,tN]=size(Link(1).p);
delta_q=zeros(5,1,tN);
for n=1:N
    J=CalJacobian(ray);   %%%雅可比矩阵每次循环改变，问题：这样计算出的关节角修正量准确？
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
%%%Jacbion
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
            ccc=repmat(Link(ray(i)).a',[1,1,tN]);
            bsxfun(@times,Link(ray(i)).R,repmat(Link(ray(i)).a',[1,1,tN]));
            a1=sum(bsxfun(@times,Link(ray(i)).R,repmat(Link(ray(i)).a',[1,1,tN])),1);
            a=sum(bsxfun(@times,Link(ray(i)).R,repmat(Link(ray(i)).a',[1,1,tN])),2); %a要为行向量
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
        tR1=repmat(R1,3,1);
        tR2=permute(R2,[2 1 3]);
        [m,~,tN]=size(tR2);
        idx=1:m;
        idx=idx(ones(1,3),:);
        ttR2=tR2(idx(:),:,:);
        tem=sum(bsxfun(@times,tR1,ttR2),2);
        R=reshape(tem,[3,3,tN]);
    end
end
function [bp,bodyR,bodyv,bodyw]=trans_bp(bodyp,t,N)    %%参量维度转换
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
            tR1=repmat(R1,3,1);
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

%%%%计算反馈量
function [ del_com,del_p ] = cal_del_x( zmpreal,t,zmpdesign,Zc)
del_p=zmpdesign-zmpreal;
%--------------------------A

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


