function [yx,zmp,totZMP,bodyp1,vcomx,com,t] = Singlemass()
% ,Rp,Lp,RR,bodypx,bodypy,comx,comy,comz
palse=1;
num=6; T_r=4;  Ts=4;  dt=0.02; stepx=100; Zc=710; %%依次为循环周期、步行周期时间、单足相时间、采样间隔、步长、倒立摆质心高度

stepwidth=220;stepwidth_r=220;%%stepwidth为髋部宽度，stepwidth_r为zmp轨迹宽度，以双足边缘间距为参考


h_m=20; %%摆动足部最大抬脚高度

global Link N1;
global bodyp;

[zmp,comx,comy,comz,~,Sx]=COM_ZMP_planning(T_r,Ts,dt,stepx,stepwidth_r,num,Zc);
[Lp,LR,Lv,Lw,Rp,RR,Rv,Rw,t_goal]=Selfplanning_Foot_position(T_r,Ts,dt,stepx,stepwidth,num,h_m,Sx);

bodyp1=[comx;comy;comz];

% vcomx=Diff(t_goal,comx);
% vcomy=Diff(t_goal,comy);
% vcomz=Diff(t_goal,comz);
% 
% acomx=Diff(t_goal,vcomx);
% acomy=Diff(t_goal,vcomy);
% acomz=Diff(t_goal,vcomz);
% 
% jcomx=Diff(t_goal,acomx);
% 
% px=[comx;vcomx;acomx];
% py=[comy;vcomy;acomy];
% pz=[comz;vcomz;acomz];
% zmpx_r=zmp(1,:);
% zmpy_r=zmp(2,:);
% 
%          
% % X_pr=comx';     X_vr=Diff(t_goal',X_pr);    X_ar=Diff(t_goal',X_vr);
% % Y_pr=comy';     Y_vr=Diff(t_goal',Y_pr);    Y_ar=Diff(t_goal',Y_vr);
% % Z_pr=comz';     Z_vr=Diff(t_goal',Z_pr);    
% 
% N=round(Ts/dt);
% % x_k=(x_p;x_v;x_a);
% % y_k=(y_p;y_v;y_a);
% % z_k=(z_p;z_v;z_a);
% 
% A=[1,dt,dt^2;0,1,dt;0,0,1];
% B=[dt^3/6;dt^2/2;dt];
% 
% Cp=[1,0,0]; Cv=[0,1,0]; Ca=[0,0,1];
% P_ps=zeros(N,3); P_pu=zeros(N,N); P_vs=zeros(N,3);
% P_vu=zeros(N,N); P_as=zeros(N,3); P_au=zeros(N,N);
% % X_p=zeros(N,1);  Y_p=zeros(N,1);  Z_p=zeros(N,1);
% % X_j=zeros(N,1);
% for i=1:1:N
%     P_ps(i,:)=Cp*A^(i);
%     P_vs(i,:)=Cv*A^(i);
%     P_as(i,:)=Ca*A^(i);
%     for j=1:1:i
%         if j<=i
%             P_pu(i,j)=Cp*A^(i-j)*B;
%             P_vu(i,j)=Cv*A^(i-1)*B;
%             P_au(i,j)=Ca*A^(i-1)*B;
%         else
%             P_pu(i,j)=0;
%             P_vu(i,j)=0;
%             P_au(i,j)=0;            
%         end
%     end
% end
% 
% C_j=[0;0;0];
% ap=100; av=100; az=10^8;  aa=100;
% Aa=Ca*A;    Av=Cv*A;    Ap=Cp*A;    Ba=Ca*B;    Bv=Cv*B;    Bp=Cp*B;
% td=T_r-Ts;
% n_t_d=round(td/dt);
% kt=3*T_r/dt+n_t_d;
% 
% xk=px(:,kt);  yk=py(:,kt);  zk=pz(:,kt);
% g=9800;
% ZMPx=zmpx_r(1,kt+round(Ts/dt)/2);
% ZMpy=zmpy_r(1,kt+round(Ts/dt)/2);
% pxma=40+ZMPx;
% pxmi=-40+ZMPx;
% pyma=20+ZMpy;
% pymi=-20+ZMpy;
% for i=1:1:N    
%     %xi=[xj;yj;zj]
% %     X_p(i)=P_ps*x_k+P_pu*X_j;
% %     X_v(i)=P_vs*x_k+P_vu*X_j;
% %     X_a(i)=P_as*x_k+P_au*X_j;
% %     px(i),py(i)
% %     x_p(k+i)=A*x_p(k+i-1)+B*X_j;
% %     y_p(k+i)=A*z_p(k+i-1)+B*Y_j;
% %     z_p(k+i)=A*z_p(k+i-1)+B*Z_j;
%     % min: (P-Pr)^2+ P_v^2+P_j^2+(Z-Z_(k-1))^2   
%     %min av*C_v(i)^2+C_j(i)^2+ap(Dpx(i)^2+Dpy(i)^2)+az(Zp(i)-Zp(i-1))^2   
%     x_r=px(1,kt+i);
%     y_r=py(1,kt+i);
%     Q1=ap*(Bp)'*Bp+av*(Bv)'*Bv+aa;
%     Q2=av*(Bv)'*Bv+az*(Bp)'*Bp+aa;
%     [m,n]=size(Q1);
% %     I1=eye(n);
% %     Q1=Q1+I1;
%     Q=[Q1,0,0;0,Q1,0;0,0,Q2];
%     
%         
%     p1=ap*xk'*(Ap)'*Bp-ap*2*x_r'*Bp+av*xk'*(Av)'*Bv;
%     p2=ap*yk'*(Ap)'*Bp-ap*2*y_r'*Bp+av*yk'*(Av)'*Bv;
%     p3=az*zk'*(Ap)'*Bp+av*zk'*(Av)'*Bv-az*2*(Cp*zk)'*Bp;  %%zk是绝对位置，一维数据
%     
%     q=[p1;p2;p3];
%     
% %     f=X'*Q*X+q'*X;
%         
%     a1=Aa*zk*Bp+g*Ba-Ap*(zk)*Ba;
%     b1=Ap*xk*Ba-Aa*xk*Bp-pxma*Bp;
%     c1=pxma*Ap*zk+Ap*zk*Aa*xk-Aa*zk*Ap*xk;
%     
%     a2=Aa*zk*Bp+g*Ba-Ap*zk*Ba;
%     b2=Ap*xk*Ba-Aa*xk*Bp-pxmi*Bp;
%     c2=pxmi*Ap*zk+Ap*zk*Aa*xk-Aa*zk*Ap*xk;
%     
%     a3=Aa*zk*Bp+g*Ba-Ap*zk*Ba;
%     b3=Ap*yk*Ba-Aa*yk*Bp-pyma*Bp;
%     c3=pyma*Ap*zk+Ap*zk*Aa*yk-Aa*zk*Ap*yk;
%     
%     a4=Aa*zk*Bp+g*Ba-Ap*zk*Ba;
%     b4=Ap*yk*Ba-Aa*yk*Bp-pymi*Bp;
%     c4=pymi*Ap*zk+Ap*zk*Aa*yk-Aa*zk*Ap*yk;
%     
%     c5=750-Ap*zk;
%     c6=9800-Aa*zk;
%     
%     P=[a1 0 b1;-a2 0 -b2;0 a3 b3;0 -a4 -b4;0 0 Bp;0 0 Ba];
%     p=[c1;-c2;c3;-c4;c5;c6];
% %     P*[xj,yj,zj]'<=C;
%     cj=quadprog(2*Q,q,P,p);
%     
%     C_j(:,i)=cj;
%     X_j=C_j(1,i);
%     Y_j=C_j(2,i);
%     Z_j=C_j(3,i);
%     x_p(:,i)=A*xk+B*X_j;
%     y_p(:,i)=A*yk+B*Y_j;
%     z_p(:,i)=A*zk+B*Z_j;
%     xk=x_p(:,i);
%     yk=y_p(:,i);
%     zk=z_p(:,i);
% end

% X_j=Cj(1,:)';
% Y_j=Cj(1,:)';
% Z_j=Cj(1,:)';
% comx1=P_ps*x_k+P_pu*X_j;
% comy1=P_ps*y_k+P_pu*y_j;
% comz1=P_ps*y_k+P_pu*y_j;





body=bodyp1';
% zmp_ref=zmp_eq;
t=t_goal;
N1=length(t);
t1=t';
vcomx=Diff(t1,comx)/1000;
v0=zeros(200,1);
vcomx=[v0;vcomx];

if palse==1
    Initial(N1);
    for pii=1:2 %这个for循环用来修正误差ZMP，循环次数为修正次数
        if isempty(bodyp)
            %             bodyp2=bodyp1';
            bodyp=reshape(bodyp1,[3,1,N1]);
        end
        [bp,bR,bv,bw]=trans_bp(bodyp,t,N1);%trans输出的bodyp依然是三维数组
        Link(1).p=bp;Link(1).R=bR;Link(1).v=bv;Link(1).w=bw;
        R_RefpR.v=Rv;R_RefpR.w=Rw;R_RefpR.p=Rp;R_RefpR.R=RR;
        L_RefpR.v=Lv;L_RefpR.w=Lw;L_RefpR.p=Lp;L_RefpR.R=LR;
        Forward_Kinematics(2,N1);
        %         InverseKinematics(6,R_RefpR);
        %         InverseKinematics(11,L_RefpR);
        %         com=calCom();%计算总体质心
        %         [totZMP,~]=cal_linkZMP(com,R_RefpR,L_RefpR,t);
        %         [del_com]=cal_del_x(totZMP,zmp_ref,t,Zc,mb,ml,dt);
        %         del_com2=del_com/2;
        %         bodyp=bodyp+reshape(del_com2,3,1,[]);%delc-com原为二维矩阵，转换为三维数组
        [JR]=InverseKinematics(6,R_RefpR);
        [JL]=InverseKinematics(11,L_RefpR);
        com=calCom();%计算总体质心
        [totZMP,bodyZMP]=cal_linkZMP(com,R_RefpR,L_RefpR,t);
        
        %%%-------------反馈修正ZMP反馈量-------
        zmp_ref = zmp(1:2,:);
        %         n_ini = round(T_ref/dt);
        %         [zmp_x,zpm_y]=size(bodyZMP);
        %         zmp_ref = zeros(zmp_x,zpm_y);
        %         zmp_ref(1,1:n_t_ini) = bodyZMP(1,1:n_t_ini);
        %         zmp_ref(1,n_t_ini+1:end) = zmp(1,n_ini+1:end);
        %         zmp_ref(2,1:n_t_ini) = bodyZMP(2,1:n_t_ini);
        %         zmp_ref(2,n_t_ini+1:end) = zmp(2,n_ini+1:end);
        [del_com,del_p]=cal_del_x(bodyZMP,t,zmp_ref);
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
% COM=reshape(bodyp,3,[]);
% com=COM';
com=com';
zmp=zmp';

clear global bodyp;
clear global Link;
y=reshape(y,[],10);
yx=[t',y];
yx=data_Interpolation(yx);
yx=data_Interpolation2(yx);
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
% b_Init={[-4 -128 -64],[0 0 -78],[0 0 -230],[0 -2 -230],[0 0 -78],...
%     [-4 128 -64],[0 0 -78],[0 0 -230],[0 -2 -230],[0 0 -78]};
b_Init={[-4 -110 -64],[0 0 -78],[0 0 -230],[0 -2 -230],[0 0 -78],...
    [-4 110 -64],[0 0 -78],[0 0 -230],[0 -2 -230],[0 0 -78]};

c_Init={[-0.84,0,-7.28],[1.62,0.01,-6.13],[0,-0.85,-58.72],...
    [0,-0.4,-114.99],[1.62,0.01,-71.87],[9.9,0,-67.4],...
    [1.62,-0.01,-6.13],[0,0.85,-58.72],[0,0.4,-114.99],...
    [1.62,-0.01,-71.87],[9.9,0,-67.4]};
m_Init=[4.68,1.600,1.783,2.591,1.600,1.414,1.600,1.783,2.591,1.600,1.414];
dq_Init=zeros(1,11);%link(1)的dq没有意义。

I1=[254889.296,0.00,601.553;0.00,40146.121,0.00;6011.553,0.00,251158.462];
% I1=[54889.296,0.00,141.553;0.00,8146.121,0.00;141.553,0.00,61158.462];
I2=[3230.591,-0.121,70.560;-0.121,1547.051,46.679;70.560,46.679,3029.561];
I3=[13916.918,0.00,0.00;0.00,12449.332,971.630;0.00,971.630,2283.117];
I4=[41751.921,0.00,0.00;0.00,41282.239,117.776;0.00,117.776,1316.250];
I5=[11439.657,-0.121,-272.800;-0.121,9756.118,-48.452;-272.800,-48.452,3029.561];
I6=[10206.247,0.00,-1920.541;0.00,13657.739,0.00;-1920.541,0.00,6971.306];
I7=[3230.591,0.121,70.560;0.121,1547.051,-46.679;70.560,-46.679,3029.561];
I8=[13916.918,0.00,0.00;0.00,12449.332,-971.630;0.00,-971.630,2283.117];
I9=[41751.921,0.00,0.00;0.00,41282.239,-117.776;0.00,-117.776,1316.250];
I10=[11439.657,0.121,-272.800;0.121,9756.118,48.452;-272.800,48.452,3029.561];
I11=[10206.247,0.00,-1920.541;0.00,13657.739,0.00;-1920.541,0.00,6971.306];
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
    Link(j).q=repmat(q(j),[1,1,tN]);
    if j>1
        Link(j).a=a_Init{1,j-1}';
        Link(j).b=b_Init{1,j-1}';
    end
end
end

function [zmp,comx,comy,comz,tN,Sx]=COM_ZMP_planning(T_r,Ts,dt,stepx,stepwidth,num,Zc)

[zmpx,~,~,tN,Sx]=Self_parabol_planning_ZMP_eq_x(T_r,Ts,dt,stepx,stepwidth,num);
[zmpy,~,~,~]   = Self_parabol_planning_ZMP_eq_y(T_r,Ts,dt,stepx,stepwidth,num);

[comx,comy,comz]=LIPM_cal_ZMP_COM(zmpx,zmpy,tN,dt,Zc);
zmp=[zmpx';zmpy'];

end

function [zmpx,Lfoot_x,Rfoot_x,t,Sx]=Self_parabol_planning_ZMP_eq_x(T_r,Ts,dt,stepx,stepwidth,num)
%%%%--这里只用到了zmpx，单足向ZMP轨迹为抛物线形，X和Y轴上投影单独求解，后面同理-----%%%%%%%%
% num=10; T_r=4;  Ts=3.2;  dt=0.01; stepx=80; stepwidth=256; ml=0.67; mb=8.69;
Td=T_r-Ts;    DSP=Td;
TN=num*T_r;  t=dt:dt:TN;

[~,XN]=size(t);
% MM=mb+2*ml;Czs=Zc/2; Ex=0;
t_s = Ts;   t_d = Td;   t_f = T_r;  %t_s、t_d、t_f为单足相、双足相、整周期时间
n_t_s = round(t_s/dt);
n_t_d = round(t_d/dt);
n_t_f = round(t_f/dt);

Rfoot_x=zeros(XN,1);
Lfoot_x=zeros(XN,1);
zmpx=zeros(XN,1);         %规划坐标系X轴方向ZMP轨迹,参数调用前先初始化，

an=1; %+1前向，-1反向


Wx=40; Wy=20;  %%Wx、Wy为单足向zmp轨迹在X和Y方向的范围

length=128;     width=10;     %足部前半足长度 足部半宽度
dl=zeros(10,1);
Sx=stepx*ones(num,1);  Sy=stepwidth*ones(num,1);
% Sx(3,1)=stepx+100;
% Sx(5,1)=stepx-40;
% dl(3)=(40)/120;
% dl(5)=(-40)/120;
Sx(6,1)=0; %%最后一步两足并行

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
function [zmpy,Lfoot_y,Rfoot_y,t,Sy]=Self_parabol_planning_ZMP_eq_y(T_r,Ts,dt,stepx,stepwidth,num)%%生成zmp、足部位置Y方向的位移
% num=10; T_r=4;  Ts=3.2;  dt=0.01; stepx=80; stepwidth=256; ml=0.67; mb=8.69;
Td=T_r-Ts;    DSP=Td;
TN=num*T_r;  t=dt:dt:TN;

[~,XN]=size(t);

t_s = Ts;   t_d = Td;   t_f = T_r;  %t_s、t_d、t_f为单足相、双足相、整周期时间
n_t_s = round(t_s/dt);
n_t_d = round(t_d/dt);
n_t_f = round(t_f/dt);

Wx=20; Wy=40;

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
    ttt_i=[(0)^3,(0)^2,0,1; (t_s/2)^3,(t_s/2)^2,t_s/2,1; (t_s)^3,(t_s)^2,t_s,1; 3*(t_s/2)^2,2*(t_s/2),1,0];
    ay_i=ttt_i\[py_ii;py_if;py_id;vpy_if];
    pay0=ay_i(4);pay1=ay_i(3);pay2=ay_i(2);pay3=ay_i(1);  %%注意，ZMP规划Z坐标为0，Y坐标为四次多项式插值，保证单足相轨迹为抛物线
    
    for i=1:1:n_t_s
        ti=dt*i;
        zmpy(n_t_n+n_t_d+i,1)=pay3*ti^3+pay2*ti^2+pay1*ti+pay0+footy(nxx);
    end
    n_t_n=n_t_n+n_t_f;
end

%%%=========ZMP规划=================%%
k_n=0;
ay=zeros(num,6);
Ay=zeros(num,6);
for i=1:1:num
    
    %%%===================================双足相ZMP轨迹插值==============================%%%%
    %%单足相ZMP轨迹有规划，根据双足相始末位置、速度，进行五次多项式插值%%
    if i==1
        py_ini=0;   py_end=zmpy(n_t_d+1,1);  vpy_end=(zmpy(n_t_d+2,1)-zmpy(n_t_d+1,1))/dt;  vpy_ini=vpy_end;   apy_ini=0;  apy_end=0; %py_end=footy(i)+width*(-1)^(i+1)+delta*(-1)^;
        ay(i,:)=[1,0,0,0,0,0; 1,DSP,(DSP)^(2),(DSP)^(3),(DSP)^(4),(DSP)^(5);...
            0,1,0,0,0,0; 0,1,2*DSP,3*(DSP)^(2),4*(DSP)^(3),5*(DSP)^(4);...
            0,0,2,0,0,0; 0,0,2,6*DSP,12*(DSP)^(2),20*(DSP)^(3)  ]\[py_ini;py_end;vpy_ini;vpy_end;apy_ini;apy_end];
        
        a0y=ay(i,1);  a1y=ay(i,2);a2y=ay(i,3);a3y=ay(i,4);a4y=ay(i,5);a5y=ay(i,6);
        
        for k=1:n_t_d
            zmpy(k_n+k,1)= a0y+a1y*(k*dt)+a2y*(k*dt)^2+a3y*(k*dt)^3+a4y*(k*dt)^4+a5y*(k*dt)^5;   %注意已经是全局坐标系
        end
    else
        %footy(i-1)+width*(-1)^(i)+delta*(-1)^(i+1);footy(i)+width*(-1)^(i+1)+delta*(-1)^(i)
        py_ini=zmpy(k_n,1);  py_end=zmpy(k_n+n_t_d+1);  vpy_init=(zmpy(k_n,1)-zmpy(k_n-1,1))/dt; vpy_end=vpy_init; apy_ini=0; apy_end=0;
        ay(i,:)=[1,0,0,0,0,0; 1,DSP,(DSP)^(2),(DSP)^(3),(DSP)^(4),(DSP)^(5);...
            0,1,0,0,0,0; 0,1,2*DSP,3*(DSP)^(2),4*(DSP)^(3),5*(DSP)^(4);...
            0,0,2,0,0,0; 0,0,2,6*DSP,12*(DSP)^(2),20*(DSP)^(3)  ]\[py_ini;py_end;vpy_ini;vpy_end;apy_ini;apy_end];
        
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


% Rfoot_p=[Rfoot_pi,Rfoot_p];
% Lfoot_p=[Lfoot_pi,Lfoot_p];
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

%% adaptive foot location + body inclination with time-varying height
function [xk,yk,zk,footx_ref,footy_ref]=NMPC_20190424_use_ownmodel
%%% the term of CoM position error has been replaced with ZMPy Position
%% reference foot location parameter
Tn = 10;
hcom = 0.76;
g = 9.8;
dt = 0.05;
Nh = floor(1.5/dt);

%%% modified the first steps
steplength = 0.2*ones(Tn,1);
steplength(1) = 0;
steplength(3) = 0.2;
% %% change velocity: speed up ---> speed down
% steplength(3) = 0.7;
% steplength(6) = 0.7;
% steplength(7) = 0.7;


stepwidth = 0.256*ones(Tn,1);
stepwidth(1) = stepwidth(1)/2;

%%% stage switch
stepheight =  0*zeros(Tn,1);  %% no stage
%%% stage
% stepheight =  0.1*ones(Tn,1);
% stepheight(4,:) =  0;
% stepheight(3,:) =  0;
% stepheight(6:9,:) =  -0.1*ones(4,1);
% stepheight(9:end,:) = zeros(Tn-8,1);

footx_ref = zeros(Tn,1);
footy_ref = zeros(Tn,1);
footz_ref = zeros(Tn,1);
% footy_ref(1) = stepwidth(1)/2;
for i=2:Tn       %%% singular period ===> right support
    footx_ref(i) = footx_ref(i-1)+steplength(i-1);
    footy_ref(i) = footy_ref(i-1)+(-1)^(i)*stepwidth(i-1);
    footz_ref(i) = footz_ref(i-1)+stepheight(i-1);
end



Ts = 1*ones(Tn,1); Td = 0.0*ones(Tn,1);
n_t_s = round(Ts/dt);
n_t_d = round(Td/dt);
Tx = zeros(Tn,1);
for i=2:Tn
    Tx(i) = Tx(i-1)+Ts(i-1)+Td(i-1);
end
%%%%%time interval also determine the robustness
t = dt:dt:Tx(end);
Nsum = length(t);
nT = round((Ts(1)+Td(1))/dt);
nTs = floor(Ts(1)/dt);
nTd = nT-nTs;

for i=1:Tn
    footx_r(nT*(i-1)+1:nT*i,:) = footx_ref(i);
    footy_r(nT*(i-1)+1:nT*i,:) = footy_ref(i);
    if i==1
    COMXcenter_r(nT*(i-1)+1:nT*(i-1)+nTd,:) = (footx_ref(i))/2;
    COMXcenter_r(nT*(i-1)+nTd+1:nT*(i-1)+nT,:) =  footx_ref(i);
    COMYcenter_r(nT*(i-1)+1:nT*(i-1)+nTd,:) = (footy_ref(i))/2;
    COMYcenter_r(nT*(i-1)+nTd+1:nT*(i-1)+nT,:) =  footy_ref(i);
    else
    COMXcenter_r(nT*(i-1)+1:nT*(i-1)+nTd,:) = (footx_ref(i-1) + footx_ref(i))/2;
    COMXcenter_r(nT*(i-1)+nTd+1:nT*(i-1)+nT,:) =  footx_ref(i);
    COMYcenter_r(nT*(i-1)+1:nT*(i-1)+nTd,:) = (footy_ref(i-1) + footy_ref(i))/2;
    COMYcenter_r(nT*(i-1)+nTd+1:nT*(i-1)+nT,:) =  footy_ref(i);
    end
end

ZMPx_real = zeros(1,Nsum);  ZMPy_real = zeros(1,Nsum);
comx = zeros(1,Nsum);  comvx = zeros(1,Nsum);     comax = zeros(1,Nsum);
comy = zeros(1,Nsum);  comvy = zeros(1,Nsum);     comay = zeros(1,Nsum);
comz = zeros(1,Nsum);  comvz = zeros(1,Nsum);     comaz = zeros(1,Nsum);


%% CoM+angular state and control input
xk = zeros(3,Nsum);    x_vacc_k = zeros(1,Nsum);
yk = zeros(3,Nsum);    y_vacc_k = zeros(1,Nsum);
zk = zeros(3,Nsum);    z_vacc_k = zeros(1,Nsum);

yk(1,:) = footy_ref(1)*ones(1,Nsum);
zk(1,:) = hcom*ones(1,Nsum);
%% pamameters for MPC

A = [1,dt,dt^2/2; 0, 1,dt;0,0,1];
B = [dt^3/6; dt^2/2;dt];
C = [1,0,-hcom/g];
Cp = [1,0,0];
Cv = [0,1,0];
Ca = [0,0,1];


%%  vertical height range _initialize
%%%%fixed height
% Z_max = 0.03*ones(Nsum,1);
%% time-varying constraints
% Z_max = 0.000*ones(Nsum,1);
% for iiii = 1:Nsum
%     Z_max(iiii,:) = Z_max(iiii,:)-0.013*cos(2*pi/Tx(2)*(iiii-1)*dt);
% end
%% time-varying constraints+stage
% for iiii = 4*nT:7*nT
%     Z_max(iiii,:) = Z_max(iiii,:)-abs(0.013*cos(2*pi/Tx(2)*(iiii-1)*dt));
% end

Z_max = 0.1*ones(Nsum,1);

Z_min = -0.2;


%% footz reference :just the height of footz location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zsc = zeros(Nsum,1);
for i = 1:Nsum
    TPX = find(i*dt>=Tx);
    tpx1 = TPX(end);
    Zsc(i) = footz_ref(tpx1);
    
end

yk(1,:) = footy_ref(1)*ones(1,Nsum);
zk(1,:) = hcom*ones(1,Nsum);

%% predicitve model
% %%% for CoP
% Pzs = matrix_ps(A,Nh,C);
% Pzu = matrix_pu(A,B,Nh,C);
%%% for comy
Pps = matrix_ps(A,Nh,Cp);
Ppu = matrix_pu(A,B,Nh,Cp);
%% for comvy
Pvs = matrix_ps(A,Nh,Cv);
Pvu = matrix_pu(A,B,Nh,Cv);
%% for comay
Pas = matrix_ps(A,Nh,Ca);
Pau = matrix_pu(A,B,Nh,Ca);

footx_real = zeros(Tn,1);
footx_real_next = zeros(Nsum,1);
footy_real = zeros(Tn,1);
footy_real_next = zeros(Nsum,1);
footz_real = zeros(Nsum,1);
footz_real_next = zeros(Nsum,1);

Footx_max = 0.6;
Footx_min = 0.0;
Footy_max = 0.6;
Footy_min = 0.2;


%% robot parameters
%%% physical
M_f = 100;
rad = 0.4;
J_ini = M_f*rad^2;

%% constraints-boundaries parameters initilize
%%% ZMP BOUNDARIES-UPPER BOUNDARIES AND LOWER BOUNDARIES

%%% based on current support leg
ZMPx_ub = 0.08*ones(Nsum,1);     ZMPx_lb = -0.03*ones(Nsum,1);
ZMPy_ub = 0.033*ones(Nsum,1);    ZMPy_lb = -0.033*ones(Nsum,1);



%%% com-support range
comx_max = 0.06;
comx_min = -0.04;
comy_max = 0.6;
comy_min = 0.02;

%%% angle range
thetax_max = 13*pi/180;
thetax_min = -10*pi/180;

thetay_max = 10*pi/180;
thetay_min = -10*pi/180;


%%% torque range
torquex_max = 160;
torquex_min = -160;

torquey_max = 160;
torquey_min = -160;
torquex_max = torquex_max/J_ini;
torquex_min = torquex_min/J_ini;

torquey_max = torquey_max/J_ini;
torquey_min = torquey_min/J_ini;


%%% swing foot maximal velocity
footx_vmax = 3;
footx_vmin = -1;
footy_vmax = 1;
footy_vmin = -0.3;

%% solution preparation
%%%%%  optimization variables/state solution: be careful that the number of optimal variables is changing:
V_optimal = zeros(3*Nh+3*2,Nsum);
%%% n_vis flag
flag = zeros(1,Nsum);

% % Rx = 1; alphax = 100; beltax = 1000; gamax =10000000; %%%%% balance
% % Ry = 1; alphay = 100; beltay = 1000; gamay =10000000; %%%%% balance#
% Rx = 1; alphax = 10; beltax = 300000; gamax =30000000; %%%%% balance
% Ry = 1; alphay = 10; beltay = 100000; gamay =30000000; %%%%% balance
% Rz = 1; alphaz = 100; beltaz =10000000; gamaz =200; %%%%% balance
% Rthetax = 1; alphathetax = 10; beltathetax = 1000000;  %%%%% balance
% Rthetay = 1; alphathetay = 10;  beltathetay = 1000000;  %%%%% balance

% Rx = 1; alphax = 10; beltax = 300000; gamax =30000000; %%%%% balance
% Ry = 1; alphay = 10; beltay = 100000; gamay =30000000; %%%%% balance
% Rz = 1; alphaz = 10; beltaz = 100000000;gamaz =200; %%%%% balance
% Rthetax = 1; alphathetax = 1; beltathetax = 300000;  %%%%% balance
% Rthetay = 1; alphathetay = 1; beltathetay = 300000;  %%%%% balance

% Rx = 1; alphax = 10; beltax = 300000; gamax =10000000; %%%%% balance
% Ry = 1; alphay = 10; beltay = 100000; gamay =10000000; %%%%% balance
% Rz = 1; alphaz = 10; beltaz = 30000000;gamaz =200; %%%%% balance
% Rthetax = 1; alphathetax = 10; beltathetax = 300000;  %%%%% balance
% Rthetay = 1; alphathetay = 10; beltathetay = 100000;  %%%%% balance

% Rx = 1; alphax = 10; beltax = 300000; gamax =10000000; %%%%% balance
% Ry = 1; alphay = 10; beltay = 100000; gamay =10000000; %%%%% balance
% Rz = 1; alphaz = 10; beltaz = 30000000;gamaz =200; %%%%% balance
% Rthetax = 1; alphathetax = 10; beltathetax = 300000;  %%%%% balance
% Rthetay = 1; alphathetay = 10; beltathetay = 300000;  %%%%% balance

% Rx = 1; alphax = 10; beltax = 300000; gamax =10000000; %%%%% balance
% Ry = 1; alphay = 10; beltay = 100000; gamay =10000000; %%%%% balance
% Rz = 1; alphaz = 10;beltaz = 20000000;gamaz =200; %%%%% balance
% Rthetax = 1; alphathetax = 10; beltathetax = 300000;  %%%%% balance
% Rthetay = 1; alphathetay = 10; beltathetay = 230000;  %%%%% balance

Rx = 1; alphax = 10; beltax = 3000; gamax =10000000; %%%%% balance
Ry = 1; alphay = 10; beltay = 3000; gamay =10000000; %%%%% balance
Rz = 1; alphaz = 10; beltaz = 20000000;gamaz =200; %%%%% balance
Rthetax = 1; alphathetax = 10; beltathetax = 200000;  %%%%% balance
Rthetay = 1; alphathetay = 10; beltathetay = 200000;  %%%%% balance


%% test time consumption:
tcpu=zeros(1,Nsum-Nh-1);

%% predictive control_tracking with time-varying height
for i=1:1:Nsum-Nh
    t1=cputime;
    %%%% current period:
    bj1x = find(i*dt>=Tx);
    bjxx = bj1x(end); %%起始位置
    
    %% COM_center_ref = ZMP_center_ref = v_i*f + V_i*L_ref %%%% ref-foot
    %solve the following steps: 1.3s may couve 2 or three  steps, so  one/two following steps
    t_f = (i+1):1:(i+Nh);       %%预测时域长度
    bj1 = find(t_f(1)*dt>Tx);
    bjx1 = bj1(end);             %%预测开始周期
    bj2 = find(t_f(end)*dt>Tx);
    bjx2 = bj2(end);
    mx = bjx2-bjx1+1;        %%预测时域包含的周期
    bjx = bjx1:1:bjx2;
    tnx = zeros(mx,1);
    record_mx(:,i)=mx;
    
    
    
    for j =1:mx-1
        ccc = find(abs(Tx(bjx(j+1))-t_f*dt)<=0.001);  %%找到预测时域结束的周期
        tnx(j) = ccc(1);            %%0~C1 落在第一个周期，C1~C2落在第二个周期，C2~Nh落在第三个周期
    end
    %     v_i = zeros(Nh,1);
    %     if tnx(1) == nT+1          %%%% mx ==2; two following steps
    %         V_i = zeros(Nh,mx);
    %         for jj = 1:mx
    %             if jj==1
    %                 xxx1 =tnx(1)-1;
    %                 V_i(1:xxx1,jj)=ones(xxx1,1);
    %             else
    %                 xxx =Nh-tnx(1)+1;
    %                 V_i(tnx(1):end,jj)=ones(xxx,1);
    %             end
    %         end
    %     else
    %         xxx =tnx(1)-1;
    %         v_i(1:tnx(1)-1,:) = ones(xxx,1);
    %         V_i = zeros(Nh,mx-1);
    %         if mx ==2                              %%% one following steps
    %             V_i(tnx(1):Nh,1)=ones(Nh-tnx(1)+1,1);
    %         else                                    %%% mx==3, two following steps
    %             xxx1 = Nh-tnx(2)+1;
    %             V_i(tnx(1):tnx(2)-1,1)=ones(nT,1);
    %             V_i(tnx(2):Nh,2)=ones(xxx1,1);
    %         end
    %     end
    
    v_i = zeros(Nh,1);
    
    vv_i = zeros(Nh,mx+1);
 V_i = zeros(Nh,1);   
%     footx_p = footx_ref(bjx1:bjx2+1);
%     
%     if mx == 2        
%         C1 = find(abs(Tx(bjx2)-t_f*dt)<=0.001);
%         C2 = Nh-C1;
%         
%         V_i(C1+1:Nh,:)=1;
%         v_i(1:C1) = 1;
%         if C1<= n_t_s
%             vv_i(1:C1,2) = 1;
%             if C2 <=n_t_d
%                 vv_i(C1+1:Nh,2) = 1/2;
%                 vv_i(C1+1:Nh,3) = 1/2;
%             else
%                 vv_i(C1+1:C1+n_t_d,2) = 1/2;
%                 vv_i(C1+1:C1+n_t_d,3) = 1/2;
%                 
% %                 vv_i(C1+n_t_d+1:Nh,1) = 1/2;
%                 vv_i(C1+n_t_d+1:Nh,3) = 1;
%             end
%         else
%             vv_i(1:C1-n_t_s,1) = 1/2;
%             vv_i(1:C1-n_t_s,2) = 1/2;    %%DSP
%             
%             vv_i(C1-n_t_s+1:C1,2) = 1;   %%SSp
%             
%             if C2 <=n_t_d
%                 vv_i(C1+1:Nh,2) = 1/2;
%                 vv_i(C1+1:Nh,3) = 1/2;
%             else
%                 vv_i(C1+1:C1+n_t_d,2) = 1/2;
%                 vv_i(C1+1:C1+n_t_d,3) = 1/2;
% 
%                 vv_i(C1+n_t_d+1:Nh,3) = 1;
%             end
%         end
%     else
%         if mx == 3
%             C1 = find(abs(Tx(bjx1+1)-t_f*dt)<=0.001);
%             C3 = find(abs(Tx(bjx2+1)-t_f*dt)<=0.001);
%             v_i(1:C1) = 1;
%             V_i = zeros(Nh,2);
%             V_i(C1+1:C1+nT,1) = 1;
%             V_i(C1+nT+1:Nh,2) = 1;
%             if C1<= n_t_s 
%                 vv_i(1:C1,1) = 1;
%                 if C3 < n_t_d
%                     vv_i(C1+1:C1+nT,2) = 1;
%                     vv_i(C1+nT+1:end,3) = 1/2 ;
%                 elseif  C3 >n_t_d
%                     vv_i(C1+1:C1+nT,2) = 1;
%                     vv_i(C1+nT+1:C1+nT+n_t_d,3) = 1/2 ;
%                     vv_i(C1+nT+n_t_d+1:end,4) = 1/2 ;
%                 end       
%             else
%                 vv_i(1:C1-n_t_s,1) = 1/2;
%                 vv_i(C1-n_t_s+1:C1,2) = 1;
%                 if C3 <n_t_d
%                     vv_i(C1+1:C1+nT,3) = 1;
%                     vv_i(C1+nT+1:end,4) = 1/2 ;
%                 elseif  C3 >n_t_d
%                     vv_i(C1+1:C1+nT,3) = 1;
%                     vv_i(C1+nT+1:C1+nT+n_t_d,4) = 1/2 ;
%                     vv_i(C1+nT+n_t_d+1:end,5) = 1/2 ;
%                 end       
%             end
%         end
%     end
    
%     Center_x_r = vv_i*footx_p;
    
    [~,n_vis]=size(V_i);   %%% the number of following steps: one or two:    
    
    Lx_ref = zeros(n_vis,1);
    Ly_ref = zeros(n_vis,1);
    Lz_ref = zeros(n_vis,1);
    if n_vis == 1      %%%只预测一步
        Lx_ref(1) = footx_ref(bjx2,:);
        Ly_ref(1) = footy_ref(bjx2,:);
        Lz_ref(1) = footz_ref(bjx2,:);
    else                %%%预测两步
        Lx_ref(1) = footx_ref(bjx2-1,:);
        Lx_ref(n_vis) = footx_ref(bjx2,:);
        Ly_ref(1) = footy_ref(bjx2-1,:);
        Ly_ref(n_vis) = footy_ref(bjx2,:);
        Lz_ref(1) = footz_ref(bjx2-1,:);
        Lz_ref(n_vis) = footz_ref(bjx2,:);
    end
        
    flag(i) = n_vis;    
%     Nt = 3*Nh+3*n_vis;   %%%X_j, Y_j Z_j  thelt_x, thelt_y,  footx, footy,footz
    Nt = 3*Nh;
    %% hot start
    V_ini = zeros(Nt,1);
    V_ini2 = zeros(3*Nh,1);
    
%     COMx_center_ref = footx_r(i+1:i+Nh,:);
%     COMy_center_ref = footy_r(i+1:i+Nh,:);
    COMx_center_ref = COMXcenter_r(i+1:i+Nh,:);
    COMy_center_ref = COMYcenter_r(i+1:i+Nh,:);
    %     footx_ref =
    
    %%%v_i记录当前支撑足位置，用于预测时域内初始足部参考位置表达
    %%%Lx_ref为预测时域内后一个/两个周期参考位置，V_i为取出矩阵
    
    COMz_center_ref = Zsc(i+1:i+Nh,:)+hcom;
    
    COMx_ref(i,:) = COMx_center_ref(1,:);
    COMy_ref(i,:) = COMy_center_ref(1,:);
    if i == Nsum-Nh-1
        COMx_ref(i+1:i+Nh,:) = COMx_center_ref;
        COMy_ref(i+1:i+Nh,:) = COMy_center_ref;
    end
    
    %% SEQUENCE QUADARTIC PROGRAMMING
    for xxx = 1:3      %%% definition of Ns:
        %% SQP MODELS
        % model coefficient matrix
        %%%optimization number is = Nh*3+ n_vis*3
        %%%%% optimal programme formulation:
        %%%%% x0 = xk(:,i);y0 = yk(:,i);z0 = zk(:,i);
        WX = Rx/2*eye(Nh) + alphax/2*(Pvu'*Pvu) + beltax/2*(Ppu'*Ppu);
        WY = Ry/2*eye(Nh) + alphay/2*(Pvu'*Pvu) + beltay/2*(Ppu'*Ppu);
        WZ = Rz/2*eye(Nh) + alphaz/2*(Pvu'*Pvu) + beltaz/2*(Ppu'*Ppu);
        
        PHIX = gamax/2*eye(n_vis);
        PHIY = gamay/2*eye(n_vis);
        PHIZ = gamaz/2*eye(n_vis);
        
        Q_goal = blkdiag(WX,WY,WZ,PHIX,PHIY,PHIZ);
        
        
        q_goal = [alphax*Pvu'*Pvs*xk(:,i)+beltax*Ppu'*Pps*xk(:,i)-beltax*Ppu'* COMx_center_ref;...
            alphay*Pvu'*Pvs*yk(:,i)+beltay*Ppu'*Pps*yk(:,i)-beltay*Ppu'* COMy_center_ref;...
            alphaz*Pvu'*Pvs*zk(:,i)+beltaz*Ppu'*Pps*zk(:,i)-beltaz*Ppu'* COMz_center_ref;...
            -gamax*Lx_ref;...
            -gamay*Ly_ref;...
            -gamaz*Lz_ref];
        
        
%         Q_goal1 = 2 * Q_goal;
%         q_goal1 = 2 * Q_goal * V_ini + q_goal;
        
        
        Q = blkdiag(WX,WY,WZ);
        q_goal22 = [alphax*Pvu'*Pvs*xk(:,i)+beltax*Ppu'*Pps*xk(:,i)-beltax*Ppu'* COMx_center_ref;...
                    alphay*Pvu'*Pvs*yk(:,i)+beltay*Ppu'*Pps*yk(:,i)-beltay*Ppu'* COMy_center_ref;...
                    alphaz*Pvu'*Pvs*zk(:,i)+beltaz*Ppu'*Pps*zk(:,i)-beltaz*Ppu'* COMz_center_ref];
        
        Q_goal2 = 2 * Q;
        q_goal2 = 2 * Q * V_ini2 + q_goal22;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% boundary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Sjx = zeros(Nh,Nt);  Sjy = zeros(Nh,Nt); Sjz = zeros(Nh,Nt);
        Sfx = zeros(n_vis,Nt);  Sfy = zeros(n_vis,Nt); Sfz = zeros(n_vis,Nt);
        Sjx(:,1:Nh) = eye(Nh);                 Sjy(:,Nh+1:2*Nh) = eye(Nh);                    Sjz(:,2*Nh+1:3*Nh) = eye(Nh);
        Sfx(:,3*Nh+1:3*Nh+n_vis) = eye(n_vis); Sfy(:,3*Nh+n_vis+1:3*Nh+2*n_vis) = eye(n_vis); Sfz(:,3*Nh+2*n_vis+1:3*Nh+3*n_vis) = eye(n_vis);
        
        
        % ZMP boundary preparation:
        H_q_upx = zeros(Nh,Nt);
        F_zmp_upx = zeros(Nh,1);
        H_q_lowx = zeros(Nh,Nt);
        F_zmp_lowx = zeros(Nh,1);
        H_q_upy = zeros(Nh,Nt);
        F_zmp_upy = zeros(Nh,1);
        H_q_lowy = zeros(Nh,Nt);
        F_zmp_lowy = zeros(Nh,1);        
        
        % CoM height boundary perparation:
        H_h_upz = zeros(Nh,Nt);
        F_h_upz = zeros(Nh,1);
        H_h_lowz = zeros(Nh,Nt);
        F_h_lowz = zeros(Nh,1);
        
        % CoM height accereation boundary preparation:
        H_hacc_lowz = zeros(Nh,Nt);
        F_hacc_lowz = zeros(Nh,1);
        
        
        %         Si*V_i*Sfx* W_jerk
        
        for j = 1:Nh
            Si = zeros(1,Nh); Si(j) = 1;
            %% ZMP constraints
            %% x-ZMP upper boundary
            phi_i_x_up =  Sjx'*Ppu'*Si'*Si*Pau*Sjz  ...                
                - Sjx'*Pau'*Si'*Si*Ppu*Sjz ;
            phi_i_x_up = M_f*(phi_i_x_up+phi_i_x_up')/2;
            
            zmpxxx_ref = Si*COMx_center_ref;
            zmpyyy_ref = Si*COMy_center_ref;
            

            p_i_x_t_up = (xk(:,i)'*Pps'*Si'*Si*Pau*Sjz  + zk(:,i)'*Pas'*Si'*Si*Ppu*Sjx  ...
                +  g*Si*Ppu*Sjx ...
                - (zk(:,i)'*Pps'*Si'*Si*Pau*Sjx  + xk(:,i)'*Pas'*Si'*Si*Ppu*Sjz) ...
                + Zsc(i+j,:)*Si*Pau*Sjx ...
                -zmpxxx_ref*Si*Pau*Sjz ...
                - ZMPx_ub(i+j)*Si*Pau*Sjz)';            
            p_i_x_t_up = M_f*p_i_x_t_up ;
            
            del_i_x_up = xk(:,i)'*Pps'*Si'*(zk(:,i)'*Pas'*Si')' ...
                + g*Si*Pps*xk(:,i) ...
                - xk(:,i)'*Pas'*Si'*(zk(:,i)'*Pps'*Si')' ...
                + xk(:,i)'*Pas'*Si'*  Zsc(i+j,:)...
                - g* zmpxxx_ref -zmpxxx_ref*Si*Pas*zk(:,i)...
                - ZMPx_ub(i+j)*Si*Pas*zk(:,i) ...
                - g*ZMPx_ub(i+j);
            del_i_x_up = M_f*del_i_x_up ;
            
            
            H_q_upx(j,:) = (2*phi_i_x_up*V_ini2 + p_i_x_t_up)';
            F_zmp_upx(j,:) = -(V_ini2'*phi_i_x_up*V_ini2 + p_i_x_t_up'*V_ini2+ del_i_x_up);
            
            
            %% x-ZMP low boundary
            
            phi_i_x_low =  Sjx'*Ppu'*Si'*Si*Pau*Sjz - Sjx'*Pau'*Si'*Si*Ppu*Sjz ;            
            phi_i_x_low = M_f*(phi_i_x_low+phi_i_x_low')/2;
            
            p_i_x_t_low = (xk(:,i)'*Pps'*Si'*Si*Pau*Sjz  + zk(:,i)'*Pas'*Si'*Si*Ppu*Sjx  ...
                +  g*Si*Ppu*Sjx ...
                - (zk(:,i)'*Pps'*Si'*Si*Pau*Sjx  + xk(:,i)'*Pas'*Si'*Si*Ppu*Sjz) ...
                + Zsc(i+j,:)*Si*Pau*Sjx ...
                -zmpxxx_ref*Si*Pau*Sjz ...
                - ZMPx_ub(i+j)*Si*Pau*Sjz)';            
            p_i_x_t_low = M_f*p_i_x_t_low ;
  
            del_i_x_low = xk(:,i)'*Pps'*Si'*(zk(:,i)'*Pas'*Si')' ...
                + g*Si*Pps*xk(:,i) ...
                - xk(:,i)'*Pas'*Si'*(zk(:,i)'*Pps'*Si')' ...
                + xk(:,i)'*Pas'*Si'*  Zsc(i+j,:)...
                - g* zmpxxx_ref -zmpxxx_ref*Si*Pas*zk(:,i)...
                - ZMPx_lb(i+j)*Si*Pas*zk(:,i) ...
                - g*ZMPx_lb(i+j);      
            del_i_x_low = M_f*del_i_x_low ;
            
            
            H_q_lowx(j,:) = -(2*phi_i_x_low*V_ini2 + p_i_x_t_low)';
            F_zmp_lowx(j,:) = (V_ini2'*phi_i_x_low*V_ini2 + p_i_x_t_low'*V_ini2+ del_i_x_low);
            
            %% y-ZMP upper boundary
            phi_i_y_up =  Sjy'*Ppu'*Si'*Si*Pau*Sjz  ...                
                - Sjy'*Pau'*Si'*Si*Ppu*Sjz ;
            phi_i_y_up =  M_f*(phi_i_y_up+phi_i_y_up')/2;
            
            
            p_i_y_t_up = (yk(:,i)'*Pps'*Si'*Si*Pau*Sjz  + zk(:,i)'*Pas'*Si'*Si*Ppu*Sjy  ...
                +  g*Si*Ppu*Sjy ...
                - (zk(:,i)'*Pps'*Si'*Si*Pau*Sjy  + yk(:,i)'*Pas'*Si'*Si*Ppu*Sjz) ...
                + Zsc(i+j,:)*Si*Pau*Sjy ...
                -zmpyyy_ref*Si*Pau*Sjz ...
                - ZMPy_ub(i+j)*Si*Pau*Sjz)';            
            p_i_y_t_up = M_f*p_i_y_t_up ;
                      
            
            del_i_y_up = yk(:,i)'*Pps'*Si'*(zk(:,i)'*Pas'*Si')' ...
                + g*Si*Pps*yk(:,i) ...
                - yk(:,i)'*Pas'*Si'*(zk(:,i)'*Pps'*Si')' ...
                + yk(:,i)'*Pas'*Si'*  Zsc(i+j,:)...
                - g* zmpyyy_ref -zmpyyy_ref*Si*Pas*zk(:,i)...
                - ZMPy_ub(i+j)*Si*Pas*zk(:,i) ...
                - g*ZMPy_ub(i+j);         
            del_i_y_up = M_f*del_i_y_up ;
            
            
            H_q_upy(j,:) = (2*phi_i_y_up*V_ini2 + p_i_y_t_up)';
            F_zmp_upy(j,:) = -(V_ini2'*phi_i_y_up*V_ini2 + p_i_y_t_up'*V_ini2+ del_i_y_up);            
            
            %% y-ZMP lower boundary

            phi_i_y_low =  Sjy'*Ppu'*Si'*Si*Pau*Sjz - Sjy'*Pau'*Si'*Si*Ppu*Sjz ;
            
            phi_i_y_low =  M_f*(phi_i_y_low+phi_i_y_low')/2;
            
            p_i_y_t_low = (yk(:,i)'*Pps'*Si'*Si*Pau*Sjz  + zk(:,i)'*Pas'*Si'*Si*Ppu*Sjy  ...
                +  g*Si*Ppu*Sjy ...
                - (zk(:,i)'*Pps'*Si'*Si*Pau*Sjy  + yk(:,i)'*Pas'*Si'*Si*Ppu*Sjz) ...
                + Zsc(i+j,:)*Si*Pau*Sjy ...
                -zmpyyy_ref*Si*Pau*Sjz ...
                - ZMPy_lb(i+j)*Si*Pau*Sjz)';            
            p_i_y_t_low = M_f*p_i_y_t_low ;
       
            del_i_y_low = yk(:,i)'*Pps'*Si'*(zk(:,i)'*Pas'*Si')' ...
                + g*Si*Pps*yk(:,i) ...
                - yk(:,i)'*Pas'*Si'*(zk(:,i)'*Pps'*Si')' ...
                + yk(:,i)'*Pas'*Si'*  Zsc(i+j,:)...
                - g* zmpyyy_ref -zmpyyy_ref*Si*Pas*zk(:,i)...
                - ZMPy_lb(i+j)*Si*Pas*zk(:,i) ...
                - g*ZMPy_lb(i+j);                       
            del_i_y_low = M_f*del_i_y_low ;
            
            H_q_lowy(j,:) = -(2*phi_i_y_low*V_ini2 + p_i_y_t_low)';
            F_zmp_lowy(j,:) = (V_ini2'*phi_i_y_low*V_ini2 + p_i_y_t_low'*V_ini2+ del_i_y_low);
            
            %% body height constraints
%             P_footz_up = Si*Ppu*Sjz;
%             delta_footz_up = Si*Pps*zk(:,i) - COMz_center_ref(j) - Z_max(i+j,:);
%             H_h_upz(j,:) = P_footz_up;
%             F_h_upz(j,:) = -(P_footz_up*V_ini + delta_footz_up);
%             
%             P_footz_up = Si*Ppu*Sjz;
%             delta_footz_lowp = Si*Pps*zk(:,i) - COMz_center_ref(j) - Z_min;
%             H_h_lowz(j,:)  = -P_footz_up;
%             F_h_lowz(j,:)  = (P_footz_up*V_ini + delta_footz_lowp);
%             
            
            
            %% body height acceleration
%             P_footzacc_up = Si*Pau*Sjz;
%             delta_footzacc_up = Si*Pas*zk(:,i) -(-g);
%             H_hacc_lowz(j,:)  = -P_footzacc_up;
%             F_hacc_lowz(j,:)  = (P_footzacc_up*V_ini + delta_footzacc_up);
%             
        end
        
        %% CoM-supportl leg constraints
        %%% next time moment support leg: s1*(v_i*fx + V_i*Sfw*V_ini);
        
        
        
        %% Foot vertical loction-equality constraints
%         H_q_footz = Sfz;
%         F_footz = -(Sfz*V_ini-Lz_ref);
        
        %% fixed height
%         h_h = Ppu*Sjz;
%         hhhx = -(Ppu*Sjz*V_ini+Pps*zk(:,i)-hcom*ones(Nh,1));
%         
%         H_q_footz1 = [H_q_footz;h_h];
%         F_footz1 = [F_footz;hhhx];
        
        %% fixed inclined angle

        A_q2 = [H_q_upx;  H_q_lowx;   H_q_upy;   H_q_lowy];        
        b_q2 = [F_zmp_upx;F_zmp_lowx; F_zmp_upy; F_zmp_lowy];        
                       
        X_jerk = quadprog(Q_goal2,q_goal2,A_q2,b_q2);
        
        V_ini2 = V_ini2 + X_jerk;
    end
    
    
    %%
    x_jerk(i) = V_ini2(1);
    y_jerk(i) = V_ini2(Nh+1);
    z_jerk(i) = V_ini2(2*Nh+1);
    xk(:,i+1) = A*xk(:,i)+B*x_jerk(i);
    yk(:,i+1) = A*yk(:,i)+B*y_jerk(i);
    zk(:,i+1) = A*zk(:,i)+B*z_jerk(i);
    if i == Nsum-Nh-1
        xxk = Pps*xk(:,i)+Ppu*Sjx*V_ini2;
        vxk = Pvs*xk(:,i)+Pvu*Sjx*V_ini2;
        axk = Pas*xk(:,i)+Pau*Sjx*V_ini2;
        xxxk = [xxk' ;vxk' ;axk'];
        xk(:,i+1:i+Nh) = xxxk;
        
        yyk = Pps*yk(:,i)+Ppu*Sjy*V_ini2;
        vyk = Pvs*yk(:,i)+Pvu*Sjy*V_ini2;
        ayk = Pas*yk(:,i)+Pau*Sjy*V_ini2;
        yyyk = [yyk' ;vyk' ;ayk'];
        yk(:,i+1:i+Nh) = yyyk;    
        xk(:,end) = xk(:,end-1);
        yk(:,end) = yk(:,end-1);
        
    end
    
    ZMPx_real(i) = xk(1,i)-(zk(1,i)/(g+zk(3,i)))*xk(3,i);
    ZMPy_real(i) = yk(1,i)-(zk(1,i)/(g+zk(3,i)))*yk(3,i);
    t2=cputime;
    tcpu(i)=t2-t1;
end


end

%% matrix solution
function [y]=matrix_ps(A,N,C)
y = zeros(N,3);
for i =1:N
    y(i,:) = C*A^i;
end

end
function [y]=matrix_pu(A,B,N,C)
y = zeros(N,N);
for i = 1:N
    for j = 1:i
        y(i,j) = C*A^(i-j)*B;
    end
end


end


