## -*- encoding: UTF-8 -*-
import time
import numpy as np
import math
import string
import scipy.interpolate as itp
from naoqi import ALProxy
import argparse
import matplotlib.pyplot as plt

##---本地子函数简记------
local_sin = math.sin
local_cos = math.cos
local_dot = np.dot
local_floor = math.floor
uni_matrix = np.eye(3)
local_array = np.array
local_zeros = np.zeros
local_sinh = math.sinh
local_cosh = math.cosh
local_tanh = math.tanh
local_ones = np.ones




##-----------全局变量声明----------
##----默认浮点数为float64
##这里强行调整为float32,同时去除掉默写不必要的全局变量声明。
mylink = []
Link = []
eps = 2.2204e-008

NNx = 40
##%总时间
dtx =  0.02
Tx = 2
Tsupx = 1.6
Tdoux = 0.4
stepwidthx = 100.0
stepxx = 40.0
gx = 9800
Zcx = 310.0

flagxn =[]
total_n = int(local_floor(NNx * Tx/dtx))

Input = open(r'E:\python_unknown_slope\nao_20160504straight5_mod_2s_100_40_angle.txt','r')
data = []
for i in range(0,10):
    a = Input.readline()
    a = string.split(a)
    data.append(a)
data = map(map,[float,float,float,float,float,float,float,float,float,float], data)
datax = local_array(data).T
data_angle0=data[0];  data_angle1=data[1]; data_angle2=data[2]; data_angle3=data[3];
data_angle4=data[4];  data_angle5=data[5]; data_angle6=data[6]; data_angle7=data[7];
data_angle8=data[8];  data_angle9=data[9];
data_angle0 =local_array(data_angle0);  data_angle1 =local_array(data_angle1);  data_angle2 =local_array(data_angle2);
data_angle3 =local_array(data_angle3);  data_angle4 =local_array(data_angle4);  data_angle5 =local_array(data_angle5);
data_angle6 =local_array(data_angle6);  data_angle7 =local_array(data_angle7);  data_angle8 =local_array(data_angle8);           
data_angle9 =local_array(data_angle9);

data =np.vstack((data_angle0,data_angle1,data_angle2,data_angle3,data_angle4,data_angle5,
                        data_angle6,data_angle7,data_angle8,data_angle9))
data_anglex = np.float32(data.T)

data = []
Input = open(r'E:\python_unknown_slope\nao_20160504straight5_mod_preview_para.txt','r')
for i in range(0,1):
    a = Input.readline()
    a = string.split(a)
    data.append(a)
data = map(map,[float], data)
F_pre = np.float32(np.reshape(local_array(data),[100,1]))
sumde1_p = local_zeros([1,3])

data = []
Input = open(r'E:\python_unknown_slope\nao_20160504straight5_mod_2s_100_40_zmp.txt','r')
for i in range(0,3):
    a = Input.readline()
    a = string.split(a)
    data.append(a)
data = map(map,[float,float,float], data)
data = np.vstack((local_array(data[0]),local_array(data[1]),local_array(data[2])))
data = np.float32(data.T)
zmp_ref = data[0:total_n,:]
zmp_offline = data[0:total_n,:]

zmp_real = np.float32(local_zeros([total_n,3]))
com_real = np.float32(local_zeros([total_n,3]))

data = []
Input = open(r'E:\python_unknown_slope\nao_20160504straight5_mod_2s_100_40_bodyp2.txt','r')
for i in range(0,3):
    a = Input.readline()
    a = string.split(a)
    data.append(a)
data = map(map,[float,float,float], data)
data =np.vstack((local_array(data[0]),local_array(data[1]),local_array(data[2])))
data =data.T
com_ref =np.float32(data[0:total_n,:])

zmp_ref[:,2] = com_ref[:,2]
zmp_offline[:,2] = com_ref[:,2]


t_land = 0.1*local_ones([NNx,1])
t_off = 0.1*local_ones([NNx,1])
Mea = local_zeros([total_n,21])
angle = local_zeros([total_n,10])
flag = local_zeros([NNx,1])
wf = local_zeros([total_n,3])
wf[:,1] = 55 * local_ones([total_n])
zmp_g = local_zeros([total_n,3])
com_g = local_zeros([total_n,3])
rfoot_g = local_zeros([total_n,3])
lfoot_g = local_zeros([total_n,3])
n_period = local_zeros([1,NNx])
n_period_N  = local_zeros([1,NNx])

flagx = local_zeros([total_n,2])
T_initial = 6
DSP = np.float32(Tdoux * local_ones([NNx,1]))
SSP = np.float32(Tsupx * local_ones([NNx,1]))
R_sensor = local_array([[0.07025,0.0231,-45.19],[0.07025,-0.0299,-45.19],[-0.03025,0.0191,-45.19],[-0.02965,-0.0299,-45.19]])
L_sensor = local_array([[0.07025,0.0299,-45.19],[0.07025,-0.0231,-45.19],[-0.03025,0.0299,-45.19],[-0.02965,-0.0191,-45.19]])

LZMPy0x = 50.0
tN1x = int(local_floor((2*Tx-dtx)*100/dtx/100))                                           

##    %两个周期，四段
x_init = local_zeros([4,1])                  
x_end = local_zeros([4,1])
y_init = local_zeros([4,1])
y_end = local_zeros([4,1])
z_init = local_zeros([4,1])
z_end = local_zeros([4,1])
Vx_init = local_zeros([4,1])
Vx_end = local_zeros([4,1])
Vy_init = local_zeros([4,1])
Vy_end = local_zeros([4,1])
Vz_init = local_zeros([4,1])
Vz_end = local_zeros([4,1])

##    %-----各个时刻的ZMP位置（相对于全局坐标系）====
PX1= local_zeros([tN1x,1])
PY1 = local_zeros([tN1x,1])
PZ1 = Zcx *local_ones([tN1x,1])

##    %步行参数序列,
Sx = stepxx * local_ones([NNx,1])
##    %---第一步迈右脚，步长序列---
Sx[0,0] =20.0
Sx[1,0] =30.0

##    %---第一步迈右脚，步宽序列----
Sy = stepwidthx * local_ones([NNx,1])
Sy[0,0] = LZMPy0x + stepwidthx/2
widtxh_iNNxer = 5 
Sy[1,0] = Sy[1,0]-widtxh_iNNxer
Sy[2:,:] = (Sy[1,0]-widtxh_iNNxer) * local_ones([NNx-2,1])

Sz = Zcx * local_ones([NNx,1])
z_h = local_zeros([NNx,1])


slope_angle = local_zeros([NNx,1])
for i in range(0,NNx):
    z_h[i,0] = -Sx[i,0] * math.tan(slope_angle[i,0])

yaw_ref = local_zeros([total_n,1])
f_angle_roll = local_zeros([total_n,1]) 
f_angle_pitch =local_zeros([total_n,1]) 
f_angle_yaw =local_zeros([total_n,1])
body_Lroll =[]                                                      
body_Lpitch =[]                                                     
body_Rroll =[]                                                       
body_Rpitch =[]
body_yaw =[]
pitch_err_intergral =[]                                              
roll_err_intergral =[]
yaw_err_intergral =[]

##    %规划的落脚点位置，也就是各个时期局部坐标系的初始位置（相对于全局坐标系）-----
footx = local_zeros([NNx,1])
footy = local_zeros([NNx,1])
footz = local_zeros([NNx,1])
##    %初始化
footx[0,0] = 0
footy[0,0] = LZMPy0x
footx[1,0] = footx[0,0] + Sx[0,0]
footy[1,0] = footy[0,0] - Sy[0,0]
footx[2,0] = footx[1,0] + Sx[1,0]
footy[2,0] = footy[1,0] + Sy[1,0]
footx[3,0] = footx[2,0] + Sx[2,0]
footy[3,0] = footy[2,0] - Sy[2,0]

footy2 = local_zeros([NNx,1])

##计算 n_period
for i in range(0,NNx):
    if i==0:
        n_period[:,i] = T_initial
    else:
        n_period[:,i]= n_period[:,i-1]+Tx


for i in range(0,NNx):
    if i==0:
        n_period_N[:,i] = int(local_floor(T_initial/dtx))
    else:
        n_period_N[:,i]= n_period_N[:,i-1]+int(local_floor(Tx/dtx))
  

################===================================================#######################
################===================================================#######################
################===================================================#######################
###                左右足部踝关节俯仰自由度初始化，其实是由footpR函数生成的
def posRef(): 
    class posRef(object):
                def __init__(self,p,R,v,w):           
                    self.p = p
                    self.R = R
                    self.v = v
                    self.w = w
                    
    p_init = local_array([[0],[0],[0]])
    R_init = np.eye(3)
    v_init = local_array([[0],[0],[0]])
    w_init = local_array([[0],[0],[0]])

    R_RefpR = posRef(p_init,R_init,v_init,w_init)
    L_RefpR = posRef(p_init,R_init,v_init,w_init)    
    
    return R_RefpR,L_RefpR



def Initial_once():
    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones  
    class Link(object):
                def __init__(self,name,brotherID,childID,motherID,c,m,dq,I,v,w,q,a,b,p,R):
                        self.name = name
                        self.brotherID = brotherID
                        self.childID = childID
                        self.motherID = motherID
                        self.c = c
                        self.m = m
                        self.dq = dq
                        self.I = I
                        self.v = v
                        self.w = w
                        self.q = q
                        self.a = a
                        self.b = b
                        self.p = p
                        self.R = R
    
    Name=['Body','Rleg_hip_r','Rleg_hip_p','Rleg_knee','Rleg_ank_p','Rleg_ank_r','Lleg_hip_r','Lleg_hip_p','Lleg_knee','Lleg_ank_p','Lleg_ank_r']
    BrotherID=[0,7,0,0,0,0,0,0,0,0,0]
    ChildID=[2,3,4,5,6,0,8,9,10,11,0]
    MotherID=[0,1,2,3,4,5,1,7,8,9,10]
#//    多生成一个，保持结构体一致
    a_Init = local_array([[[1], [0], [0]],[[1], [0], [0]],[[0], [1], [0]],[[0], [1], [0]],
                       [[0],[1], [0]],[[1], [0], [0]],[[1], [0], [0]],[[0], [1], [0]],
                       [[0], [1], [0]],[[0], [1], [0]],[[1], [0], [0]]])
    b_Init = local_array([[[0], [0], [0]],[[0], [-50], [-85]],[[0], [0], [0]],[[0], [0], [-100]],
                       [[0], [0], [-102.9]],[[0], [0], [0]],[[0], [50], [-85]],[[0], [0], [0]],
                       [[0], [0], [-100]],[[0], [0], [-102.9]],[[0], [0], [0]]])
    
    c_Init = local_array([[-1.44,-1.29,54.92],[-15.49,-0.29,-5.15],[1.38,-2.21,-53.73],
                       [4.53,-2.25,-49.36],[0.45,-0.29,6.85],[25.42,-3.3,-32.39],
                       [-15.49,0.29,-5.15],[1.38,2.21,-53.73],[4.53,2.25,-49.36],
                       [0.45,0.29,6.85],[25.42,3.3,-32.39]])
    
    m_Init=[2.841,0.141,0.390,0.301,0.134,0.172,0.141,0.390,0.301,0.134,0.172]
    dq_Init=[0,0,0,0,0,0,0,0,0,0,0]
    I1 = local_array([[13953.66,6.19,-198.09],[6.19,13318.88,-196.16],[-198.09,-196.16,2682.42]])
    I2 = local_array([[2.76,-0.02,-4.11],[-0.02,9.83,0.00],[-4.11,0.00,8.81]])
    I3 = local_array([[1637.48,-0.92,85.88],[-0.92,1592.21,-39.18],[85.88,-39.18,303.98]])
    I4 = local_array([[1182.83,-0.90,28.00],[-0.90,1128.28,-38.48],[28.00,-38.48,191.45]])
    I5 = local_array([[38.51,-0.06,3.87],[-0.06,74.31,-0.00],[3.87,-0.00,54.91]])
    I6 = local_array([[269.30,5.88,139.13],[5.88,643.47,-18.85],[139.13,-18.85,525.03]])
    I7 = local_array([[2.76,-0.02,-4.08],[-0.02,9.83,-0.00],[-4.08,-0.00,8.81]])
    I8 = local_array([[1637.20,-0.92,85.31],[-0.92,1591.07,-38.36],[85.31,-38.36,303.74]])
    I9 = local_array([[1182.08,0.63,36.50],[0.63,1128.65,39.50],[36.50,39.50,193.22]])
    I10 = local_array([[38.51,-0.03,3.86],[-0.03,74.27,-0.02],[3.86,-0.02,54.87]])
    I11 = local_array([[269.44,-5.70,139.38],[-5.70,644.43,18.74],[139.38,18.74,525.76]])
    I_Init = [I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11]
    vw_Init = local_zeros((3,1))

    q = local_array([[0],[-0.024654],[-0.0282059],[0.0556197],[-0.0273161],[0.024654],
                  [0.024654],[-0.0282059],[0.0556197],[-0.0273161],[-0.024654]])
    p_Init=local_array([[[0], [0], [333.09]],[[0], [0], [0]],[[0], [0], [0]],[[0], [0], [0]],[[0], [0], [0]],
                     [[0], [0], [0]],[[0], [0], [0]],[[0], [0], [0]],[[0], [0], [0]],[[0], [0], [0]],[[0], [0], [0]]])
 
    R_Init=np.eye(3)
    global mylink
    for i in range(0,11):    
        link = Link(Name[i],BrotherID[i],ChildID[i],MotherID[i],c_Init[i],
                m_Init[i],dq_Init[i],I_Init[i],vw_Init,vw_Init,
               q[i],a_Init[i],b_Init[i],p_Init[i],R_Init)
        mylink.append(link)

##  x作为标志位，x为0时存在等于号，等于1的时候绝对大于。
def find0(arr,obj,x,nn):
    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones     
    tem = []
    j=0
    if x == 1:      
        for i in range(0,nn):
            if ((arr[0][i]-obj)>0.0000001):            
                tem.append(i)
                return tem        
    else:
        if x ==0:
            for i in range(0,nn):
                if arr[0][i]>=obj:            
                    tem.append(i)
                    tem.append(1)
                    return tem                
        else:
            if x == -1:
                for i in range(0,nn):
                    if abs(arr[0][i]-obj)<=0.0001:            
                        tem.append(i)
                return tem   
            else:               
                for i in range(0,nn):
                    if arr[0][i]>=obj:            
                        tem.append(i)
                return tem            
            

##  Rodrigure公式
##  符号计算法写出结果
def Rodrigues_x_once(a,q,local_sin,local_cos,local_dot,uni_matrix):
    a_skew = np.array([[0,  -a[2,0], a[1,0]],
                       [a[2,0],  0,  -a[0,0]],
                       [-a[1,0], a[0,0], 0]])
    a_skewpu = local_dot(a_skew,a_skew)
    
    ##已经有相应的算法了
    sq = local_sin(q)
    cq = 1 - local_cos(q)
    factor2 = a_skew * sq
    factor3 = a_skewpu * cq
    e_out = uni_matrix+ factor2 +factor3    
##  完全写成符号表达，结果时间反而延长      
    return e_out

###=========================正运动学公式
def Forward_Kinematics_x(j,local_sin,local_cos,local_dot,uni_matrix,mylink1):
    if j == 0:
        return
    if j != 1:        
        mylink1x = mylink1[j-1]
        i = mylink1x.motherID
        mylink2x = mylink1[i-1]             
        mylink1x.p = mylink2x.p+ local_dot(mylink2x.R, mylink1x.b)        
        Rod = Rodrigues_x_once(mylink1x.a, mylink1x.q,local_sin,local_cos,local_dot,uni_matrix)        
        mylink1x.R = local_dot(mylink2x.R ,Rod)        
        mylink1[j-1] = mylink1x     
       
    Forward_Kinematics_x(mylink1x.brotherID,local_sin,local_cos,local_dot,uni_matrix,mylink1)
    Forward_Kinematics_x(mylink1x.childID,local_sin,local_cos,local_dot,uni_matrix,mylink1)




####==========================================================================================
###========================反馈量计算   

def calCom(mylink1):
    #该程序返回总体质心，格式为3*N
    global Link
    global mc
    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones 
    
    Link = mylink1
    
    def calmc(j):
        if j == 0:
            mc = 0
            mc = local_zeros([3,1])
        else:
            i = j-1
            cx = np.reshape(Link[i].c,[3,1])
            Rx = Link[i].R
            mc = Link[i].p + local_dot(Rx,cx)
            mx = Link[i].m             
            mc = mx * mc
            mc = mc + calmc(Link[i].brotherID)+calmc(Link[i].childID)
        return mc
    
    def TotalMass(j):
        if j == 0:
            M = 0
        else:
            i = j-1
            m = Link[i].m
            M = m + TotalMass(Link[i].brotherID)+TotalMass(Link[i].childID)
        return M
    
    Mc = calmc(1)    #1->j
    com = Mc/TotalMass(1)
    return com

def Rea1_com(mylink1):
    Forward_Kinematics_x(2,local_sin,local_cos,local_dot,uni_matrix,mylink1)
    real_com=calCom(mylink1)
    return real_com  

    
def inverse_Eular(Rrpy):
    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones     
    y=-math.asin(Rrpy[2,0]);
    try:
        x=math.asin(Rrpy[2,1]/math.cos(y))
    except ZeroDivisionError as e:
        x=math.atan(Rrpy[2,1]/Rrpy[2,2])
    
    try:        
        z=math.asin(Rrpy[1,0]/math.cos(y))
    except ZeroDivisionError as e:
        z=math.atan(Rrpy[1,0]/Rrpy[0,0])
    return x,y,z
    
def ground_inclination_estimation(pick,xtn,mea,R_lf,L_lf,Rverti,Lverti,Lf,Rf):
    global Link
    global f_angle_roll
    global f_angle_pitch
    global f_angle_yaw
    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones 
    slope_estimation = local_zeros(3)
    ##首先从足部的旋转矩阵求解足部平板与地面的夹角。（目前不考虑绕z轴的偏转。）
    ##足部角度估计：
    R_R = Link[5].R
    L_R = Link[10].R
    [R_roll,R_pitch,R_yaw] = inverse_Eular(R_R)
    [L_roll,L_pitch,L_yaw] = inverse_Eular(L_R)

    if np.remainder(xtn,2)==0:       ##右脚支撑:
        f_angle_yaw[pick-1] = R_yaw  
    else:
        f_angle_yaw[pick-1] = L_yaw
        
    eps=0.0001
    ##地面角度估计：
    ##(左、右脚支撑时大部分时间时相同的)。
    if pick==1:        
        f_angle_roll[pick-1]=0
        f_angle_pitch[pick-1]=0
        f_angle_yaw[pick-1]=0
    else:
        if (Lf<=eps)and(Rf<=eps):     ##双脚离地    #@@@@@@eps不知道怎么表达
            f_angle_roll[pick-1] = f_angle_roll[pick-2]
            f_angle_pitch[pick-1] = f_angle_pitch[pick-2]
            f_angle_yaw[pick-1] = f_angle_yaw[pick-2]
        else:
            if Lf<=eps:          ##只有左脚悬空               
                if (((mea[13]>0)and(mea[14]>0)and(mea[15]>0)and(mea[16]>0))
                    or ((mea[13]<=0)and(mea[14]>0)and(mea[15]>0)and(mea[16]>0))
                    or ((mea[14]<=0)and(mea[13]>0)and(mea[15]>0)and(mea[16]>0))
                    or ((mea[15]<=0)and(mea[14]>0)and(mea[13]>0)and(mea[16]>0))
                    or ((mea[16]<=0)and(mea[14]>0)and(mea[13]>0)and(mea[15]>0))):
                    f_angle_roll[pick-1] = R_roll
                    f_angle_pitch[pick-1] = R_pitch                  
                else:
                    if (((mea[13]<=0)and(mea[14]<=0)and(mea[15]<=0)and(mea[16]>0))
                        or ((mea[13]<=0)and(mea[14]<=0)and(mea[16]<=0)and(mea[15]>0))
                        or ((mea[13]<=0)and(mea[15]<=0)and(mea[16]<=0)and(mea[14]>0))
                        or ((mea[14]<=0)and(mea[15]<=0)and(mea[16]<=0)and(mea[13]>0))):                        
                        f_angle_roll[pick-1] = f_angle_roll[pick-2]
                        f_angle_pitch[pick-1] = f_angle_pitch[pick-2]                     
                    else:
                        if ((mea[13]>0)and(mea[14]>0)) or ((mea[15]>0)and(mea[16]>0)):      ##roll角度与足部一致。
                            f_angle_roll[pick-1] = R_roll               
                            f_angle_pitch[pick-1] = f_angle_pitch[pick-2]                          
                        else:
                            if (((mea[13]>0)and(mea[15]>0))or((mea[14]>0)and(mea[16]>0))):   ##pitch角度
                                f_angle_roll[pick-1] = f_angle_roll[pick-2]               
                                f_angle_pitch[pick-1] = R_pitch                                
                            else:
                                f_angle_roll[pick-1] = R_roll
                                f_angle_pitch[pick-1] = R_pitch                              
                            
            else:
                if Rf<=eps:                                                             ##只有右脚悬空
                    if ((((mea[17]>0)and(mea[18]>0)and(mea[19]>0)and(mea[20]>0)))
                        or ((mea[17]<=0)and(mea[18]>0)and(mea[19]>0)and(mea[20]>0))
                        or ((mea[18]<=0)and(mea[17]>0)and(mea[19]>0)and(mea[20]>0))
                        or ((mea[19]<=0)and(mea[18]>0)and(mea[17]>0)and(mea[20]>0))
                        or ((mea[20]<=0)and(mea[18]>0)and(mea[17]>0)and(mea[19]>0))):                        
                        f_angle_roll[pick-1] = L_roll
                        f_angle_pitch[pick-1] = L_pitch
                    else:
                        if (((mea[17]<=0)and(mea[18]<=0)and(mea[19]<=0)and((mea[20]>0)))
                            or ((mea[17]<=0)and(mea[18]<=0)and(mea[20]<=0)and((mea[19]>0)))
                            or ((mea[17]<=0)and(mea[19]<=0)and(mea[20]<=0)and((mea[18]>0)))
                            or ((mea[18]<=0)and(mea[19]<=0)and(mea[20]<=0)and((mea[17]>0)))):                            
                            f_angle_roll[pick-1] = f_angle_roll[pick-2]
                            f_angle_pitch[pick-1] = f_angle_pitch[pick-2]                        
                        else:
                            if ((mea[17]>0)and(mea[18]>0))or((mea[19]>0)and(mea[20]>0)):      ##roll角度与足部一致。
                                f_angle_roll[pick-1] = L_roll               
                                f_angle_pitch[pick-1] = f_angle_pitch[pick-2]                           
                            else:
                                if ((mea[17]>0)and(mea[19]>0))or((mea[18]>0)and(mea[20]>0)):  ##pitch角度
                                    f_angle_roll[pick-1] = f_angle_roll[pick-2]               
                                    f_angle_pitch[pick-1] = L_pitch                                
                                else:
                                    f_angle_roll[pick-1] = L_roll
                                    f_angle_pitch[pick-1] = L_pitch                  
                               
                else:                                                           ## Lf>eps且Rf>eps
                    if ((((mea[17]<=0)and(mea[18]<=0)and(mea[19]<=0)and((mea[20]>0)))
                         or((mea[17]<=0)and(mea[18]<=0)and(mea[20]<=0)and((mea[19]>0)))
                         or((mea[17]<=0)and(mea[19]<=0)and(mea[20]<=0)and((mea[18]>0)))
                         or((mea[18]<=0)and(mea[19]<=0)and(mea[20]<=0)and((mea[17]>0))))
                        and(((mea[13]<=0)and(mea[14]<=0)and(mea[15]<=0)and((mea[16]>0)))
                            or((mea[13]<=0)and(mea[14]<=0)and(mea[16]<=0)and((mea[15]>0)))
                            or((mea[13]<=0)and(mea[15]<=0)and(mea[16]<=0)and((mea[14]>0)))
                            or(mea[14]<=0)and(mea[15]<=0)and(mea[16]<=0)and((mea[13]>0)))):

                        r_lf = np.reshape(R_lf,(1,4))
                        l_lf = np.reshape(L_lf,(1,4))
                        rx = find0(r_lf,0,1,4)
                        lx = find0(l_lf,0,1,4)                                       ##@@@@@find
                        rxx = Rverti[rx,:]
##                        rxx=list(rxx)
                        lxx = Lverti[lx,:]
##                        lxx=list(lxx)
                        if np.remainder(xtn,2)==0:                                      ##右脚支撑
                            det_r = lxx-rxx
                            f_angle_roll[pick-1] = math.atan(det_r[0,2]/det_r[0,1])
                            f_angle_pitch[pick-1] = -math.atan(det_r[0,2]/det_r[0,0])
                        else:
                            det_r = rxx-lxx
                            f_angle_roll[pick-1] = math.atan(det_r[0,2]/det_r[0,1])
                            f_angle_pitch[pick-1] = -math.atan(det_r[0,2]/det_r[0,0])                          
                                     
                    else:                                                       ##每个脚都有至少两个力传感器接触到地面:假设左右两脚的地面倾斜度一致。
                        if (((mea[13]>0)and(mea[14]>0)and(mea[15]<=0)and(mea[16]<=0))
                            or((mea[15]>0)and(mea[16]>0)and(mea[13]<=0)and(mea[14]<=0))):
                            
                            if (((mea[17]>0)and(mea[18]>0)and(mea[19]<=0)and(mea[20]<=0))
                                or((mea[19]>0)and(mea[20]>0)and(mea[17]<=0)and(mea[18]<=0))):                        
                                f_angle_pitch[pick-1] = f_angle_pitch[pick-2]
                                f_angle_roll[pick-1] = (L_roll+R_roll)/2                                         
                            else:
                                if (((mea[17]>0)and(mea[19]>0)and(mea[18]<=0)and(mea[20]<=0))
                                or((mea[18]>0)and(mea[20]>0)and(mea[17]<=0)and(mea[19]<=0))):
                                    f_angle_roll[pick-1] = R_roll
                                    f_angle_pitch[pick-1] = L_pitch
                                else:
                                    f_angle_roll[pick-1] = (R_roll+L_roll)/2
                                    f_angle_pitch[pick-1] = L_pitch                             
                              
                        else:
                            if (((mea[13]>0)and(mea[15]>0)and(mea[14]<=0)and(mea[16]<=0))
                            or((mea[14]*mea[16]>0)and(mea[13]<=0)and(mea[15]<=0))):
                                if (((mea[17]>0)and(mea[18]>0)and(mea[19]<=0)and(mea[20]<=0))
                                    or((mea[19]>0)and(mea[20]>0)and(mea[17]<=0)and(mea[18]<=0))):
                                    f_angle_roll[pick-1]=L_roll
                                    f_angle_pitch[pick-1]=R_pitch                                  
                                else:
                                    if (((mea[17]>0)and(mea[19]>0)and(mea[18]<=0)and(mea[20]<=0))
                                        or((mea[18]>0)and(mea[20]>0)and(mea[17]<=0)and(mea[19]<=0))):
                                        f_angle_roll[pick-1]=f_angle_roll[pick-2]
                                        f_angle_pitch[pick-1]=(L_pitch+R_pitch)/2                                                        
                                    else:
                                        f_angle_roll[pick-1]=L_roll
                                        f_angle_pitch[pick-1]=(L_pitch+R_pitch)/2                                   
                                  
                            else:
                                if (((mea[17]>0)and(mea[18]>0)and(mea[19]<=0)and(mea[20]<=0))
                                    or((mea[19]>0)and(mea[20]>0)and(mea[17]<=0)and(mea[18]<=0))):
                                        f_angle_roll[pick-1]=(L_roll+R_roll)/2
                                        f_angle_pitch[pick-1]=R_pitch                                    
                                else:
                                    if (((mea[17]>0)and(mea[19]>0)and(mea[18]<=0)and(mea[20]<=0))
                                        or((mea[18]>0)and(mea[20]>0)and(mea[17]<=0)and(mea[19]<=0))):
                                        
                                        f_angle_roll[pick-1]=R_roll
                                        f_angle_pitch[pick-1]=(L_pitch+R_pitch)/2                                                                      
                                    else:
                                        f_angle_roll[pick-1]=(L_roll+R_roll)/2
                                        f_angle_pitch[pick-1]=(L_pitch+R_pitch)/2
                                        
                                    
    slope_estimation[0] = f_angle_roll[pick-1]
    slope_estimation[1] = f_angle_pitch[pick-1]
    slope_estimation[2] = f_angle_yaw[pick-1]
    return slope_estimation



## pick:采样时间
## xtn:当前时刻对应的周期
## mea:反馈数据
## real_time:质心坐标
## Lf和Rf:左右脚接触力之和
def cal_foot_zmp_nao(pick,xtn,mea,real_com,Lf,Rf,mylink1):
    global Link
    global wf
    global zmp_g
    global com_g
    global rfoot_g
    global lfoot_g
    global n_period_N
    global R_sensor
    global L_sensor
    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones 
    Link = mylink1
    
    ##取出左脚运动链末端和右脚运动链末端的左边和姿态矩阵，并且计算出四个力传感器的坐标。
    l6 = Link[5].p + local_dot(Link[5].R, local_array([[0],[0],[-45.19]]))                           ##右足底坐标
    l11 = Link[10].p + local_dot(Link[10].R, local_array([[0],[0],[-45.19]]))                        ##左足底
    l6 = np.reshape(l6,(1,3))
    l11 = np.reshape(l11,(1,3))
    ##----------右脚
    Rtemp_verti = np.kron(Link[5].p,local_ones((1,4)))+ local_dot(Link[5].R,R_sensor.T)
    Rverti = Rtemp_verti.T
    ##----------左脚
    Ltemp_verti = np.kron(Link[10].p,local_ones((1,4))) + local_dot(Link[10].R,L_sensor.T)
    Lverti = Ltemp_verti.T
    ##左右四个力传感器的坐标
    L1 = Lverti[0,:]
    L2 = Lverti[1,:] 
    L3 = Lverti[2,:]
    L4 = Lverti[3,:]
    R1 = Rverti[0,:]
    R2 = Rverti[1,:]
    R3 = Rverti[2,:]
    R4 = Rverti[3,:]

    R_lf = local_array([mea[13],mea[14],mea[15],mea[16]])
    L_lf = local_array([mea[17],mea[18],mea[19],mea[20]])
    Rxfz = R1[0]*mea[13]+R2[0]*mea[14]+R3[0]*mea[15]+R4[0]*mea[16]
    Ryfz = R1[1]*mea[13]+R2[1]*mea[14]+R3[1]*mea[15]+R4[1]*mea[16]
    Lxfz = L1[0]*mea[17]+L2[0]*mea[18]+L3[0]*mea[19]+L4[0]*mea[20]
    Lyfz = L1[1]*mea[17]+L2[1]*mea[18]+L3[1]*mea[19]+L4[1]*mea[20]

    eps=0.0001
    if Rf<eps:
        Rzmpx=0
        Rzmpy=0
    else:
        Rzmpx=Rxfz/Rf
        Rzmpy=Ryfz/Rf
        
    if Lf<eps:
        Lzmpx=0
        Lzmpy=0
    else:
        Lzmpx=Lxfz/Lf
        Lzmpy=Lyfz/Lf
    
    if (Lf<eps)and(Rf<eps):
        ZMPx=0
        ZMPy=0
    else:
        ZMPx=(Rzmpx*Rf+Lzmpx*Lf)/(Rf+Lf)
        ZMPy=(Rzmpy*Rf+Lzmpy*Lf)/(Rf+Lf)
   

    ##求解全局坐标的算法其实是最关键的，这里不考虑支撑脚与地面的相对滑动。
    ##这里起步阶段xtn为1，第二步xtn为2;;;;但是对于python对应的下标则要比matlab少1，这里要加上1.

    xtn = xtn +1    
    if np.remainder(xtn,2)==0:         ##右脚支撑:
        ##相对于支撑足部踝关节坐标系
        zmpx_g = ZMPx-l6[0,0]
        zmpy_g = ZMPy-l6[0,1]
        zmp_g[pick-1,:] = local_array([zmpx_g,zmpy_g,0])
        com_g[pick-1,:] = local_array([real_com[0,0]-l6[0,0],real_com[1,0]-l6[0,1],real_com[2,0]-l6[0,2]])
        lfoot_g[pick-1,:] = l11-l6
        rfoot_g[pick-1,:] = local_array([0,0,0])
    else:
        
        zmpx_g = ZMPx-l11[0,0]
        zmpy_g = ZMPy-l11[0,1]
        zmp_g[pick-1,:] = local_array([zmpx_g,zmpy_g,0])
        com_g[pick-1,:] = local_array([real_com[0,0]-l11[0,0],real_com[1,0]-l11[0,1],real_com[2,0]-l11[0,2]])
        rfoot_g[pick-1,:] = l6-l11
        lfoot_g[pick-1,:] = local_array([0,0,0])
    

    ##每个支撑腿中心的全局坐标。
    if np.remainder(xtn,2)==0:      ##右脚支撑:
        if xtn>1:
            nx = n_period_N[0,xtn-2]
            if (pick-nx)==1:
                wf[pick-1,:] = rfoot_g[nx-1,:]+wf[pick-2,:]
            else:
                wf[pick-1,:] = wf[pick-2,:]
           
    else:
        if xtn>1:
            nx = n_period_N[0,xtn-2]
            if (pick-nx)==1:
                wf[pick-1,:] = lfoot_g[nx-1,:]+wf[pick-2,:]
            else:
                wf[pick-1,:] = wf[pick-2,:]

    

    slope_estimation = ground_inclination_estimation(pick,xtn,mea,R_lf,L_lf,Rverti,Lverti,Lf,Rf)

    wff = wf[pick-1,:]

    if (Lf<eps) and (Rf<eps):
        if pick>1:
            zmpf = zmp_g[pick-2,:]
        else:
            zmpf = zmp_g[pick-1,:]
    else:
        zmpf = zmp_g[pick-1,:]

    comf = com_g[pick-1,:]
    rfootf = rfoot_g[pick-1,:]
    lfootf = lfoot_g[pick-1,:]
    L11 = L1-l11
    L21 = L2-l11
    L31 = L3-l11
    L41 = L4-l11
    R11 = R1-l6
    R21 = R2-l6
    R31 = R3-l6
    R41 = R4-l6

    return wff,zmpf,comf,rfootf,lfootf,slope_estimation,L11,L21,L31,L41,R11,R21,R31,R41
            


####==========================================================================================
####========================================================================================
###  样条插值微分
def Diff(x,y):
    spl = itp.splrep(x, y,k=3)
    dspl = itp.splder(spl)
    dy = itp.splev(x, dspl)
    return dy


### ZMP到COM
def ZMP2COM(rt,dt,Zc,g,px,py,pz):
    #长度调整算法。
    tsh = 120
    tn =len(rt)-tsh
    rt = rt[tsh/2:(-tsh/2)]
    px1 = local_array(px)
    py1 = local_array(py)
    pz1 = local_array(pz)
    px2 = px1[tsh/2:(-tsh/2),0]    
    py2 = py1[tsh/2:(-tsh/2),0]
    pz2 = pz1[tsh/2:(-tsh/2),0]
    
    a = -Zc/(g*(dt**2))
    b = 2*Zc/(g*(dt**2))+1
    c = -Zc/(g*(dt**2))
    diag_1 = local_ones((tn))*b
    diag_2 = local_ones((tn-1))*c
    diag_3 = local_ones((tn-1))*a
    diag_1[0] = a+b
    diag_1[-1]= b+c
    A = np.diag(diag_2,1)+np.diag(diag_1)+np.diag(diag_3,-1)
    invsA = np.linalg.inv(A)   

    px2.reshape(tn)
    py2.reshape(tn)
    pz2.reshape(tn)
    
    vx = Diff(rt,px2)
    px2[0] = px2[0] + a*vx[0]*dt
    px2[-1] = px2[-1]-a*vx[-1]*dt
    px2 = px2.reshape(tn,1)
    com_x = local_dot(invsA,px2)

    vy = Diff(rt,py2)
    py2[0] = py2[0] + a*vy[0]*dt
    py2[-1] = py2[-1]-a*vy[-1]*dt
    py2 = py2.reshape(tn,1)    
    com_y = local_dot(invsA,py2)
##  
    vz = Diff(rt,pz2)
    pz2[0] = pz2[0] + a*vz[0]*dt
    pz2[-1] = pz2[-1]-a*vz[-1]*dt
    pz2 = pz2.reshape(tn,1)    
    com_z = local_dot(invsA,pz2)

    return com_x,com_y,com_z,tsh       ###n*3矩阵，R为3*3矩阵


## -----------------------------walking modes transition-------------------------------------
##---walking on the unslope
def optimal_upslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN):
    global x_init
    global x_end
    global Vx_init
    global Vx_end 
    global y_init
    global y_end
    global Vy_init
    global Vy_end
    global z_init
    global z_end
    global Vz_init
    global Vz_end
##  PX1组存储实时ZMP轨迹  
    global PX1
    global PY1
    global PZ1     
    
    global footx
    global footy
    global footz
    global DSP
    global SSP
    global Sx
    global Sy
    global Sz    
    global z_h
    global slope_angle
    global flag
    global NNx

##  为普通的全局变量。    
    px =local_zeros([XN,1])
    py =local_zeros([XN,1])
    pz =local_zeros([XN,1])

    
    flag[n_stnum-1,0] = -1

    n_global = n_stnum
##  t_stnum为0时，会出现下标溢出
    if (t_stnum - dt>0.01):
        t_yud = T-t_stnum
        n_t_yud = int(local_floor(t_yud/dt)+1)
        px[0:(nx+n_t_yud),0] = PX1[1:(nx+n_t_yud+1),0]
        py[0:(nx+n_t_yud),0] = PY1[1:(nx+n_t_yud+1),0]
        pz[0:(nx+n_t_yud),0] = PZ1[1:(nx+n_t_yud+1),0]
                    
    ##%步长调整：判断应该是左脚支撑还是右脚支撑
##      左脚支撑
        if np.remainder(n_stnum,2)==0 :                       
           Sx[n_int-1,0] = outx[9]                   
           Sy[n_int-1,0] = abs(outx[10])                                                 
           z_h[n_int-1,0] = abs(outx[11])                                              
        else:                                                          
           Sx[n_int-1,0] = outx[12]
           Sy[n_int-1,0] = abs(outx[13])             
           z_h[n_int-1,0] = abs(outx[14])

        if (isinstance(Sx[n_int-1,0], float) ):
            if math.isnan(Sx[n_int-1,0]):
                Sx[n_int-1,0] = (Sx[n_int-2,0]+Sx[n_int-1,0])/2
                if  Sx[n_int-1,0]>=80:
                   Sx[n_int-1,0] = 80
                else:
                    if Sx[n_int-1,0]<=30:
                        Sx[n_int-1,0] = 30
            else:
                Sx[n_int-1,0]=30
        else:
            Sx[n_int-1,0]=30
                
        if (isinstance(Sy[n_int-1,0], float) ):
            if math.isnan(Sy[n_int-1,0]):
                Sy[n_int-1,0] = (Sy[n_int-2,0]+Sy[n_int-1,0])/2
                if Sy[n_int-1,0]>=120:
                   Sy[n_int-1,0]=120
                else:
                    if Sy[n_int-1,0]<=100:
                        Sy[n_int-1,0]=100
            else:
                Sy[n_int-1,0]=100
        else:
            Sy[n_int-1,0]=100
                
##      初始状态判断;
        if n_int<=8:
            slope_angle[n_global,0] = -abs(outx[16])/5   
        else:
            if n_int<=19:
                slope_angle[n_global,0] = -abs(outx[16])
            else:
                slope_angle[n_global,0] = -2*abs(outx[16])
##                slope_angle[n_global,0] = -abs(outx[16])
                
        z_h[n_int-1,0] = Sx[n_int-1,0]* math.tan(abs(slope_angle[n_global,0]))
        
        for i in range(n_global+1,NNx):
            slope_angle[i,0] = slope_angle[n_global,0]

        for i in range(n_int,NNx):
             Sx[i,0] = Sx[n_int-1,0]
             Sy[i,0] = Sy[n_int-1,0]                    
             z_h[i,0] = z_h[n_int-1,0]
             Sz[i,0] = Sz[i-1,0]+z_h[i-1,0]
             footz[i,0] = footz[i-1,0]+z_h[i-1,0]             
                               

    ##%%%%%由当前时刻所在步行周期的位置计算现在意义以后的质心。         
    ##%%%%%%计算下一周期的步长调整量。       
    ##%调整下一周期的双组相 
        SSP[n_int,0] = SSP[n_int - 1,0]-5*dt
        DSP[n_int,0] = T-SSP[n_int,0]

        
        if DSP[n_int,0]<=0.05*T:
            DSP[n_int,0]= 0.05*T
        else:
            if DSP[n_int,0]>=0.5*T:
                DSP[n_int,0]= 0.5*T
                
        SSP[n_int,0] = T - DSP[n_int,0]

##  判断是否出现0.99999之类的情况,防止双组相插值计算过程中出现数值跳变
        n_ssp = int(1000*SSP[n_int - 1,0])
        n_dsp = int(1000*DSP[n_int - 1,0])
        if ((n_ssp+1)- 1000*SSP[n_int - 1,0])<0.2:
            SSP[n_int - 1,0] = (n_ssp+1)/ 1000
            
        if ((n_dsp+1)- 1000*DSP[n_int - 1,0])<0.2:
            DSP[n_int - 1,0] = (n_dsp+1)/ 1000

##      接下来所有周期
        for i in range(n_int+1,NNx):
           SSP[i,0] = SSP[n_int,0]
           DSP[i,0] = DSP[n_int,0] 

                         
    ##    %=======单组相,最小能量控制，求解析解==
        Kssp = int(local_floor(SSP[n_int-1,0]/dt))
        px0 = local_zeros([2,1])
        py0 = local_zeros([2,1])
        pz0 = local_zeros([2,1])
        
    ##%双组相的起始位置
        for i in range(Kssp-1,Kssp+1):
            px0[i+1-Kssp,0] = footx[n_int-1,0]
            py0[i+1-Kssp,0] = footy[n_int-1,0]
            pz0[i+1-Kssp,0] = Sz[n_int-1,0]

            

    ##    %下一周期的位置和速度。
        n_int = n_intx
        footx[n_int-1,0] = footx[n_int-2,0]+Sx[n_int-2,0]
        footy[n_int-1,0] = footy[n_int-2,0]+(-1)**(n_int-1)*Sy[n_int-2,0]
        footx[n_int,0] = footx[n_int-1,0]+Sx[n_int-1,0]
        footy[n_int,0] = footy[n_int-1,0]+(-1)**(n_int)*Sy[n_int-1,0]
        footz[n_int-1,0] = footz[n_int-2,0]+z_h[n_int-1,0]
        Sz[n_int-1,0] = Sz[n_int-2,0]+z_h[n_int-2,0]
        footz[n_int,0] = footz[n_int-1,0]+z_h[n_int,0]
        Sz[n_int,0] = Sz[n_int-1,0]+z_h[n_int-1,0] 
       

##  判断是否出现0.99999之类的情况,防止双组相插值计算过程中出现数值跳变
        n_ssp = int(1000*SSP[n_int - 1,0])
        n_dsp = int(1000*DSP[n_int - 1,0])
        if ((n_ssp+1)- 1000*SSP[n_int - 1,0])<0.2:
            SSP[n_int - 1,0] = (n_ssp+1)/ 1000
            
        if ((n_dsp+1)- 1000*DSP[n_int - 1,0])<0.2:
            DSP[n_int - 1,0] = (n_dsp+1)/ 1000
            
        if DSP[n_int,0]<=0.05*T:
            DSP[n_int,0]= 0.05*T
        else:
            if DSP[n_int,0]>=0.5*T:
                DSP[n_int,0]= 0.5*T
                
        SSP[n_int,0] = T - DSP[n_int,0]            
             
    ##    %=======单组相,最小能量控制，求解析解==
        pxx = local_zeros([2,1])
        pyy = local_zeros([2,1])
        pzz = local_zeros([2,1])
        
    ##%双组相的起始位置
        for i in range(1,3):
            pxx[i-1,0] = footx[n_int-1,0]
            pyy[i-1,0] = footy[n_int-1,0]
            pzz[i-1,0] = Sz[n_int-1,0]             

    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式。。
    ##%!!!!!!!!!!!!!!!!!!!注意，这时候要调整双足相期间质心的高度位移
    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式。。

        x_init[0,0]=px0[1,0]
        Vx_init[0,0]=(px0[1,0]-px0[0,0])/dt
        apx_init=0
        x_end[0,0]=pxx[0,0]
        Vx_end[0,0]=(pxx[1,0]-pxx[0,0])/dt
        apx_end=0
        y_init[0,0]=py0[1,0]
        Vy_init[0,0]=(py0[1,0]-py0[0,0])/dt
        apy_init=0
        y_end[0,0]=pyy[0,0]
        Vy_end[0,0]=(pyy[1,0]-pyy[0,0])/dt
        apy_end=0
        z_init[0,0]=pz0[1,0]
        Vz_init[0,0]=(pz0[1,0]-pz0[0,0])/dt
        apz_init=0
        z_end[0,0]=pzz[0,0]
        Vz_end[0,0]=(pzz[1,0]-pzz[0,0])/dt
        apz_end=0
       
        ds3= (DSP[n_int-1,0])**3
        ds4= (DSP[n_int-1,0])**4
        ds5= (DSP[n_int-1,0])**5
        det_px = px0[1,0]-pxx[0,0]
        det_py = py0[1,0]-pyy[0,0]
        det_pz = pz0[1,0]-pzz[0,0]
        
        ax = (np.array([px0[1,0],0,0,-10*det_px/ds3,15*det_px/ds4,-6*det_px/ds5])).reshape([1,6])
        ay = (np.array([py0[1,0],0,0,-10*det_py/ds3,15*det_py/ds4,-6*det_py/ds5])).reshape([1,6])
        az = (np.array([pz0[1,0],0,0,-10*det_pz/ds3,15*det_pz/ds4,-6*det_pz/ds5])).reshape([1,6])
        
        K = int(local_floor(DSP[n_int-1,0]/dt))
##        if ((int(100*(K+1)*dt) - 100*DSP[n_int-1,0])<= 0.2):
##            K = K + 1
            
        t_end = t[-1]-(n_int-1)*T
        if t_end > DSP[n_int-1,0]-0.0001:
            for i in range(1,K+1):
                t_yu = i*dt
                px[nx+n_t_yud+i-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5  
                py[nx+n_t_yud+i-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
                pz[nx+n_t_yud+i-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5
     
            t_yu1 = t_end-DSP[n_int-1,0] 
            k1 = int(local_floor((t_yu1)/dt))
            if ((int(100*(k1+1)*dt) - 100*t_yu1)<= 0.2):
                k1 = k1 + 1            
            for i in range(1,k1+1):
                px[nx+n_t_yud+K+i-1,0] = footx[n_int-1,0]
                py[nx+n_t_yud+K+i-1,0] = footy[n_int-1,0] 
                pz[nx+n_t_yud+K+i-1,0] = Sz[n_int-1,0]
        else:
            if abs(t_end - DSP[n_int-1,0])<0.0001:
                for i in range(1,K+1):
                    t_yu = i*dt
                    px[nx+n_t_yud+i-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5  
                    py[nx+n_t_yud+i-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
                    pz[nx+n_t_yud+i-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5
            else:
                K1 = int(local_floor(t_end/dt))
                if ((int(100*(K1+1)*dt) - 100*t_end)<= 0.2):
                    K1 = K1 + 1   
                for i in range(1,K1+1):
                    t_yu = i*dt
                    px[nx+n_t_yud+i-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5  
                    py[nx+n_t_yud+i-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
                    pz[nx+n_t_yud+i-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5
    else:
        px[0:XN-1,0] = PX1[1:XN,0]
        py[0:XN-1,0] = PY1[1:XN,0]
        pz[0:XN-1,0] = PZ1[1:XN,0]
        px[-1,0] = PX1[-1,0]
        py[-1,0] = PY1[-1,0]
        pz[-1,0] = PZ1[-1,0]                      
           
    PX1 = px
    PY1 = py
    PZ1 = pz

    return px,py,pz

##---walking on the downslope
def optimal_downslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN):
    global x_init
    global x_end
    global Vx_init
    global Vx_end 
    global y_init
    global y_end
    global Vy_init
    global Vy_end
    global z_init
    global z_end
    global Vz_init
    global Vz_end    
 
    global footx
    global footy
    global footz
    global DSP
    global SSP
    global Sx
    global Sy
    global Sz    
    global z_h
    global slope_angle
    global flag    
    global NNx

##  px0组存储实时ZMP轨迹  
    global PX1
    global PY1
    global PZ1


##  为普通的全局变量。    
    px =local_zeros([XN,1])
    py =local_zeros([XN,1])
    pz =local_zeros([XN,1])
    
    flag[n_stnum-1,:] = 1

    n_global = n_stnum
    
    if (t_stnum - dt>0.01):
        
        t_yud = T-t_stnum        
        n_t_yud = int(local_floor(t_yud/dt)+1)
        px[0:(nx+n_t_yud),0] = PX1[1:(nx+n_t_yud+1),0]
        py[0:(nx+n_t_yud),0] = PY1[1:(nx+n_t_yud+1),0]
        pz[0:(nx+n_t_yud),0] = PZ1[1:(nx+n_t_yud+1),0]
                    
    ##%步长调整：判断应该是左脚支撑还是右脚支撑
        if np.remainder(n_stnum,2)==0 :                       
           Sx[n_int-1,0] = outx[9]                   
           Sy[n_int-1,0] = abs(outx[10])                                                 
           z_h[n_int-1,0] = -abs(outx[11])/2                                              
        else:                                                          
           Sx[n_int-1,0] = outx[12]
           Sy[n_int-1,0] = abs(outx[13])             
           z_h[n_int-1,0] = -abs(outx[14])/2


        if (isinstance(Sx[n_int-1,0], float) ):
            if math.isnan(Sx[n_int-1,0]):
                Sx[n_int-1,0] = (Sx[n_int-2,0]+Sx[n_int-1,0])/2
                if  Sx[n_int-1,0]>=80:
                   Sx[n_int-1,0] = 80
                else:
                    if Sx[n_int-1,0]<=30:
                        Sx[n_int-1,0] = 30
            else:
                Sx[n_int-1,0]=50
        else:
            Sx[n_int-1,0]=50
                
        if (isinstance(Sy[n_int-1,0], float) ):
            if math.isnan(Sy[n_int-1,0]):
                Sy[n_int-1,0] = (Sy[n_int-2,0]+Sy[n_int-1,0])/2
                if Sy[n_int-1,0]>=120:
                   Sy[n_int-1,0]=120
                else:
                    if Sy[n_int-1,0]<=100:
                        Sy[n_int-1,0]=100
            else:
                Sy[n_int-1,0]=110
        else:
            Sy[n_int-1,0]=110

        slope_angle[n_global,0] = abs(outx[16])/2
        for i in range(n_global+1,NNx):
            slope_angle[i,0] = slope_angle[n_global,0]

        z_h[n_int-1,0] = -Sx[n_int-1,0]*math.tan(abs(slope_angle[n_global,0]))
        
        for i in range(n_int,NNx):
             Sx[i,0] = Sx[n_int-1,0]
             Sy[i,0] = Sy[n_int-1,0]                    
             z_h[i,0] = z_h[n_int-1,0]
             Sz[i,0] = Sz[i-1,0]+z_h[i-1,0]
             footz[i,0] = footz[n_int-1,0]+z_h[i-1,0]
##  判断是否出现0.99999之类的情况,防止双组相插值计算过程中出现数值跳变
        n_ssp = int(1000*SSP[n_int - 1,0])
        n_dsp = int(1000*DSP[n_int - 1,0])
        if ((n_ssp+1)- 1000*SSP[n_int - 1,0])<0.2:
            SSP[n_int - 1,0] = (n_ssp+1)/ 1000
            
        if ((n_dsp+1)- 1000*DSP[n_int - 1,0])<0.2:
            DSP[n_int - 1,0] = (n_dsp+1)/ 1000
             
    ##    %=======单组相,最小能量控制，求解析解==
        Kssp = int(local_floor(SSP[n_int-1,0]/dt))
        px0 = local_zeros([2,1])
        py0 = local_zeros([2,1])
        pz0 = local_zeros([2,1])
        
    ##%双组相的起始位置
        for i in range(Kssp-1,Kssp+1):
            px0[i+1-Kssp,0] = footx[n_int-1,0]
            py0[i+1-Kssp,0] = footy[n_int-1,0]
            pz0[i+1-Kssp,0] = Sz[n_int-1,0]            

    ##    %下一周期的位置和速度。
        n_int = n_intx
        footx[n_int-1,0] = footx[n_int-2,0]+Sx[n_int-2,0]
        footy[n_int-1,0] = footy[n_int-2,0]+(-1)**(n_int-1)*Sy[n_int-2,0]
        footx[n_int,0] = footx[n_int-1,0]+Sx[n_int-1,0]
        footy[n_int,0] = footy[n_int-1,0]+(-1)**(n_int)*Sy[n_int-1,0]
        
        footz[n_int-1,0] = footz[n_int-2,0]+z_h[n_int-1,0]
        Sz[n_int-1,0] = Sz[n_int-2,0]+z_h[n_int-2,0]
        footz[n_int,0] = footz[n_int-1,0]+z_h[n_int,0]
        Sz[n_int,0] = Sz[n_int-1,0]+z_h[n_int-1,0] 
        
##  判断是否出现0.99999之类的情况,防止双组相插值计算过程中出现数值跳变
        n_ssp = int(1000*SSP[n_int - 1,0])
        n_dsp = int(1000*DSP[n_int - 1,0])
        if ((n_ssp+1)- 1000*SSP[n_int - 1,0])<0.2:
            SSP[n_int - 1,0] = (n_ssp+1)/ 1000
            
        if ((n_dsp+1)- 1000*DSP[n_int - 1,0])<0.2:
            DSP[n_int - 1,0] = (n_dsp+1)/ 1000
            
    ##    %=======单组相,最小能量控制，求解析解==
        pxx = local_zeros([2,1])
        pyy = local_zeros([2,1])
        pzz = local_zeros([2,1])

        
    ##%双组相的起始位置
        for i in range(1,3):
            pxx[i-1,0] = footx[n_int-1,0]
            pyy[i-1,0] = footy[n_int-1,0]
            pzz[i-1,0] = Sz[n_int-1,0]             

    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式。。
    ##%!!!!!!!!!!!!!!!!!!!注意，这时候要调整双足相期间质心的高度位移
    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式。。

        x_init[0,0]=px0[1,0]
        Vx_init[0,0]=(px0[1,0]-px0[0,0])/dt
        apx_init=0
        x_end[0,0]=pxx[0,0]
        Vx_end[0,0]=(pxx[1,0]-pxx[0,0])/dt
        apx_end=0
        y_init[0,0]=py0[1,0]
        Vy_init[0,0]=(py0[1,0]-py0[0,0])/dt
        apy_init=0
        y_end[0,0]=pyy[0,0]
        Vy_end[0,0]=(pyy[1,0]-pyy[0,0])/dt
        apy_end=0
        z_init[0,0]=pz0[1,0]
        Vz_init[0,0]=(pz0[1,0]-pz0[0,0])/dt
        apz_init=0
        z_end[0,0]=pzz[0,0]
        Vz_end[0,0]=(pzz[1,0]-pzz[0,0])/dt
        apz_end=0
       
        ds3= (DSP[n_int-1,0])**3
        ds4= (DSP[n_int-1,0])**4
        ds5= (DSP[n_int-1,0])**5
        det_px = px0[1,0]-pxx[0,0]
        det_py = py0[1,0]-pyy[0,0]
        det_pz = pz0[1,0]-pzz[0,0]
        
        ax = (np.array([px0[1,0],0,0,-10*det_px/ds3,15*det_px/ds4,-6*det_px/ds5])).reshape([1,6])
        ay = (np.array([py0[1,0],0,0,-10*det_py/ds3,15*det_py/ds4,-6*det_py/ds5])).reshape([1,6])
        az = (np.array([pz0[1,0],0,0,-10*det_pz/ds3,15*det_pz/ds4,-6*det_pz/ds5])).reshape([1,6])
        
        K = int(local_floor(DSP[n_int-1,0]/dt))
##        if ((int(100*(K+1)*dt) - 100*DSP[n_int-1,0])<= 0.2):
##            K = K + 1
            
        t_end = t[-1]-(n_int-1)*T
        if t_end > DSP[n_int-1,0]-0.0001:
            for i in range(1,K+1):
                t_yu = i*dt
                px[nx+n_t_yud+i-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5  
                py[nx+n_t_yud+i-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
                pz[nx+n_t_yud+i-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5
     
            t_yu1 = t_end-DSP[n_int-1,0] 
            k1 = int(local_floor((t_yu1)/dt))
            if ((int(100*(k1+1)*dt) - 100*t_yu1)<= 0.2):
                k1 = k1 + 1            
            for i in range(1,k1+1):
                px[nx+n_t_yud+K+i-1,0] = footx[n_int-1,0]
                py[nx+n_t_yud+K+i-1,0] = footy[n_int-1,0] 
                pz[nx+n_t_yud+K+i-1,0] = Sz[n_int-1,0]
        else:
            if abs(t_end - DSP[n_int-1,0])<0.0001:
                for i in range(1,K+1):
                    t_yu = i*dt
                    px[nx+n_t_yud+i-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5  
                    py[nx+n_t_yud+i-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
                    pz[nx+n_t_yud+i-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5
            else:
                K1 = int(local_floor(t_end/dt))
                if ((int(100*(K1+1)*dt) - 100*t_end)<= 0.2):
                    K1 = K1 + 1   
                for i in range(1,K1+1):
                    t_yu = i*dt
                    px[nx+n_t_yud+i-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5  
                    py[nx+n_t_yud+i-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
                    pz[nx+n_t_yud+i-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5
    else:
        px[0:XN-1,0] = PX1[1:XN,0]
        py[0:XN-1,0] = PY1[1:XN,0]
        pz[0:XN-1,0] = PZ1[1:XN,0]
        px[-1,0] = PX1[-1,0]
        py[-1,0] = PY1[-1,0]
        pz[-1,0] = PZ1[-1,0]                      



    PX1 = px
    PY1 = py
    PZ1 = pz

    return px,py,pz      

##---walking on the flatground
def optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN):
    global x_init
    global x_end
    global Vx_init
    global Vx_end 
    global y_init
    global y_end
    global Vy_init
    global Vy_end
    global z_init
    global z_end
    global Vz_init
    global Vz_end     
    global footx
    global footy
    global footz
    global DSP
    global SSP
    global Sx
    global Sy
    global Sz    
    global z_h
    global slope_angle
    global flag    

##  px0组存储实时ZMP轨迹  
    global PX1
    global PY1
    global PZ1
    
##  为普通的全局变量。    
    px =local_zeros([XN,1])
    py =local_zeros([XN,1])
    pz =local_zeros([XN,1])

    flag[n_stnum-1,0] = 0

    if (t_stnum - dt>0.01):
    
        px[0:-1,0] = PX1[1:,0]
        py[0:-1,0] = PY1[1:,0]
        pz[0:-1,0] = PZ1[1:,0]

        n_int = n_intx
        n_int = n_int-1 
        footx[n_int-1,0] = footx[n_int-2,0]+Sx[n_int-2,0]
        footy[n_int-1,0] = footy[n_int-2,0]+(-1)**(n_int-1)*Sy[n_int-2,0]              


    ##  判断是否出现0.99999之类的情况,防止双组相插值计算过程中出现数值跳变
        n_ssp = int(1000*SSP[n_int - 1,0])
        n_dsp = int(1000*DSP[n_int - 1,0])
        if ((n_ssp+1)- 1000*SSP[n_int - 1,0])<0.2:
            SSP[n_int - 1,0] = (n_ssp+1)/ 1000
            
        if ((n_dsp+1)- 1000*DSP[n_int - 1,0])<0.2:
            DSP[n_int - 1,0] = (n_dsp+1)/ 1000
                         
    ##    %=======单组相,最小能量控制，求解析解==
        Kssp = int(local_floor(SSP[n_int-1,0]/dt))
        px0 = local_zeros([2,1])
        py0 = local_zeros([2,1])
        pz0 = local_zeros([2,1])
        
    ##%双组相的起始位置
        for i in range(Kssp-1,Kssp+1):        
            px0[i+1-Kssp,0] = footx[n_int-1,0]
            py0[i+1-Kssp,0] = footy[n_int-1,0]
            pz0[i+1-Kssp,0] = Sz[n_int-1,0]            

    ##    %下一周期的位置和速度。
        n_int = n_intx
        footx[n_int-1,0] = footx[n_int-2,0]+Sx[n_int-2,0]
        footy[n_int-1,0] = footy[n_int-2,0]+(-1)**(n_int-1)*Sy[n_int-2,0]       
        footx[n_int,0] = footx[n_int-1,0]+Sx[n_int-1,0]
        footy[n_int,0] = footy[n_int-1,0]+(-1)**(n_int)*Sy[n_int-1,0]
        
        footz[n_int-1,0] = footz[n_int-2,0]+z_h[n_int-1,0]
        Sz[n_int-1,0] = Sz[n_int-2,0]+z_h[n_int-2,0]
        footz[n_int,0] = footz[n_int-1,0]+z_h[n_int,0]
        Sz[n_int,0] = Sz[n_int-1,0]+z_h[n_int-1,0]        

    ##  判断是否出现0.99999之类的情况,防止双组相插值计算过程中出现数值跳变
        n_ssp = int(1000*SSP[n_int - 1,0])
        n_dsp = int(1000*DSP[n_int - 1,0])
        if ((n_ssp+1)- 1000*SSP[n_int - 1,0])<0.2:
            SSP[n_int - 1,0] = (n_ssp+1)/ 1000
                
        if ((n_dsp+1)- 1000*DSP[n_int - 1,0])<0.2:
            DSP[n_int - 1,0] = (n_dsp+1)/ 1000
                
                
    ##    %=======单组相,最小能量控制，求解析解==
        pxx = local_zeros([2,1])
        pyy = local_zeros([2,1])
        pzz = local_zeros([2,1])
        
    ##%双组相的起始位置
        for i in range(1,3):
            pxx[i-1,0] = footx[n_int-1,0]
            pyy[i-1,0] = footy[n_int-1,0]
            pzz[i-1,0] = Sz[n_int-1,0]                      


    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式。。
    ##%!!!!!!!!!!!!!!!!!!!注意，这时候要调整双足相期间质心的高度位移
    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式。。  
        ds3= (DSP[n_int-1,0])**3
        ds4= (DSP[n_int-1,0])**4
        ds5= (DSP[n_int-1,0])**5
        det_px = px0[1,0]-pxx[0,0]
        det_py = py0[1,0]-pyy[0,0]
        det_pz = pz0[1,0]-pzz[0,0]
        
        ax = (np.array([px0[1,0],0,0,-10*det_px/ds3,15*det_px/ds4,-6*det_px/ds5])).reshape([1,6])
        ay = (np.array([py0[1,0],0,0,-10*det_py/ds3,15*det_py/ds4,-6*det_py/ds5])).reshape([1,6])
        az = (np.array([pz0[1,0],0,0,-10*det_pz/ds3,15*det_pz/ds4,-6*det_pz/ds5])).reshape([1,6])
        
        t_yu = t[-1]-(n_int-1)*T  #预测
        if (t_yu <= DSP[n_int-1,0]+0.001) and (t_yu >= dt-0.0001):
            px[-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5    ##最后一个值？
            py[-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
            pz[-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5
        else:
            if (t_yu > DSP[n_int-1,0] + 0.01):
                px[-1,0] = footx[n_int-1,0]
                py[-1,0] = footy[n_int-1,0]
                pz[-1,0] = Sz[n_int-1,0]
    else:
        px[0:XN-1,0] = PX1[1:XN,0]
        py[0:XN-1,0] = PY1[1:XN,0]
        pz[0:XN-1,0] = PZ1[1:XN,0]
        px[-1,0] = PX1[-1,0]
        py[-1,0] = PY1[-1,0]
        pz[-1,0] = PZ1[-1,0] 

            
    PX1 = px
    PY1 = py
    PZ1 = pz
    
    return px,py,pz
        
        
##--------------------------------------------------------------------------------------------------------
##======================================================================================
##-------------------------- 在线规划com
## -------最小能量法求解ZMP与com轨迹
##便于处理,也就是12+dt开始处理
def position_ZMP_COM_energymin_footplace(Tc,stnum,dt,t,XN,n_int,Zc,g,outx,pick,T,nx):
##def position_ZMP_COM_energymin_footplace():    
    global x_init
    global x_end
    global Vx_init
    global Vx_end 
    global y_init
    global y_end
    global Vy_init
    global Vy_end
    global z_init
    global z_end
    global Vz_init
    global Vz_end
##  px0组存储实时ZMP轨迹  
    global PX1
    global PY1
    global PZ1
    
    global footx
    global footy
    global footz
    global DSP
    global SSP
    global Sx
    global Sy
    global Sz    
    global z_h
    global flagx
    global T_initial
    global t_land 
    global Mea
    global flag
    global zmp_ref
    global com_ref
    global zmp_offline

    global flagxn

    length = 10
    width = 10    

    if stnum == 8+dt:
        
    ##  为普通的全局变量。    
        px =local_zeros([XN,1])
        py =local_zeros([XN,1])
        pz =local_zeros([XN,1])
        
        flagx[pick-1,0] = 5
        flagx[pick-1,1] = 0
##        %第一计算，前半部分去离线的结果
        nxxx1 = int(local_floor(T/dt))
        pic = int(local_floor((stnum-T+dt+T_initial-T)/dt))
        px[0:nx,0] = zmp_offline[pic:pic+nx,0]
        py[0:nx,0] = zmp_offline[pic:pic+nx,1]
        pz[0:nx,0] = zmp_offline[pic:pic+nx,2]
        
        n_int=n_int-1
        footx[n_int-1,0] = footx[n_int-2,0]+Sx[n_int-2,0]
        footy[n_int-1,0] = footy[n_int-2,0]+(-1)**(n_int-1)*Sy[n_int-2,0]
        footz[n_int-1,0] = footz[n_int-2,0]+z_h[n_int-1,0]
        Sz[n_int-1,0] = Sz[n_int-2,0]+z_h[n_int-2,0]        
        
             
    ##    %=======单组相,最小能量控制，求解析解==
        Kssp = int(local_floor(SSP[n_int-1,0]*100/dt/100))
        px0 = local_zeros([2,1])
        py0 = local_zeros([2,1])
        pz0 = local_zeros([2,1])        
    ##%双组相的起始位置:充分利用单组的零位移特性。
        for i in range(Kssp-1,Kssp+1):           
            px0[i+1-Kssp,0] = footx[n_int-1,0]
            py0[i+1-Kssp,0] = footy[n_int-1,0]
            pz0[i+1-Kssp,0] = Sz[n_int-1,0]                                
 

##%当前周期的局部坐标系的位置------
        n_int = n_int + 1
        footx[n_int-1,0] = footx[n_int-2,0]+Sx[n_int-2,0]
        footy[n_int-1,0] = footy[n_int-2,0]+(-1)**(n_int-1)*Sy[n_int-2,0]
        footx[n_int,0] = footx[n_int-1,0]+Sx[n_int-1,0]
        footy[n_int,0] = footy[n_int-1,0]+(-1)**(n_int)*Sy[n_int-1,0]
        
        footz[n_int-1,0] = footz[n_int-2,0]+z_h[n_int-1,0]
        Sz[n_int-1,0] = Sz[n_int-2,0]+z_h[n_int-2,0]
        footz[n_int,0] = footz[n_int-1,0]+z_h[n_int,0]
        Sz[n_int,0] = Sz[n_int-1,0]+z_h[n_int-1,0]       
              
    ##    %=======单组相,最小能量控制，求解析解==
        K = int(local_floor(DSP[n_int-1,0]*100/dt/100))
        Kssp = int(local_floor(SSP[n_int-1,0]*100/dt/100))        
        pxxn = int(local_floor(T/dt))
        pxx = local_zeros([pxxn,1])
        pyy = local_zeros([pxxn,1])
        pzz = local_zeros([pxxn,1])

        
    ##%双组相的起始位置
        for i in range(1,Kssp+1):                       
            pxx[i-1,0] = footx[n_int-1,0]
            pyy[i-1,0] = footy[n_int-1,0]
            pzz[i-1,0] = Sz[n_int-1,0] 
            px[nx+K+i-1,0] = pxx[i-1,0]
            py[nx+K+i-1,0] = pyy[i-1,0]
            pz[nx+K+i-1,0] = pzz[i-1,0]                

        
    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式。。
    ##%!!!!!!!!!!!!!!!!!!!注意，这时候要调整双足相期间质心的高度位移
    ##%---双足相始、末ZMP位置和速度, 双组相ZMP的五次多项式
        ##  考虑到单组相期间zero_ZMP位移，边界条件中速度和加速度都为0
    
        ds3= (DSP[n_int-1,0])**3
        ds4= (DSP[n_int-1,0])**4
        ds5= (DSP[n_int-1,0])**5
        det_px = px0[1,0]-pxx[0,0]
        det_py = py0[1,0]-pyy[0,0]
        det_pz = pz0[1,0]-pzz[0,0]
        ax = (np.array([px0[1,0],0,0,-10*det_px/ds3,15*det_px/ds4,-6*det_px/ds5])).reshape([1,6])
        ay = (np.array([py0[1,0],0,0,-10*det_py/ds3,15*det_py/ds4,-6*det_py/ds5])).reshape([1,6])
        az = (np.array([pz0[1,0],0,0,-10*det_pz/ds3,15*det_pz/ds4,-6*det_pz/ds5])).reshape([1,6])        
        for i in range(1,K+1):
            t_yu = i*dt
            px[nx+i-1,0] = ax[:,0]+ax[:,1]*(t_yu)+ax[:,2]*(t_yu)**2+ax[:,3]*(t_yu)**3+ax[:,4]*(t_yu)**4+ax[:,5]*(t_yu)**5  
            py[nx+i-1,0] = ay[:,0]+ay[:,1]*(t_yu)+ay[:,2]*(t_yu)**2+ay[:,3]*(t_yu)**3+ay[:,4]*(t_yu)**4+ay[:,5]*(t_yu)**5
            pz[nx+i-1,0] = az[:,0]+az[:,1]*(t_yu)+az[:,2]*(t_yu)**2+az[:,3]*(t_yu)**3+az[:,4]*(t_yu)**4+az[:,5]*(t_yu)**5

        PX1 = px
        PY1 = py
        PZ1 = pz
    
    else:      
        n_intx = n_int
##        %模式判断：
##        %落地碰撞检测
        n_stnum = int(np.floor((stnum-dt)/T))
        t_stnum = stnum-n_stnum*T
        n_int = n_stnum+1              
        m1 = Mea[pick-1,:]
        m2 = Mea[pick-2,:]
        m3 = Mea[pick-3,:]
        m4 = Mea[pick-4,:]
        m5 = Mea[pick-5,:]
        m6 = Mea[pick-6,:]
        m7 = Mea[pick-7,:]
        m8 = Mea[pick-8,:]
        m9 = Mea[pick-9,:]
        m10 = Mea[pick-10,:]
        m11 = Mea[pick-11,:]
        m12 = Mea[pick-16,:]
        Lf1= m1[17]+m1[18]+m1[19]+m1[20];Lf2= m2[17]+m2[18]+m2[19]+m2[20];Lf3= m3[17]+m3[18]+m3[19]+m3[20];
        Lf4= m4[17]+m4[18]+m4[19]+m4[20];Lf5= m5[17]+m5[18]+m5[19]+m5[20];Lf6= m6[17]+m6[18]+m6[19]+m6[20];
        Lf7= m7[17]+m7[18]+m7[19]+m7[20];Lf8= m8[17]+m8[18]+m8[19]+m8[20];
        Lf9= m9[17]+m9[18]+m9[19]+m9[20];Lf10= m10[17]+m10[18]+m10[19]+m10[20];
        Lf11= m11[17]+m11[18]+m11[19]+m11[20];Lf12= m12[18]+m12[19]+m12[19]+m12[20];    
        Rf1= m1[13]+m1[14]+m1[15]+m1[16];Rf2= m2[13]+m2[14]+m2[15]+m2[16];Rf3= m3[13]+m3[14]+m3[15]+m3[16];
        Rf4= m4[13]+m4[14]+m4[15]+m4[16];Rf5= m5[13]+m5[14]+m5[15]+m5[16];Rf6= m6[13]+m6[14]+m6[15]+m6[16];
        Rf7= m7[13]+m7[14]+m7[15]+m7[16];Rf8= m8[13]+m8[14]+m8[15]+m8[16];
        Rf9= m9[13]+m9[14]+m9[15]+m9[16];Rf10= m10[13]+m10[14]+m10[15]+m10[16];    
        Rf11= m11[13]+m11[14]+m11[15]+m11[16];Rf12= m12[13]+m12[14]+m12[15]+m12[16];          

        eps =0.001
        flagx[pick-1,0] = n_int
        flagx1 = flagx[(pick-100):(pick),0]
        flagx1 = np.reshape(flagx1,[1,100])
        
        flagx2=flagx[(pick-100):(pick),1]
        flagxn=0
        anx = find0(flagx1,n_int,2,100)
        ana =len(anx)
        for i in range(0,ana):
            x = anx[i]
            flagxn=flagxn+abs(flagx2[x])

####      无调整，有预观
##            
##        flagx[pick-1,1] = 0
##        [px,py,pz] = optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
##        
##    %判断是否已经调整过当前周期，确保只进行一次调整；
##     根据实际情况，应该对足底接触力的条件进行一次修正。这里感觉对每只脚分别判断会比较好。但是阈值不好处理。
        if (flagxn>=1):
            flagx[pick-1,1] = 0
            [px,py,pz] = optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)    
        else:            
            if ((Lf1*Rf1>eps)and (Lf2*Rf2<=eps)and(Lf3*Rf3<=eps)and(Lf4*Rf4<=eps)and(Lf5*Rf5<=eps)):                
                if (t_stnum<=(T-t_land[n_int-1,:]))and (t_stnum>=DSP[n_int-1,:]):
                    if np.remainder(n_stnum,2)==0 :
                        z_hx = outx[11]
                    else:
                        z_hx = outx[14]

                    if ((z_hx<0) and (t_stnum> 3*T/4)and (outx[16]>0)):
                        flagx[pick-1,1] = 1
                        DSP[n_int,:] = 3*DSP[n_int-1,:]/2                           
                        SSP[n_int,:] = T-DSP[n_int,:]
                        
                        if DSP[n_int,0]<=0.05*T:
                            DSP[n_int,0]= 0.05*T
                        else:
                            if DSP[n_int,0]>=0.5*T:
                                DSP[n_int,0]= 0.5*T
                                
                        SSP[n_int,0] = T - DSP[n_int,0]
                ##      接下来所有周期
                        for i in range(n_int+1,NNx):
                           SSP[i,0] = SSP[n_int,0]
                           DSP[i,0] = DSP[n_int,0]                                                
                                
                        [px,py,pz] = optimal_downslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
                    else:
                        if (flag[n_stnum-2,:]==1)and(flag[n_stnum-3,:]==1):
                            flagx[pick-1,1] = 0
                            [px,py,pz] = optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
                        else:
                            flagx[pick-1,1] = -1
                            [px,py,pz] = optimal_upslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
                else:
                    if (t_stnum>=(T-1*t_land[n_int-1,:]/2))or(t_stnum<=t_off[n_int-1,:]):
                        if np.remainder(n_stnum,2)==0 :
                            z_hx = outx[11]
                        else:
                            z_hx = outx[14]

##                        if ((z_hx>0)and(outx[16]<0)):
                        if ((z_hx>0)and(outx[16]<0)):                            
                            flagx[pick-1,1] = -1
                            [px,py,pz] = optimal_upslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)                            
                            
                        else:
                            if (flag[n_stnum-2,:]==-1)and(flag[n_stnum-3,:]==-1):
                                flagx[pick-1,1] = 0
                                [px,py,pz] = optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)                            
                            else:
                                flagx[pick-1,1] = 1
                                if (t_stnum>=(T-1*t_land[n_int-1,:]/2)):
                                    DSP[n_int,:] = 3*DSP[n_int-1,:]/2
                                    SSP[n_int,:] = T-DSP[n_int,:]                                
                                else:
                                    DSP[n_int,:] = DSP[n_int-1,:]-t_stnum
                                    SSP[n_int,:] = T-DSP[n_int,:]
                                            
                                if DSP[n_int,0]<=0.05*T:
                                    DSP[n_int,0]= 0.05*T
                                else:
                                    if DSP[n_int,0]>=0.5*T:
                                        DSP[n_int,0]= 0.5*T
                                                
                                SSP[n_int,0] = T - DSP[n_int,0]
                        ##      接下来所有周期
                                for i in range(n_int+1,NNx):
                                   SSP[i,0] = SSP[n_int,0]
                                   DSP[i,0] = DSP[n_int,0]                       
                                                
                                [px,py,pz] = optimal_downslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
                    else:
                        flagx[pick-1,1] = 0
                        [px,py,pz] = optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
            else:               
                if (((t_stnum==T)or (t_stnum==T-dt))and(Lf1*Rf1<=eps)and (Lf2*Rf2<=eps)and(Lf3*Rf3<=eps)and(Lf4*Rf4<=eps)and(Lf5*Rf5<=eps)
                    and(Lf6*Rf6<=eps)and(Lf7*Rf7<=eps)and(Lf8*Rf8<=eps)and(Lf9*Rf9<=eps)and(Lf10*Rf10<=eps)):

                    if np.remainder(n_stnum,2)==0 :
                        z_hx = outx[11]
                    else:
                        z_hx = outx[14]
                    if ((z_hx>0)and(outx[16]<0)):                            
                        flagx[pick-1,1] = -1
                        [px,py,pz] = optimal_upslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)                            
                            
                    else:                    
                        if (flag[n_stnum-2,:]==-1)and(flag[n_stnum-3,:]==-1):
                            flagx[pick-1,1] = 0
                            [px,py,pz] = optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
                        else:
                            flagx[pick-1,1] = 1
                            DSP[n_int,:] = 3*DSP[n_int-1,:]/2
                            SSP[n_int,:] = T-DSP[n_int,:]
                            
                            if DSP[n_int,0]<=0.05*T:
                                DSP[n_int,0]= 0.05*T
                            else:
                                if DSP[n_int,0]>=0.5*T:
                                    DSP[n_int,0]= 0.5*T
                                    
                            SSP[n_int,0] = T - DSP[n_int,0]
                        ##      接下来所有周期
                            for i in range(n_int+1,NNx):
                                SSP[i,0] = SSP[n_int,0]
                                DSP[i,0] = DSP[n_int,0]                                    
                            [px,py,pz] = optimal_downslope(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)
                else:
                    flagx[pick-1,1] = 0
                    [px,py,pz] = optimal_flatground(n_stnum,t_stnum,n_intx,n_int,Tc,dt,t,Zc,g,outx,pick,nx,T,length,width,XN)                       

##   异常处理
    nn = len(pz)

##    for i in range(0,len(tempz)):
##        pz[tempz[i]] = pz[tempz[i]-1]

    tempz = find0(np.array([pz]),0,-1,nn)    
    for i in range(0,len(tempz)):
        if tempz[i]==0:
            pz[tempz[i]] = pz[tempz[i]+1]
        else:
            pz[tempz[i]] = pz[tempz[i]-1]
            
    tempx = find0(np.array([px]),0,-1,nn)
    for i in range(0,len(tempx)):
        if tempx[i]==0:
            px[tempx[i]] = px[tempx[i]+1]
        else:
            px[tempx[i]] = px[tempx[i]-1]
            
    tempy = find0(np.array([py]),0,-1,nn)
    for i in range(0,len(tempy)):
        if tempy[i]==0:
            py[tempy[i]] = py[tempy[i]+1]
        else:
            py[tempy[i]] = py[tempy[i]-1]    

    px1 = list(px)
    py1 = list(py)
    pz1 = list(pz)
    
    [comx,comy,comz,tsh] = ZMP2COM(t,dt,Zc,g,px1,py1,pz1)

    return comx,comy,comz,tsh
##-------------------------------------------------------------------------
    
    

##--------------------预观控制（ZMP跟踪控制）
def bp_previewcontrol(stnum,N0,t,p_real,p_ref,picst,picen,picst1,picen1,bodyp_1,dt,x):
    global k
    global F_pre
    global sumde1_p
        
    N1 = np.size(t)
    N1x = picen1 - picst1 +1
    U = local_zeros([N1x,3])
##  该函数应该离线计算好
##    [Ks,Kx,F,L,Ao,b,C]=calclkf_new(dt,N0)

    Ks =  -25.9144
    Kx = local_array([-666.3026,-121.4437,119.2303])
    F = F_pre[0:N0,:]
##    F = np.reshape(F,[])
    L = local_array([0.0828,-0.1666,1.0000])
    Ao = local_array([[1.0000,0.0200,0],[0.6323,1.0000,-0.6323],[0,0,1.0000]])
    b = local_array([[0],[0],[0.02]])
    C = local_array([0,0,1])
   
    del_p = p_ref - p_real
    X = local_ones([3*N1,3])
    sump_ref = local_zeros([1,3])   

    for k in range(1,N1+1):
        X[(3*k-3):(3*k),:] = local_array([x[picst1+k-2,:],(x[picst1+k-2,:]-x[picst1+k-3,:])/dt,p_ref[picst1+k-2,:]])
        if stnum == 12+dt:
            for i1 in range(1,picst1+k+1):
                sumde1_p = sumde1_p + del_p[i1 - 1,:]
        else:
            sumde1_p = sumde1_p + del_p[picst1,:]
            
        for j1 in range(picst1+k,picst1+k+N0):
            sump_ref = sump_ref + F[j1-(picst1+k)] * p_ref[j1-1,:]

        U[picst+k-1,:] =Ks*sumde1_p - local_dot(Kx, X[(3*k-3):(3*k),:]) + sump_ref

    bodypx = local_ones([3*N1x,3])
    bodypx[0:3,:] = X[0:3,:]
    bodyp = bodyp_1

    for k in range(1,N1+1):
        bodypx[3*k:(3*k+3),:] = local_dot(Ao,X[(3*k-3):3*k,:])-local_dot(L.T,(p_real[picst1+k-2,:]-local_dot(C,X[(3*k-3):3*k,:])))+b*U[picst+k-1,:]

    for k in range(1,N1+1):
        bodyp[k-1,:] = bodypx[3*k-3,:]
        
##  bodyp为5*3矩阵
    return bodyp
##----------------------------------------------------------------------------------


##  -----------------------------------规划出足部轨迹-----------------------------
def FootpR(dt,t_des,N_foot,stepwidth,T,pick):
    global T_initial 
    global DSP
    global SSP
    global slope_angle
    global Sy
    global footx
    global footy
    global footz
    global footy2
    global flagxn     
    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones
    
    T_ini = T_initial
    b1 = 48.3
##    %前半足长度。
    b2 = 102.3
    NX = footy2.shape[0]
    #落地时，足于地面的夹角
    th10 = 2 * np.pi / 180
    th1 = 2 * np.pi / 180 * local_ones([NX,1]) 
    
    #离地时，足与地面的夹角。
    th20 = 2 * np.pi / 180
    th2 = 2 * np.pi / 180 * local_ones([NX,1]) 

    #落地时间
    t1 = 0.1
    #离地时间
    t2 = 0.1

    t_th1 = 0.1 * local_ones([NX,1])
    t_th2 = 0.1 * local_ones([NX,1])
    h = 20
    h0 = 20;

##    #调整步宽
##  这里不进行调整  
    for i in range(3,NX):       
##        footy2[i-1,0] = footy[i-1,0]+(-1)**(i-1)*5
        footy2[i-1,0] = footy[i-1,0]
        
##        footy2[i-1,0] = footy[i-1,0]        
    ## /*注意position_ZMP_COM_energymin_footplace.m
    ## 足部落脚的序列中第三个点是左脚的第一步落脚点，第四个点是右脚的第二步落脚点 */

    ## 注意这里的下标与
    i= N_foot
    if DSP[i-2,0]<=0:
        DSP[i-2,0] = 0.05*T
    if DSP[i-1,0]<=0:
        DSP[i-1,0] = 0.05*T
    if DSP[i,0]<=0:
        DSP[i,0] = 0.05*T
    SSP[i-1,0] = T-DSP[i-1,0]
    SSP[i-2,0] = T-DSP[i-2,0]
    SSP[i,0] = T-DSP[i,0]
##      考虑上下坡模式的切换：下坡要加速

    ax = local_array([[0],[1],[0]])
    bx = local_array([[0],[0],[-45.19]])
        
    if np.remainder(i,2) != 0:     #奇数为左脚支撑，偶数为右脚支撑
        Rgai_t =  local_array([T_ini+local_floor((2*i-3)/2)*T, T_ini+local_floor((2*i-3)/2)*T+DSP[i-2,0], T_ini+local_floor((2*i-3)/2)*T+DSP[i-2,0]+t_th2[i-2,0],
                                      T_ini+local_floor((2*i-3)/2)*T+DSP[i-2,0]+SSP[i-2,0]/2, T_ini+local_floor((2*i-1)/2)*T-t_th1[i-1,0], T_ini+local_floor((2*i-1)/2)*T,
                                      T_ini+local_floor((2*i-1)/2)*T+DSP[i-1,0] ])            
        if (footx[i,0] - footx[i-1,0]>0):
            Ra = local_array([0+slope_angle[i-2,0],0+slope_angle[i-2,0],th2[i,0]+slope_angle[i-2,0],th2[i,0]/2,-th1[i,0]+slope_angle[i,0],0+slope_angle[i,0],0+slope_angle[i,0]])                
        else:
            if (footx[i,0] - footx[i-1,0]==0):
                if (footx[i-1,0] - footx[i-2,0]==0):
                    Ra =  local_array([0+slope_angle[i-2,0],0+slope_angle[i-2,0],0+slope_angle[i-2,0],0,0+slope_angle[i,0],0+slope_angle[i,0],0+slope_angle[i,0]])
                else:
                    if (footx[i-1,0] - footx[i-2,0]>0):
                        Ra =  local_array([+slope_angle[i-2,0],+slope_angle[i-2,0],th2[i,0]+slope_angle[i-2,0],th2[i,0]/2,-th1[i,0]+slope_angle[i,0],+slope_angle[i,0],+slope_angle[i,0]])
                    else:
                        th2[i,0]=1* np.pi/180
                        th1[i,0]=1* np.pi/180
                        Ra =  local_array([+slope_angle[i-2,0],+slope_angle[i-2,0],-th2[i,0]+slope_angle[i-2,0],-th2[i,0]/2,th1[i,0]+slope_angle[i,0],+slope_angle[i,0],+slope_angle[i,0]])           
            else:
                th2[i,0]=1* np.pi/180
                th1[i,0]=1* np.pi/180
                Ra =  local_array([+slope_angle[i-2,0],+slope_angle[i-2,0],-th2[i,0]+slope_angle[i-2,0],-th2[i,0]/2,th1[i,0]+slope_angle[i,0],+slope_angle[i,0],+slope_angle[i,0]])
        if (flagx[pick-1,1] <=0):
            Rgai_px =  local_array([footx[i-2,0],footx[i-2,0],footx[i-2,0]+ b2*(1-local_cos(Ra[-5])),footx[i-1,0],footx[i,0]-b1*(1-local_cos(Ra[-3])),footx[i,0],footx[i,0]])
            Rgai_py =  local_array([footy2[i-2,0],footy2[i-2,0],footy2[i-2,0],(footy2[i,0]+footy2[i-2,0])/2, footy2[i,0],footy2[i,0],footy2[i,0]])
            Rgai_pz =  local_array([footz[i-2,0],footz[i-2,0],footz[i-2,0]+b2*local_sin(abs(Ra[-5])),h+(footz[i-2,0]+footz[i,0])/2, b1*local_sin(abs(Ra[-3]))+footz[i,0],+footz[i,0],+footz[i,0]])
##下坡加速========================
        else:
            Rgai_px =  local_array([footx[i-2,0],footx[i-2,0],footx[i-2,0]+ b2*(1-local_cos(Ra[-5])),(footx[i-1,0]+3*footx[i,0])/4,footx[i,0]-b1*(1-local_cos(Ra[-3])),footx[i,0],footx[i,0]])
            Rgai_py =  local_array([footy2[i-2,0],footy2[i-2,0],footy2[i-2,0],(footy2[i-2,0]+3*footy2[i,0])/4, footy2[i,0],footy2[i,0],footy2[i,0]])
            Rgai_pz =  local_array([footz[i-2,0],footz[i-2,0],footz[i-2,0]+b2*local_sin(abs(Ra[-5])),h+(footz[i-2,0]+footz[i,0])/2, b1*local_sin(abs(Ra[-3]))+footz[i,0],+footz[i,0],+footz[i,0]])

        Rpbmx = itp.pchip(Rgai_t,Rgai_px)(t_des)
        Rpbmy = itp.pchip(Rgai_t,Rgai_py)(t_des)
        Rpbmz = itp.pchip(Rgai_t,Rgai_pz)(t_des)
        Raa =   itp.pchip(Rgai_t,Ra)(t_des)
        Rfoot = local_array([[Rpbmx],[Rpbmy],[Rpbmz]])

        RR = Rodrigues_x_once(ax,Raa,local_sin,local_cos,local_dot,uni_matrix)
        Rp = Rfoot- local_dot(RR,bx)

        LR = Rodrigues_x_once(ax,slope_angle[i-1,0],local_sin,local_cos,local_dot,uni_matrix)
        Lp = local_array([[footx[i-1,0]],[footy2[i-1,0]],[footz[i-1,0]]])- local_dot(LR,bx)
            
    else:     #右腿支撑
        Lgai_t =  local_array([T_ini+local_floor((2*i-3)/2)*T, T_ini+local_floor((2*i-3)/2)*T+DSP[i-2,0], T_ini+local_floor((2*i-3)/2)*T+DSP[i-2,0]+t_th2[i-2,0],
                                  T_ini+local_floor((2*i-3)/2)*T+DSP[i-2,0]+SSP[i-2,0]/2, T_ini+local_floor((2*i-1)/2)*T-t_th1[i-1,0], T_ini+local_floor((2*i-1)/2)*T,
                                  T_ini+local_floor((2*i-1)/2)*T+DSP[i-1,0] ])
        if (footx[i,0] - footx[i-1,0]>0):
            La =  local_array([0+slope_angle[i-2,0],0+slope_angle[i-2,0],th2[i,0]+slope_angle[i-2,0],th2[i,0]/2,-th1[i,0]+slope_angle[i,0],0+slope_angle[i,0],0+slope_angle[i,0]])
        else:
            if (footx[i,0] - footx[i-1,0]==0):
                if (footx[i-1,0] - footx[i-2,0]==0):
                    La =  local_array([0+slope_angle[i-2,0],0+slope_angle[i-2,0],0+slope_angle[i-2,0],0,0+slope_angle[i,0],0+slope_angle[i,0],0+slope_angle[i,0]])
                else:
                   if (footx[i-1,0] - footx[i-2,0]>0):
                       La =  local_array([+slope_angle[i-2,0],+slope_angle[i-2,0],th2[i,0]+slope_angle[i-2,0],th2[i,0]/2,-th1[i,0]+slope_angle[i,0],+slope_angle[i,0],+slope_angle[i,0]])
                   else:
                        th2[i,0]=1* np.pi/180
                        th1[i,0]=1* np.pi/180
                        La =  local_array([+slope_angle[i-2,0],+slope_angle[i-2,0],-th2[i,0]+slope_angle[i-2,0],-th2[i,0]/2,th1[i,0]+slope_angle[i,0],+slope_angle[i,0],+slope_angle[i,0]])
            else:
                th2[i,0]=1* np.pi/180
                th1[i,0]=1* np.pi/180                            
                La =  local_array([+slope_angle[i-2,0],+slope_angle[i-2,0],-th2[i,0]+slope_angle[i-2,0],-th2[i,0]/2,th1[i,0]+slope_angle[i,0],+slope_angle[i,0],+slope_angle[i,0]])
        if (flagx[pick-1,1] <=0):           
            Lgai_px =  local_array([footx[i-2,0],footx[i-2,0],footx[i-2,0]+ b2*(1-local_cos(La[-5])),footx[i-1,0],footx[i,0]-b1*(1-local_cos(La[-3])),footx[i,0],footx[i,0]])
            Lgai_py =  local_array([footy2[i-2,0],footy2[i-2,0],footy2[i-2,0],(footy2[i,0]+footy2[i-2,0])/2, footy2[i,0],footy2[i,0],footy2[i,0]])
            Lgai_pz =  local_array([footz[i-2,0],footz[i-2,0],footz[i-2,0]+b2*local_sin(abs(La[-5])),h+(footz[i-2,0]+footz[i,0])/2, b1*local_sin(abs(La[-3]))+footz[i,0],+footz[i,0],+footz[i,0]])
        else:
            Lgai_px =  local_array([footx[i-2,0],footx[i-2,0],footx[i-2,0]+ b2*(1-local_cos(La[-5])),(footx[i-1,0]+3*footx[i,0])/4,footx[i,0]-b1*(1-local_cos(La[-3])),footx[i,0],footx[i,0]])
            Lgai_py =  local_array([footy2[i-2,0],footy2[i-2,0],footy2[i-2,0],(footy2[i-2,0]+3*footy2[i,0])/4, footy2[i,0],footy2[i,0],footy2[i,0]])
            Lgai_pz =  local_array([footz[i-2,0],footz[i-2,0],footz[i-2,0]+b2*local_sin(abs(La[-5])),h+(footz[i-2,0]+footz[i,0])/2, b1*local_sin(abs(La[-3]))+footz[i,0],+footz[i,0],+footz[i,0]])       
                                 

        Lpbmx = itp.pchip(Lgai_t,Lgai_px)(t_des)
        Lpbmy = itp.pchip(Lgai_t,Lgai_py)(t_des)
        Lpbmz = itp.pchip(Lgai_t,Lgai_pz)(t_des)
        Laa =   itp.pchip(Lgai_t,La)(t_des)
        Lfoot = local_array([[Lpbmx],[Lpbmy],[Lpbmz]])    

        LR = Rodrigues_x_once(ax,Laa,local_sin,local_cos,local_dot,uni_matrix)   
        Lp = Lfoot- local_dot(LR,bx)

        RR = Rodrigues_x_once(ax,slope_angle[i-1,0],local_sin,local_cos,local_dot,uni_matrix)
        Rp = local_array([[footx[i-1,0]],[footy2[i-1,0]],[footz[i-1,0]]])- local_dot(RR,bx)
        
    return Rp,RR,Lp,LR   
##-----------------------------------------------------------------


                        
##------------------------ 在线足部轨迹和com轨迹生成-------------------------------------
def plan_com_foot_online(stnum,txx,txx1,tN1,dt,g,Zc,stepwidth,outx,pick,xn,t1,picst1,picen1,picst2,picen2,N3,T):
    global T_initial
    global zmp_ref
    global zmp_real
    global com_ref
    global com_real
    global PX1
    global PY1
    global PZ1

    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones
    
    Tc = np.sqrt(Zc/g)
    N_foot = int(local_floor((txx1[-1]-dt)*100/T/100))+1  #当前支撑时刻的下一周期       
    [comx,comy,comz,tsh]=position_ZMP_COM_energymin_footplace(Tc,txx,dt,txx1,tN1,N_foot,Zc,g,outx,pick,T,xn)
      
    n_comx = comx.shape[0]
    nn = int(np.floor(n_comx/2))
    comx1=comx[nn,0] 
    comy1=comy[nn,0]
    comz1=comz[nn,0]
##  并不使用，实时输出用于分析。
    px1=PX1[xn,0]
    py1=PY1[xn,0]
    pz1=PZ1[xn,0]

##  预观控制
########
    tnn = tN1 - tsh
    zmp_real1 = zmp_real[0:N3,:]
##    zmp_real1 = zmp_real1[0:-tsh/2,:]   
    zmp_ref1 = zmp_ref[0:N3,:]

    com_ref1 = com_ref[0:N3,:]
    com_ref1 = com_ref1[0:-tsh/2,:]
    
    pxx = list(PX1)
    pyy = list(PY1)
    pzz = list(PZ1)
    for i in range(0,tN1):
        zmp_ref1[-tN1+i,0] = pxx[i]
        zmp_ref1[-tN1+i,1] = pyy[i]
        zmp_ref1[-tN1+i,2] = pzz[i]
        
    for i in range(0,tnn):
        com_ref1[-tnn+i,0] = comx[i,0]
        com_ref1[-tnn+i,1] = comy[i,0]
        com_ref1[-tnn+i,2] = comz[i,0]       
        
##  目测好像应该是使用com_real，但是发现跳动比较大
    bodyp_1 = com_ref1[picst2:picen2+1,:]
    
####  注意，这里的picst已经是实际的下标了，所以这里要进行改变。
    picst1 = picst1+1
    picen1 = picen1+1
    picst2 = picst2+1
    picen2 = picen2+1
    xn1= int(local_floor((tN1-tsh)/2))
    bpx = bp_previewcontrol(stnum,xn1,t1,zmp_real1,zmp_ref1,picst1,picen1,picst2,picen2,bodyp_1,dt,com_ref1)
    com_tar = bpx[2,:]
    comx11 = com_tar[0]
    comy11 = com_tar[1]   
##    %角度系数运用
    comx1 = (comx1+comx11)/2
    comy1 = (comy1+comy11)/2

##  足部轨迹生活  
    n_foot = int(local_floor((txx-dt)/T))+1 #当前时刻对应的周期。
    stnum_foot = pick*dt
    [Rp,RR,Lp,LR] = FootpR(dt,stnum_foot,n_foot,stepwidth,T,pick)
    
    return comx1,comy1,comz1,px1,py1,pz1,Lp,LR,Rp,RR

####===============================================================================



####-------------------------------------------------------------------------------
###======================逆运动学计算
 
##计算雅可比矩阵
def CalJacobian(ray,local_dot,mylink1):
    
    x = ray[5]
    pn=mylink1[x-1].p
    J=local_zeros((6,5))
    
    for i,e in enumerate(ray[1:]):
        a = local_dot(mylink1[e-1].R,mylink1[e-1].a)
        b = pn-mylink1[e-1].p
        cr = [a[1,0]*b[2,0] - a[2,0]*b[1,0], a[2,0]*b[0,0] - a[0,0]*b[2,0], a[0,0]*b[1,0] - a[1,0]*b[0,0]]
        J[:,i]= local_array([cr[0],cr[1],cr[2],a[0,0],a[1,0],a[2,0]]) 
        
    return J

##姿态到角速度的转换
def R2W(dR,local_sin,uni_matrix):
    if np.linalg.norm (dR - uni_matrix) >= 0.000001:
        tht = np.sum(np.diag(dR))-1
        thtx = tht*(0.5)
        if abs(thtx)>1:
            if thtx >0:
                thtx = 0.999999999999999
            else:
                thtx = -0.999999999999999
        thta = math.acos(thtx)
        try:
            w = thta/(2 * local_sin(thta)) * local_array([[dR[2,1]-dR[1,2]],[dR[0,2]-dR[2,0]],[dR[1,0]-dR[0,1]]])
        except ZeroDivisionError as e:
            w = local_zeros((3,1))        
    else :
        w = local_zeros((3,1))

    return w

##计算位姿误差
def CalErr(TargetID,PosRef,local_dot,local_sin,uni_matrix,mylink1):
    linkx = mylink1[TargetID-1]
    dp = PosRef.p - linkx.p    
    rr = (linkx.R).T
    
    XR = local_dot(rr,PosRef.R)
    
    XRES = R2W(XR,local_sin,uni_matrix)
    
    dw = local_dot(linkx.R, XRES)
    
    err = local_array([dp[0,0],dp[1,0],dp[2,0],dw[0,0],dw[1,0],dw[2,0]]) 

    return err


##逆运动学循环求解
def Inverse_Kinematics_x(TargetID, PosRef,local_sin,local_cos,local_dot,uni_matrix,mylink1):
    ##目标干只有联众情况，这里就直接给出
    if TargetID == 6:
        ray = [1,2,3,4,5,6]
    else :
        ray = [1,7,8,9,10,11]

##  6次循环:6次误差比较大，至少8次：系数为0.65就容易跳变了  
##    lam = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
##    n_loop = 15
##  5次循环  
##    lam = [0.75,0.75,0.75,0.75,0.75,0.75,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
##    n_loop = 5
##  4次循环    
##    lam = [0.75,0.65,0.6,0.75,0.70,0.7,0.7,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
##    n_loop = 3
####    上一次角度为初值，循环3次即可达到无初值15次的结果。   
    lam = [0.55,0.55,0.55,0.75,0.70,0.7,0.7,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
    n_loop = 3
        
    for i in range(0,n_loop):        
        Forward_Kinematics_x(2,local_sin,local_cos,local_dot,uni_matrix,mylink1)        
        J = CalJacobian(ray,local_dot,mylink1)        
        err = CalErr(TargetID,PosRef,local_dot,local_sin,uni_matrix,mylink1)        
        if np.linalg.norm (err) < 0.0001:
            break            
        try:
            invers_J = np.linalg.pinv(J)
        except:
            break  
        delta_q = lam[i] *local_dot(invers_J,err)
        
        for i,e in enumerate(ray[1:]):
            mylink1[e-1].q = mylink1[e-1].q +delta_q[i]

    


####===============================================================================
####======================================main函数=================================
def nao_unknown_slope_previewcontrol(pick,mea):
    global Link
    global Mea
    global angle
    global flag
    global wf
    global zmp_g
    global com_g
    global rfoot_g
    global lfoot_g
    global n_period
    global T_initial
    global DSP
    global SSP
    global R_sensor
    global L_sensor

    global Sx 
    global Sy
    global Sz 
    global x_init
    global x_end
    global Vx_init
    global Vx_end 
    global y_init
    global y_end
    global Vy_init
    global Vy_end 
    global PX1
    global PY1
    global PZ1
    global footx
    global footy
    global footz
    global z_h 
    global slope_angle   
    global f_angle_roll 
    global f_angle_pitch 
    global f_angle_yaw 
##    %%%落地和离地时间
    global t_land 
    global t_off  
    global flagx 

##    %预观控制需要的存储器。
    global zmp_ref 
    global zmp_real 
    global com_real 
    global com_ref

    global data_anglex

    global local_sin
    global local_cos
    global local_dot
    global local_floor
    global uni_matrix
    global local_array 
    global local_zeros
    global local_sinh
    global local_cosh
    global local_tanh
    global local_ones     

    
    NN=40
    ##    总时间
    dt= 0.02
    T=2
    Tsup = np.float32(1.6)
    Tdou = np.float32(0.4)
    stepwidth = 100.0
    stepx  = 40.0
    g = 9800.0
    Zc= 310.0

##--这里pick从1开始
    palse =1
    stnum = pick * dt

    nxxx = int(local_floor((stnum+2*T)*100/dt/100))
    tx = np.linspace(dt,stnum+2*T,nxxx)    
    if stnum<= T:
        t = local_array([stnum,stnum+dt])
    else:
        t = local_array([stnum-dt,stnum,stnum+dt])
    if stnum<= T:
        picstx = int(local_floor(stnum*100/dt/100))
        picenx = picstx + 1   
    else:
        picstx = int(local_floor((stnum-dt)*100/dt/100))
        picenx = picstx + 2         
    
    t1 = local_array([stnum+dt,stnum+2*dt])
    if stnum<= dt:
        t2 = local_array([stnum+dt,stnum+2*dt,stnum+3*dt,stnum+4*dt,stnum+5*dt])
    else:
        t2 = local_array([stnum-dt,stnum,stnum+dt,stnum+2*dt,stnum+3*dt])

    N3 = int(local_floor((stnum+T-dt)*100/dt/100))
    t3 = np.linspace(dt,stnum+T-dt,N3)    
        
    if stnum<= dt:
        picst1 = 0
        picen1 = picst1 + 1
        picst2 = int(local_floor((stnum+dt)*100/dt/100)-1)
        picen2 = picst2 + 4        
    else:
        picst1 = 2
        picen1 = picst1 + 1
        picst2 = int(local_floor((stnum-dt)*100/dt/100)-1)
        picen2 = picst2 + 4 

##    LZMPy0 = 55
    txx = stnum-(T_initial)+T
    tN1 = int(local_floor((2*T-dt)*100/dt/100))

    txx1 = np.linspace(txx-T+dt,txx+T-dt,tN1)

    xn = int(local_floor(tN1/2))

##%当前时刻对应的步行周期
#### =========================================================   
    xtn = find0(n_period,stnum,0,40)
    xtn = xtn[0]
    
    Mea[pick-1,:] = mea
##  反馈量计算
##  nao实际使用 
    q = local_array([0,mea[3],mea[4],mea[5],mea[6],mea[7],mea[8],mea[9],mea[10],mea[11],mea[12]])
##  注意这里x,y,z分别为+，+，-，且都是绝对坐标        
##    x_angle = mea[0]         
##    y_angle = mea[1]        
##    z_angle = -mea[2]
    
    if pick >= 3:
####      修正三
        x_angle = mea[0] - (Mea[1,0])/2      
        y_angle = mea[1] - (Mea[1,1])/2    
        z_angle = -mea[2]- (-Mea[1,2])
        
    else:
        x_angle = mea[0]         
        y_angle = mea[1]        
        z_angle = -mea[2]        
    
    Rrpy = local_array([[local_cos(y_angle)*local_cos(z_angle),local_cos(z_angle)*local_sin(x_angle)*local_sin(y_angle) - local_cos(x_angle)*local_sin(z_angle),
                         local_sin(x_angle)*local_sin(z_angle) + local_cos(x_angle)*local_cos(z_angle)*local_sin(y_angle)],
                        [local_cos(y_angle)*local_sin(z_angle),local_cos(x_angle)*local_cos(z_angle) + local_sin(x_angle)*local_sin(y_angle)*local_sin(z_angle),
                         local_cos(x_angle)*local_sin(y_angle)*local_sin(z_angle) - local_cos(z_angle)*local_sin(x_angle)],
                        [-local_sin(y_angle),local_cos(y_angle)*local_sin(x_angle),local_cos(x_angle)*local_cos(y_angle)]])    

    Initial_once()
    [R_RefpR,L_RefpR] = posRef()
    R_RefpR.v = np.random.rand(3,1)
    R_RefpR.w = np.random.rand(3,1)
    L_RefpR.v = np.random.rand(3,1)
    L_RefpR.w = np.random.rand(3,1)
    
    mylink1 = mylink    
    for i,e in enumerate(q):
        mylink1[i].q = e
    
    mylink1[0].v = np.random.rand(3,1)
    mylink1[0].w = np.random.rand(3,1)
    mylink1[0].p = local_array([[0],[0],[0]])
    mylink1[0].R = Rrpy
    
    real_com = Rea1_com(mylink1)             ##计算质心

    Lf = mea[17]+mea[18]+mea[19]+mea[20]
    Rf = mea[13]+mea[14]+mea[15]+mea[16]
    [wff,zmpf,comf,rfootf,lfootf,slope_estimation,L11,L21,L31,L41,R11,R21,R31,R41] = cal_foot_zmp_nao(pick,xtn,mea,real_com,Lf,Rf,mylink1)

    outx = local_array([wff[0],wff[1],wff[2],zmpf[0],zmpf[1],zmpf[2],comf[0],comf[1],comf[2],rfootf[0],rfootf[1],rfootf[2],
                    lfootf[0],lfootf[1],lfootf[2],slope_estimation[0],slope_estimation[1],slope_estimation[2],
                    L11[0,0],L11[0,1],L11[0,2],L21[0,0],L21[0,1],L21[0,2],L31[0,0],L31[0,1],L31[0,2],L41[0,0],L41[0,1],L41[0,2],
                    R11[0,0],R11[0,1],R11[0,2],R21[0,0],R21[0,1],R21[0,2],R31[0,0],R31[0,1],R31[0,2],R41[0,0],R41[0,1],R41[0,2]])

    zmp_real[pick-1,:] = zmpf+wff
    com_real[pick-1,:] = comf+wff
##====================================================================================================
##  ====================================================== 
####步态生成算法
    q= local_array([0, 0, 0, 0.2, 0, 0, 0, 0, 0.2, 0, 0])
    for i in range(0,11):
        mylink1[i].q = q[i]
        
    mylink1[0].R = np.eye(3)

    if palse ==1:
#      注意，这里周期为2s,要进行修正.前三个周期要使用       
#        if stnum<=9:
        if stnum<=12:
            y =data_anglex[picstx-1:picenx,:]
            
            footx1=footx; footy1=footy; footz1=footz; SSP1=SSP; DSP1=DSP;
            t_th1=t_land; thx=slope_angle; Sx1=Sx;Sy1=Sy;Sz1=Sz;
            comx=0
            comy=0
            comz=0
            px1= PX1[99,0]
            py1= PY1[99,0]
            pz1= PZ1[99,0]
            
            Lp= local_zeros([3,1])
            LR=np.eye(3)
            Rp=local_zeros([3,1])
            RR=np.eye(3)
            if stnum<=T:
                angle[pick-1,:] = y[0,:];        
            else:
                if stnum<=12:
                    angle[pick-1,:] = y[1,:]
            
        else:
            [comx,comy,comz,px1,py1,pz1,Lp,LR,Rp,RR] = plan_com_foot_online(stnum,txx,txx1,tN1,dt,g,Zc,stepwidth,outx,pick,xn,t1,picst1,picen1,picst2,picen2,N3,T)

######可能导致程序失败，尤其是逆运动学解出现问题。
            for i in range(1,11):
                if pick<=750:
                    mylink1[i].q = data_anglex[pick-2,i-1]
                else:
                    mylink1[i].q = (angle[pick-2,i-1])

            footx1=footx; footy1=footy; footz1=footz; SSP1=SSP; DSP1=DSP;t_th1=t_land; thx=slope_angle; Sx1=Sx;Sy1=Sy;
            
            bp = local_array([[comx],[comy],[comz]])
            mylink1[0].p = bp
            R_RefpR.p=Rp;   R_RefpR.R=RR;
            L_RefpR.p=Lp;   L_RefpR.R=LR;            
            Inverse_Kinematics_x(6, R_RefpR,local_sin,local_cos,local_dot,uni_matrix,mylink1)
            Inverse_Kinematics_x(11, L_RefpR,local_sin,local_cos,local_dot,uni_matrix,mylink1)
            
            y = local_array([mylink[1].q,mylink[2].q,mylink[3].q,mylink[4].q,mylink[5].q,
                               mylink[6].q,mylink[7].q,mylink[8].q,mylink[9].q,mylink[10].q])
            if pick <750:
                angle[pick-1,:] = ((750-pick)*data_anglex[pick-1,:]+(pick-600)*y)/150
            else:
                angle[pick-1,:] = y            
    else:
        y = local_zeros([1,10])


            
##  角度异常的调理函数
    if pick >600:
        ab = angle[pick-1,:]-angle[pick-2,:]
        for i in range(0,10):
            if abs(ab[i])>=0.015:
                if ab[i]>0:
                    angle[pick-1,i] = angle[pick-2,i] +0.010
                else:
                    angle[pick-1,i] = angle[pick-2,i] -0.010
                    

    anglex = angle[pick-1,:]
    anglex = np.reshape(anglex,[1,10])
    
    footx1=footx1.T; footy1=footy1.T; footz1=footz1.T; SSP1=SSP1.T; DSP1=DSP1.T;
    t_th1=t_th1.T; thx=thx.T; Sx1=Sx1.T; Sy1=Sy1.T; flagx1=flagx[pick-1,0]; flagx2=flagx[pick-1,1];
    flagx1 = np.reshape(flagx1,[1,1])
    flagx2 = np.reshape(flagx2,[1,1])
    px1 = np.reshape(px1,[1,1])
    py1 = np.reshape(py1,[1,1])
    pz1 = np.reshape(pz1,[1,1])
    comx = np.reshape(comx,[1,1])
    comy = np.reshape(comy,[1,1])
    comz = np.reshape(comz,[1,1])
    Rpx = np.reshape(Rp,[1,3])
    Lpx = np.reshape(Lp,[1,3])
   
    outx = np.reshape(outx,[1,42])

    out = np.hstack((anglex,outx,footx1,footy1,footz1,SSP1,DSP1,t_th1,thx,Sx1,Sy1,flagx1,flagx2,px1,py1,pz1,comx,comy,comz,Rpx,Lpx))

    return out
##==========================



##只有关节角度跟踪控制，这里不使用躯干水平控制
def vertical_body_joint_tracking(SensorValue,data0,data1,data2,data3,data4,data5,data6,data7,data8,data9,dt,i,xtn):
##  躯干角度误差
    global roll_err_intergral
    global pitch_err_intergral
    global yaw_err_intergral    
    global body_Lroll                                                      
    global body_Lpitch                                                      
    global body_Rroll                                                      
    global body_Rpitch
    global body_yaw
    global yaw_ref
    
##  关节角度误差  
    global det_rhip_roll_inte
    global det_rhip_pitch_inte
    global det_rknee_pitch_inte
    global det_rankle_pitch_inte
    global det_rankle_roll_inte
    global det_lhip_roll_inte
    global det_lhip_pitch_inte
    global det_lknee_pitch_inte
    global det_lankle_pitch_inte
    global det_lankle_roll_inte   


    
    FL=SensorValue[i][19]+SensorValue[i][20] +SensorValue[i][21] +SensorValue[i][22] 
    FR=SensorValue[i][23]+SensorValue[i][24] +SensorValue[i][25] +SensorValue[i][26]

##  躯干的角度修正
    roll_err_intergral = roll_err_intergral + (0-SensorValue[i-1][10])
    pitch_err_intergral = pitch_err_intergral + (0-SensorValue[i-1][11])

    
##    det_roll = 0.2*(0-SensorValue[i][10])+0.0005*(0-SensorValue[i][10]-(0-SensorValue[i-1][10]))/(dt) + 0*roll_err_intergral
##    det_pitch = 0.005*(0-SensorValue[i][11])+0.0001*(0-SensorValue[i][11]-(0-SensorValue[i-1][11]))/(dt) + 0*pitch_err_intergral
##    det_yaw = 0.005*(0+SensorValue[i][12])+0.0001*(0+SensorValue[i][12]-(0+SensorValue[i-1][12]))/(dt)
    det_roll = 0.1*(0-SensorValue[i][10])+0.0001*(0-SensorValue[i][10]-(0-SensorValue[i-1][10]))/(dt) + 0*roll_err_intergral
    det_pitch = 0.005*(0-SensorValue[i][11])+0.0001*(0-SensorValue[i][11]-(0-SensorValue[i-1][11]))/(dt) + 0*pitch_err_intergral    
    det_yaw = 0.35*(0+SensorValue[i][12])+0.0005*(0+SensorValue[i][12]-(0+SensorValue[i-1][12]))/(dt)
    
    if np.remainder(xtn,2)==0:         ##左脚支撑:
        if FL <0.01:
            if FR <0.001:
                body_Lroll = det_roll
                body_Lpitch = det_pitch                                                      
                body_Rroll = 0                                                       
                body_Rpitch = 0
            else:
                body_Lroll = 0
                body_Lpitch = 0                                                    
                body_Rroll = det_roll                                                      
                body_Rpitch = det_pitch
        else:
                   
            body_Lroll = det_roll
            body_Lpitch = det_pitch                                  
            body_Rroll = 0                                                       
            body_Rpitch = 0                                                    
    else:
        if FR <0.01:
            if FL <0.001:
                body_Lroll = 0
                body_Lpitch = 0                                                    
                body_Rroll = det_roll                                                      
                body_Rpitch = det_pitch
            else:
                body_Lroll = det_roll
                body_Lpitch = det_pitch                                                      
                body_Rroll = 0                                                       
                body_Rpitch = 0
        else:
            body_Lroll = 0
            body_Lpitch = 0                                                      
            body_Rroll = det_roll                                                       
            body_Rpitch = det_pitch
    
    data0[i] = data0[i]+ body_Rroll
    data1[i] = data1[i]+ body_Rpitch
    data5[i] = data5[i]+ body_Lroll
    data6[i] = data6[i]+ body_Lpitch
    yaw_ref[i,0] = det_yaw

##   PID姿态跟踪
    det_rhip_roll_inte = det_rhip_roll_inte + data0[i]-SensorValue[i][0]
    det_rhip_pitch_inte = det_rhip_pitch_inte + data1[i]-SensorValue[i][1]
    det_rknee_pitch_inte = det_rknee_pitch_inte + data2[i]-SensorValue[i][2]
    det_rankle_pitch_inte = det_rankle_pitch_inte + data3[i]-SensorValue[i][3]   
    det_rankle_roll_inte = det_rankle_roll_inte + data4[i]-SensorValue[i][4]
    det_lhip_roll_inte = det_lhip_roll_inte + data5[i]-SensorValue[i][5]
    det_lhip_pitch_inte = det_lhip_pitch_inte + data6[i]-SensorValue[i][6]    
    det_lknee_pitch_inte = det_lknee_pitch_inte + data7[i]-SensorValue[i][7]
    det_lankle_pitch_inte = det_lankle_pitch_inte + data8[i]-SensorValue[i][8]  
    det_lankle_roll_inte = det_lankle_roll_inte + data9[i]-SensorValue[i][9]
    yaw_err_intergral = yaw_err_intergral + (yaw_ref[i,0]-SensorValue[i][37])


    body_Rroll = 0.05*(data0[i]-(SensorValue[i][0]))+0.0001*((data0[i]-(SensorValue[i][0]))-(data0[i-1]-(SensorValue[i-1][0])))/(dt)+0.005*det_rhip_roll_inte
    body_Rpitch = 0.3*(data1[i]-(SensorValue[i][1]))+0.0008*((data1[i]-(SensorValue[i][1]))-(data1[i-1]-(SensorValue[i-1][1])))/(dt)+0.005*det_rhip_pitch_inte
##    det_rknee_pitch =0.45*(data2[i]-(SensorValue[i][2]))+0.002*((data2[i]-(SensorValue[i][2]))-(data2[i-1]-(SensorValue[i-1][2])))/(dt)+0.008*det_rknee_pitch_inte
##    det_rankle_pitch =0.2*(data3[i]-(SensorValue[i][3]))+0.002*((data3[i]-(SensorValue[i][3]))-(data3[i-1]-(SensorValue[i-1][3])))/(dt)+0.003*det_rankle_pitch_inte          
##    det_rankle_roll =0.1*(data4[i]-(SensorValue[i][4]))+0.004*((data4[i]-(SensorValue[i][4]))-(data4[i-1]-(SensorValue[i-1][4])))/(dt)+0.001*det_rankle_roll_inte

    body_Lroll = 0.2*(data5[i]-(SensorValue[i][5]))+0.0001*((data5[i]-(SensorValue[i][5]))-(data5[i-1]-(SensorValue[i-1][5])))/(dt)+0.005*det_lhip_roll_inte
    body_Lpitch = 0.3*(data6[i]-(SensorValue[i][6]))+0.0002*((data6[i]-(SensorValue[i][6]))-(data6[i-1]-(SensorValue[i-1][6])))/(dt)+0.005*det_lhip_pitch_inte
##    det_lknee_pitch =0.5*(data7[i]-(SensorValue[i][7]))+0.002*((data7[i]-(SensorValue[i][7]))-(data7[i-1]-(SensorValue[i-1][7])))/(dt)+0.005*det_lknee_pitch_inte
##    det_lankle_pitch =0.3*(data8[i]-(SensorValue[i][8]))+0.003*((data8[i]-(SensorValue[i][8]))-(data8[i-1]-(SensorValue[i-1][8])))/(dt)+0.005*det_lankle_pitch_inte         
##    det_lankle_roll =0.05*(data9[i]-(SensorValue[i][9]))+0.004*((data9[i]-(SensorValue[i][9]))-(data9[i-1]-(SensorValue[i-1][9])))/(dt)+0.002*det_lankle_roll_inte
    
##  奇怪的是命名是右脚调整，确影响到左腿的hip关节。这里进行调整
    det_rknee_pitch =0.45*(data2[i]-(SensorValue[i][2]))+0.002*((data2[i]-(SensorValue[i][2]))-(data2[i-1]-(SensorValue[i-1][2])))/(dt)+0.001*det_rknee_pitch_inte
    det_rankle_pitch =0.35*(data3[i]-(SensorValue[i][3]))+0.002*((data3[i]-(SensorValue[i][3]))-(data3[i-1]-(SensorValue[i-1][3])))/(dt)+0.001*det_rankle_pitch_inte          
    det_rankle_roll =0.3*(data4[i]-(SensorValue[i][4]))+0.002*((data4[i]-(SensorValue[i][4]))-(data4[i-1]-(SensorValue[i-1][4])))/(dt)+0.001*det_rankle_roll_inte

    det_lknee_pitch =0.4*(data7[i]-(SensorValue[i][7]))+0.001*((data7[i]-(SensorValue[i][7]))-(data7[i-1]-(SensorValue[i-1][7])))/(dt)+0.001*det_lknee_pitch_inte
    det_lankle_pitch =0.3*(data8[i]-(SensorValue[i][8]))+0.002*((data8[i]-(SensorValue[i][8]))-(data8[i-1]-(SensorValue[i-1][8])))/(dt)+0.002*det_lankle_pitch_inte         
    det_lankle_roll =0.05*(data9[i]-(SensorValue[i][9]))+0.004*((data9[i]-(SensorValue[i][9]))-(data9[i-1]-(SensorValue[i-1][9])))/(dt)+0.002*det_lankle_roll_inte    


##    yaw = 0.05*(yaw_ref[i,0]-SensorValue[i][37])+0.0005*(yaw_ref[i,0]-SensorValue[i][37]-(yaw_ref[i-1,0]-SensorValue[i-1][37]))/(dt)+0.005*yaw_err_intergral
    yaw = 0.05*(yaw_ref[i,0]-SensorValue[i][37])+0.0005*(yaw_ref[i,0]-SensorValue[i][37]-(yaw_ref[i-1,0]-SensorValue[i-1][37]))/(dt)+0.005*yaw_err_intergral
    
    return  body_Lroll, body_Lpitch, body_Rroll, body_Rpitch, det_rknee_pitch, det_rankle_pitch, det_rankle_roll, det_lknee_pitch, det_lankle_pitch, det_lankle_roll, yaw





def main(robotIP, PORT):
    global roll_err_intergral
    global pitch_err_intergral
    
    global det_rhip_roll_inte
    global det_rhip_pitch_inte    
    global det_rknee_pitch_inte
    global det_rankle_pitch_inte
    global det_rankle_roll_inte
    global det_lhip_roll_inte
    global det_lhip_pitch_inte    
    global det_lknee_pitch_inte
    global det_lankle_pitch_inte
    global det_lankle_roll_inte
    global data_angle0
    global data_angle1
    global data_angle2
    global data_angle3
    global data_angle4
    global data_angle5
    global data_angle6
    global data_angle7
    global data_angle8
    global data_angle9
    global angle

    pitch_err_intergral= 0
    roll_err_intergral = 0  
    det_rhip_roll_inte=0
    det_rhip_pitch_inte=0    
    det_rknee_pitch_inte=0
    det_rankle_pitch_inte=0
    det_rankle_roll_inte=0
    det_lhip_roll_inte=0
    det_lhip_pitch_inte=0   
    det_lknee_pitch_inte=0
    det_lankle_pitch_inte=0
    det_lankle_roll_inte=0

    motionProxy = ALProxy("ALMotion", robotIP, PORT)
    memoryProxy = ALProxy("ALMemory", robotIP, PORT)
    dcm = ALProxy("DCM",robotIP,PORT)
    

    dcm.createAlias([
    "legmove",[
    "RHipRoll/Position/Actuator/Value",
    "RHipPitch/Position/Actuator/Value",
    "RKneePitch/Position/Actuator/Value",
    "RAnklePitch/Position/Actuator/Value",
    "RAnkleRoll/Position/Actuator/Value",
    "LHipRoll/Position/Actuator/Value",
    "LHipPitch/Position/Actuator/Value",
    "LKneePitch/Position/Actuator/Value",
    "LAnklePitch/Position/Actuator/Value",
    "LAnkleRoll/Position/Actuator/Value",
    "LHipYawPitch/Position/Actuator/Value"
    ]
    ])
   

    
    SenorList = ["Device/SubDeviceList/RHipRoll/Position/Sensor/Value","Device/SubDeviceList/RHipPitch/Position/Sensor/Value","Device/SubDeviceList/RKneePitch/Position/Sensor/Value",
                 "Device/SubDeviceList/RAnklePitch/Position/Sensor/Value","Device/SubDeviceList/RAnkleRoll/Position/Sensor/Value","Device/SubDeviceList/LHipRoll/Position/Sensor/Value",
                 "Device/SubDeviceList/LHipPitch/Position/Sensor/Value","Device/SubDeviceList/LKneePitch/Position/Sensor/Value","Device/SubDeviceList/LAnklePitch/Position/Sensor/Value",
                 "Device/SubDeviceList/LAnkleRoll/Position/Sensor/Value",                                #十个关节
                 "Device/SubDeviceList/InertialSensor/AngleX/Sensor/Value","Device/SubDeviceList/InertialSensor/AngleY/Sensor/Value","Device/SubDeviceList/InertialSensor/AngleZ/Sensor/Value",  #姿态传感器角度
                 "Device/SubDeviceList/InertialSensor/GyroscopeX/Sensor/Value","Device/SubDeviceList/InertialSensor/GyroscopeY/Sensor/Value","Device/SubDeviceList/InertialSensor/GyroscopeZ/Sensor/Value",#角加速度
                 "Device/SubDeviceList/InertialSensor/AccelerometerX/Sensor/Value","Device/SubDeviceList/InertialSensor/AccelerometerY/Sensor/Value",
                 "Device/SubDeviceList/InertialSensor/AccelerometerZ/Sensor/Value",#角加速度
                 "Device/SubDeviceList/LFoot/FSR/FrontLeft/Sensor/Value",
                 "Device/SubDeviceList/LFoot/FSR/FrontRight/Sensor/Value",
                 "Device/SubDeviceList/LFoot/FSR/RearLeft/Sensor/Value",
                 "Device/SubDeviceList/LFoot/FSR/RearRight/Sensor/Value",
                 "Device/SubDeviceList/RFoot/FSR/FrontLeft/Sensor/Value",
                 "Device/SubDeviceList/RFoot/FSR/FrontRight/Sensor/Value",
                 "Device/SubDeviceList/RFoot/FSR/RearLeft/Sensor/Value",
                 "Device/SubDeviceList/RFoot/FSR/RearRight/Sensor/Value",#八个压力传感器
                 "Device/SubDeviceList/RHipRoll/ElectricCurrent/Sensor/Value","Device/SubDeviceList/RHipPitch/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/RKneePitch/ElectricCurrent/Sensor/Value","Device/SubDeviceList/RAnklePitch/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/RAnkleRoll/ElectricCurrent/Sensor/Value","Device/SubDeviceList/LHipRoll/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/LHipPitch/ElectricCurrent/Sensor/Value","Device/SubDeviceList/LKneePitch/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/LAnklePitch/ElectricCurrent/Sensor/Value","Device/SubDeviceList/LAnkleRoll/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/LHipYawPitch/Position/Sensor/Value"                 
                 ]
    FileLength = 3000
##  加入偏航角度  
    FileWidth = 38

    SensorValue = np.zeros([FileLength,FileWidth])
    out_nao = np.zeros([FileLength,426])
    mea = np.zeros(21)
    NNx = 40
    dt = 0.02

    
    #DCM
    start = time.clock()
    for i in range(1,FileLength):
##     读取数据处理，           
        meax = SensorValue[i-1]
        mea[0:3] = meax[10:13]
        mea[3:13] = meax[0:10]
        mea[13:17] = meax[23:27]
        mea[17:] = meax[19:23]
##     阈值触发，人工判断摔倒，存储文件
        if ((abs(mea[0]>0.5))or(abs(mea[1]>0.5))):
            print "I have fallen!"
            break
        else:
            stnum = i*dt
        ##%当前时刻对应的步行周期
        #### =========================================================   
            xtn = find0(n_period,stnum,0,40)
            xtn = xtn[0]        
    ##      角度实时计算
            out = nao_unknown_slope_previewcontrol(i,mea)
            out_nao[i,:] = out
            
            SensorValue[i] = memoryProxy.getListData(SenorList)
    ##      pd姿态控制       
            if i>= 2:
                [body_Lroll,body_Lpitch,body_Rroll,body_Rpitch,det_rknee_pitch,det_rankle_pitch,det_rankle_roll,
                det_lknee_pitch,det_lankle_pitch,det_lankle_roll,yaw] = vertical_body_joint_tracking(SensorValue,out_nao[:,0],out_nao[:,1],out_nao[:,2],out_nao[:,3],out_nao[:,4],
                                                                                                  out_nao[:,5],out_nao[:,6],out_nao[:,7],out_nao[:,8],out_nao[:,9],dt,i,xtn)
##                [body_Lroll,body_Lpitch,body_Rroll,body_Rpitch,det_rknee_pitch,det_rankle_pitch,det_rankle_roll,
##                 det_lknee_pitch,det_lankle_pitch,det_lankle_roll,yaw] = vertical_body_joint_tracking(SensorValue,data_angle0,data_angle1,data_angle2,data_angle3,data_angle4,
##                                                                                              data_angle5,data_angle6,data_angle7,data_angle8,data_angle9,dt,i,xtn)                
            else:
                body_Lroll = 0; body_Lpitch = 0; body_Rroll = 0; body_Rpitch = 0;det_rknee_pitch = 0; det_rankle_pitch = 0;
                det_rankle_roll = 0; det_lknee_pitch = 0; det_lankle_pitch = 0; det_lankle_roll = 0; yaw=0      
                                
            out_nao[i,0] = angle[i-1,0]+ body_Rroll
            out_nao[i,1] = angle[i-1,1]+ body_Rpitch
            out_nao[i,5] = angle[i-1,5]+ body_Lroll
            out_nao[i,6] = angle[i-1,6]+ body_Lpitch
            out_nao[i,2] = angle[i-1,2]+ det_rknee_pitch
            out_nao[i,3] = angle[i-1,3]+ det_rankle_pitch
            out_nao[i,4] = angle[i-1,4]+ det_rankle_roll
            out_nao[i,7] = angle[i-1,7]+ det_lknee_pitch
            out_nao[i,8] = angle[i-1,8]+ det_lankle_pitch
            out_nao[i,9] = angle[i-1,9]+ det_lankle_roll      
            
    ##======收发数据                    
            t1x = time.clock()-start
            t = dcm.getTime(1000*dt*i-1000*t1x) 
            dcm.setAlias(["legmove","Merge","time-mixed",[
            [[out_nao[i,0],t]],[[out_nao[i,1],t]],
            [[out_nao[i,2],t]],[[out_nao[i,3],t]],
            [[out_nao[i,4],t]],[[out_nao[i,5],t]],
            [[out_nao[i,6],t]],[[out_nao[i,7],t]],
            [[out_nao[i,8],t]],[[out_nao[i,9],t]],
            [[yaw,t]]
            ]])
                       
            
            ## 延时8ms
            if i< 600:
                time.sleep(0.006)
        
    print "t:", time.clock()-start
    
    #程序执行完毕写文件
    np.savetxt("DCM_NAO_flat_ground_online_with_posture_control_mode_adjust_pre_上坡2_sensor1.txt", SensorValue, fmt="%f")
    np.savetxt("DCM_NAO_flat_ground_online_with_posture_control_mode_adjust_pre_上坡2_output1.txt", out_nao, fmt="%f")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ip", type=str, default="169.254.228.115",
                        help="Robot ip address")
    parser.add_argument("--port", type=int, default=9559,
                        help="Robot port number")

    args = parser.parse_args()
    main(args.ip, args.port)










    
