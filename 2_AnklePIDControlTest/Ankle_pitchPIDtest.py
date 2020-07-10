# -*- encoding: UTF-8 -*-


import time
import argparse
import string
import numpy as np
from naoqi import ALProxy

def main(robotIP, PORT):
    motionProxy = ALProxy("ALMotion", robotIP, PORT)
    memoryProxy = ALProxy("ALMemory", robotIP, PORT)
    #postureProxy = ALProxy("ALRobotPosture", robotIP, PORT)
    dcm = ALProxy("DCM",robotIP,PORT)


    
    motionProxy.wakeUp()
  
    #read file
    
    Input = open(r'C:\Users\JC_ZHOU\Desktop\NAO\MpcCs50310-600.txt','r')
    data = []
    for i in range(0,10):
        a = Input.readline()
        a = string.split(a)
        data.append(a)
    data = map(map,[float,float,float,float,float,float,float,float,float,float], data)
    data0=data[0]
    data1=data[1]
    data2=data[2]
    data3=data[3]
    data4=data[4]
    data5=data[5]
    data6=data[6]
    data7=data[7]
    data8=data[8]
    data9=data[9]
    
    dcm.createAlias([
    "legmove",[    
    "RAnklePitch/Position/Actuator/Value",
    "LAnklePitch/Position/Actuator/Value"
    ]
    ])
    
    FileLength = 1500

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
                 "Device/SubDeviceList/RHipRoll/ElectricCurrent/Sensor/Value","Device/SubDeviceList/RHipPitch/ElectricCurrent/Sensor/Value","Device/SubDeviceList/RKneePitch/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/RAnklePitch/ElectricCurrent/Sensor/Value","Device/SubDeviceList/RAnkleRoll/ElectricCurrent/Sensor/Value","Device/SubDeviceList/LHipRoll/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/LHipPitch/ElectricCurrent/Sensor/Value","Device/SubDeviceList/LKneePitch/ElectricCurrent/Sensor/Value","Device/SubDeviceList/LAnklePitch/ElectricCurrent/Sensor/Value",
                 "Device/SubDeviceList/LAnkleRoll/ElectricCurrent/Sensor/Value",     #十个关节电流
                 ]
    FileWidth = len(SenorList)
    FileLength = 800
    data0=np.zeros(FileLength)
    data4=np.zeros(FileLength)
    data5=np.zeros(FileLength)
    data9=np.zeros(FileLength)
    
    SensorValue = np.zeros([FileLength,FileWidth])

    #DCM
    dt = 0.02
    start = time.clock()
    AngleX = np.zeros(FileLength)
    AngleY = np.zeros(FileLength)
    tglobal = np.zeros(FileLength)

    #define
    pitch_err_intergral = 0
    det_pitch = 0
    
    for i in range(0,FileLength):
        SensorValue[i] = memoryProxy.getListData(SenorList)   ##读取传感器数据
        
        AngleX[i] = SensorValue[i][10]
        AngleY[i] = SensorValue[i][11]

        if i>1:
            RAnklePicth = SensorValue[i][4]
            LAnklePicth = SensorValue[i][9]
            
           ## body_Rroll = 0.05*(data0[i]-(SensorValue[i][0]))+0.0001*((data0[i]-(SensorValue[i][0]))-(data0[i-1]-(SensorValue[i-1][0])))/(dt)+0.005*det_rhip_roll_inte

            pitch_err_intergral = pitch_err_intergral + (AngleY[i]+0.01)  ##积分
            det_pitch = 0.05*(AngleY[i]+0.01)+0.0001*(AngleY[i]-0.028-(AngleY[i-1]-0.028))/dt + 0.02*pitch_err_intergral
            
            t1x = time.clock()-start
            t = dcm.getTime(1000*dt*i-1000*t1x)
            tglobal[i] = t
            RAnklePicth0 = RAnklePicth-det_pitch
            LAnklePicth0 = LAnklePicth-det_pitch
            
            dcm.setAlias(["legmove","Merge","time-mixed",[
            [[det_pitch,t]],[[det_pitch,t]]
            ]])


        time.sleep(0.02)   ##延时16ms  
        

    #写文件
    np.savetxt("SensorTest.txt", SensorValue, fmt="%f")
    np.savetxt("angleX.txt",AngleX,fmt = "%f")
    np.savetxt("angleY.txt",AngleY,fmt = "%f")
    np.savetxt("tgoal.txt",tglobal,fmt = "%f")

    print  start
    print  time.clock()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ip", type=str, default="169.254.228.115",
                        help="Robot ip address")
    parser.add_argument("--port", type=int, default=9559,
                        help="Robot port number")

    args = parser.parse_args()
    main(args.ip, args.port)
