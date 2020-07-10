/*****************************************************************************
RobotStateClass.cpp
*****************************************************************************/

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>
#include "RobotStateClass.h"

using namespace Eigen;
using namespace std;

RobotStateClass::RobotStateClass()
{

	com.setZero(); // overall COM;initialization of vectors
	dcom.setZero();
	ddcom.setZero();
	com_old.setZero();
	dcom_old.setZero();
	ddcom_old.setZero();
	Lcom .setZero();	// in pelvis local frame
	dLcom.setZero();
	Lcom_old.setZero();
	dLcom_old.setZero();
	Lft.setZero();	// in pelvis local frame
	Rft.setZero();	// in pelvis local frame

	
	htot.setZero();
	
	Fzl=0;
	Fzr=0;

	fOrigin.setZero();
	ft_com.setZero();		// COM in fOrigin frame with eye(3) orientation align with pelvis frame
	ft_com_old.setZero();
	ft_dcom.setZero();
	ft_ddcom.setZero();
	
	
	gOrigin.setZero();
	gcom.setZero();		// COM in world frame
	gcom_old = gcom;
	gdcom.setZero();
	gdcom_old.setZero();
	gddcom.setZero();
	
	ghip << 0.0, 0.0, 0.89;
	glft << 0, 0.103 ,0;		// left foot center in global/ world coordinate
	grft << 0, -0.103 ,0;		// right foot center in global/ world coordinate
	glft_old << 0, 0.103 ,0;
	grft_old << 0, -0.103 ,0;



	IMU_abs.setZero();
	IMU_abs(0,0)=1;
	IMU_abs(1,1)=1;
	IMU_abs(2,2)=1;
	
	
	IMU_Euler.setZero();	
	IMU_LinearVel.setZero();
	IMU_AngularVel.setZero();
	IMU_LinearAcc.setZero();
	IMU_AngularAcc.setZero();
	std::cout << "\n\n\n========= Constructing RobotStateClass finish=================" << std::endl;	
}


void RobotStateClass::init()
{
	gcom.setConstant(0.1);
	com.setRandom();
}