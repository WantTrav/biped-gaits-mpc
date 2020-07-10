/*****************************************************************************
RobotStateClass.h

*****************************************************************************/
#pragma once



#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>
// #include "RobotModel/RobotModelClass.h"

class RobotStateClass {
public:
	RobotStateClass();
	virtual ~RobotStateClass() {};

	void init();

	Eigen::Vector3d com;  //!< overall COM
	Eigen::Vector3d dcom; //!< COM velocity
	Eigen::Vector3d ddcom;//!< COM acceleration
	Eigen::Vector3d com_old;
	Eigen::Vector3d dcom_old;
	Eigen::Vector3d ddcom_old;
	Eigen::Vector3d Lcom; //!< overall angular momentum around COM
	Eigen::Vector3d dLcom;
	Eigen::Vector3d Lcom_old;
	Eigen::Vector3d dLcom_old;
	Eigen::Vector3d Lft;
	Eigen::Vector3d Rft;

	Eigen::Vector3d htot; //!< overall angular momentum around COM caculated by RBDL

	
	
	double Fzl,Fzr;


	Eigen::Vector3d fOrigin; //!< origin in support foot coordinate's frame
	Eigen::Vector3d ft_com; //!< com in support foot coordinate's frame
	Eigen::Vector3d ft_com_old;
	Eigen::Vector3d ft_dcom;
	Eigen::Vector3d ft_ddcom;
	
	Eigen::Vector3d gOrigin; //!< fOorigin in global frame
	Eigen::Vector3d gcom; //!< com in global frame
	Eigen::Vector3d gcom_old;
	Eigen::Vector3d gdcom;
	Eigen::Vector3d gdcom_old;
	Eigen::Vector3d gddcom;
	
	Eigen::Vector3d ghip;

	Eigen::Vector3d glft_old;
	Eigen::Vector3d grft_old;
	Eigen::Vector3d glft;
	Eigen::Vector3d grft;	
	
	Eigen::Matrix3d IMU_abs;	
	Eigen::Vector3d IMU_Euler;	
	Eigen::Vector3d IMU_LinearVel;
	Eigen::Vector3d IMU_AngularVel;
	Eigen::Vector3d IMU_LinearAcc;
	Eigen::Vector3d IMU_AngularAcc;	



private:


};


