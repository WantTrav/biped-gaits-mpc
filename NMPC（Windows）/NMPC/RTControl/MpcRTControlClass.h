/*****************************************************************************
MpcRTControlClass.h
*****************************************************************************/
#pragma once

#include "MPC\MPCClass.h"
#include "RobotState\RobotStateClass.h"
#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>

const double dt_mpc = 0.05;   // definition of sampling time of MPC solover
const double _height_offsetx = 0.05;  
const double _height_offset_time = 2; 


using namespace Eigen;
using namespace std;



class MpcRTControlClass {
public:
	MpcRTControlClass();
	~MpcRTControlClass();
/////////////////first layer
//	NLPClass nlp;
	
	/////////////second layer
	
	MPCClass mpc;
	
	RobotStateClass robotstate;
	
	int _walkdtime_max, _wal_max;	
	int _walkdtime1;
	
	Eigen::Matrix<double,27,1> WalkingReactStepping(int walkdtime); 
	
//	int File_wl(Eigen::MatrixXd Lfootz);

	///
	double GGG,TOTAL_MASS,HALF_HIP_WIDTH,LIFT_HEIGHT,HCOM;

	double sx_input, sy_input, sz_input;

	

// 	XBot::MatLogger::Ptr xbot_logger;
// 	void initLogger(int size, int interleave);
// 	void addToLog();
// 	int logger_len;
// 	int interlaver;
	
	void savedata();
	
	Eigen::MatrixXd _data_saving;
	
	
	

private:
    bool IsStartWalk;
  
  
	void StartWalking();
	void StopWalking();
	
	/// step parameters reference
	double stepwidthinput,steplengthinput,stepheightinput;

        int _refer_t_max;
	int _t_int;
	double _dtx;
	

	
	int _t_walkdtime_flag, _t_walkdtime_restart_flag;
	bool _stop_walking, _start_walking_again;
	
	double _ppx, _pix, _pdx,  _ppy,_piy, _pdy, _ppz,_piz,  _pdz, _ppthetax, _pdthetax, _pithetax,  _ppthetay, _pithetay,  _pdthetay,  _ppthetaz, _pithetaz,  _pdthetaz; 
	
	Eigen::VectorXd _error_com_position, _error_torso_angle;
	
	
	Eigen::VectorXd _flag_walkdtime;
	
	Eigen::VectorXd _stop_flag_walkdtime;
  
	

	
	Eigen::MatrixXd _COM_IN, _COM_est;
	Eigen::MatrixXd _body_IN;	
	Eigen::MatrixXd _FootR_IN;	
	Eigen::MatrixXd _FootL_IN;	

	Eigen::Matrix<double,18,1> _estimated_state;	
	Eigen::Vector3d _Rfoot_location_feedback,_Lfoot_location_feedback;
	
	

	

	
////////////////// define in this class;	
// 	time flag for online optimization
	
//         int walkdtime;	
	
	
	Eigen::Vector3d PelvisPos,body_thetax,LeftFootPosx,RightFootPosx;
	
	Eigen::Vector3d ZMPxy_realx,zmp_ref,F_R,F_L,M_L, M_R;
	
    int j_count,bjx1;
	double tx, td;

	

	


protected:
	int t_walkdtime_flag;
	int dt_sample;
	
	
	double _feedback_lamda;
	
	
// 	/****KMP based trajectory***********************/
// 	Vector6d _kmp_leg_traje;
	
	
	Eigen::Vector3d _F_r_mpc, _F_l_mpc,_M_r_mpc,_M_l_mpc;	
	
	
};

