/**
MPCClass.h

Description:	Header file of MPCClass

*/
#pragma once

#include "QP/QPBaseClass.h"
#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>
#include <Eigen/StdVector>
#include <vector> 


using namespace Eigen;
using namespace std;


/// constant variable defintion
const int _footstepsnumber = 12;        //  number of _footstepnumber
const double _dt = 0.1;                 //sampling time
const int _nh = 10;                     /// =PreviewT/_dt: number of sampling time for predictive window: <= 2*_nT; (_dt defined in MpcRTControlClass.h: dt_mpc)  	
const double _tstep = 0.6;              ///step period
const int _nT = 6;                      /// round(_tstep/_dt)  the number of one step cycle
const int _nstep = 2;                   /// maximal footstep locations where the predictive windwo covers
const int _Nt = 5*_nh + 3*_nstep;       /// _Nt = 5*_nh + 3*_nstep;  the number of the variable of the optimization problem
const int _Nt_half1 = 20;
const int _Nt_half2 = _Nt- _Nt_half1;
const int _nsum = (_footstepsnumber-1)*_nT; /// number of whole control loop


const double _height_offset = 0.05;     
const double _height_squat_time =2; 
 
 

class MPCClass : public QPBaseClass
{
public:
	MPCClass();
	virtual ~MPCClass() {};

	//void FootStepNumberInputs(int footstepsnumber);
	void FootStepInputs(double stepwidth, double steplength, double stepheight);
	
	void Initialize();
	
	void Matrix_large(int i_f);
	
	Eigen::Matrix<double,_Nt,_Nt> _matrix_large1,_matrix_large2;
	
	
	void CoM_foot_trajection_generation_local(int i, Eigen::Matrix<double,18,1> estimated_state, Eigen::Vector3d _Rfoot_location_feedback, Eigen::Vector3d _Lfoot_location_feedback,double lamda,bool _ISwalking);
	
	void solve_reactive_step_body_inclination_CoMz(); 
	void solve_reactive_step_body_inclination(); 
	void solve_reactive_step(); 	
	
	int Indexfind(double goalvari, int xyz);

	Eigen::MatrixXd Matrix_ps(Eigen::Matrix<double,3,3> a, int nh, Eigen::RowVector3d cxps);
	Eigen::MatrixXd Matrix_pu(Eigen::Matrix<double,3,3> a, Eigen::Matrix<double,3,1> b, int nh, Eigen::RowVector3d cxpu);		
	void Solve();

	void Foot_trajectory_solve(int j_index, bool _ISwalking);
	Vector3d X_CoM_position_squat(int walktime, double dt_sample);

	// current state based on the past one and two actual sampling time;
	Vector3d XGetSolution_CoM_position(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);
	Vector3d XGetSolution_Foot_positionR(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);
	Vector3d XGetSolution_Foot_positionL(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);	
	Vector3d XGetSolution_body_inclination(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);	
	
	int Get_maximal_number(double dtx);
	
	int Get_maximal_number_reference();


	int is_sigular(int num);
	
	void File_wl();	
	
	std::string _robot_name;
	double _robot_mass;
	double _lift_height;
	
	double _GGG;
	double _HCOM;
	double _half_hip_width;	
	double _foot_width;
	
	
	int _method_flag;

	int _n_end_walking;

	//////////////////////////////// substitute the 
	int _j_period;

	Eigen::Matrix<double,_footstepsnumber,1> _tx;	
	Eigen::Matrix<double,_footstepsnumber,1> _ts, _td;


	
////////// generate the COM trajectory during the DSP	
	int _nTdx;
	
	
	int _bjx1, _bjx2, _mx;	
	
	int _j_count;
	
        Vector3d _F_R, _F_L,_M_R,_M_L;
	Vector3d _ZMPxy_realx;
	
	Vector3d _comxyzx,_comvxyzx,_comaxyzx,_thetaxyx,_thetavxyx,_thetaaxyx,_Lfootxyzx,_Rfootxyzx;	
		
	void Zmp_distributor(int walktime, double dt_sample);
	
protected: 


private:    
        //parameters declaration  	
        Eigen::Matrix<double,_footstepsnumber,1> _steplength, _stepwidth, _stepheight,_lift_height_ref;	
        Eigen::Matrix<double,_footstepsnumber,1> _footx_ref, _footy_ref, _footz_ref;	
/*	Eigen::Matrix<double,_footstepsnumber,1> _ts, _td;*/	
	 
	Eigen::Matrix<double,_nsum,1> _t;
	
	Eigen::Matrix<double,1,_nsum> _zmpx_real, _zmpy_real;
	Eigen::Matrix<double,1,_nsum> _comx, _comvx, _comax;
	Eigen::Matrix<double,1,_nsum> _comy, _comvy, _comay;
	Eigen::Matrix<double,1,_nsum> _comz, _comvz, _comaz;	
	Eigen::Matrix<double,1,_nsum> _thetax, _thetavx, _thetaax;
	Eigen::Matrix<double,1,_nsum> _thetay, _thetavy, _thetaay;	
	Eigen::Matrix<double,1,_nsum> _thetaz, _thetavz, _thetaaz;	
	Eigen::Matrix<double,1,_nsum> _torquex_real, _torquey_real;


	// initial parameters for MPC
	double _hcom;
	Eigen::Matrix<double,_nh, 1> _Hcom;
	Eigen::Matrix<double,1,1> _ggg;	
	
		
	// CoM+angular momentum state and contro input
	Eigen::Matrix<double,3,_nsum> _xk,_yk,_zk,_thetaxk,_thetayk;
	Eigen::Matrix<double,1,_nsum> _x_vacc_k,_y_vacc_k,_z_vacc_k,_thetax_vacc_k,_thetay_vacc_k;
	
	
	Eigen::Matrix<double,3,3> _a;
	Eigen::Matrix<double,3,1> _b;
	Eigen::RowVector3d _cp,_cv,_ca;
		
	//predictive model
	Eigen::Matrix<double,_nh,3> _pps,_pvs,_pas;	
	Eigen::Matrix<double,_nh,_nh> _ppu,_pvu,_pau;		
	Eigen::Matrix<double,_nh,_nh> _ppu_2, _pvu_2;
	
	int xyz1;  //flag for find function 
	int xyz2;
		
	
	//footz reference: matrix operations
	Eigen::Matrix<double,_nsum,1> _Zsc;	
	
	Eigen::Matrix<double, _nh,1> _v_i; 
	Eigen::Matrix<double, _nh,_nstep> _VV_i;	
	
	
	
	//genearated foot location
	Eigen::Matrix<double, _footstepsnumber,1> _footx_real, _footy_real, _footz_real;
	Eigen::Matrix<double, _nsum,1> _footx_real_next, _footy_real_next, _footz_real_next;
	Eigen::Matrix<double, _nsum,1> _footx_real_next1, _footy_real_next1, _footz_real_next1;
	Eigen::Matrix<double, 3,_footstepsnumber> _footxyz_real;
	
	Eigen::Matrix<double, 1,_nsum> _Lfootx, _Lfooty,_Lfootz, _Lfootvx, _Lfootvy,_Lfootvz, _Lfootax, _Lfootay,_Lfootaz;
	Eigen::Matrix<double, 1,_nsum> _Rfootx, _Rfooty,_Rfootz, _Rfootvx, _Rfootvy,_Rfootvz, _Rfootax, _Rfootay,_Rfootaz;
        double _ry_left_right;
	
	
	//vertical height constraints
	double _z_max, _z_min, _detzb;	
	
	// foot step and width constraints
	double _footx_max, _footx_min,_footy_max,_footy_min, _detfootx,_detfooty;
	
	double _mass,  _rad,  _j_ini;
	
	// zmp-constraints
	double _zmpx_ub,_zmpx_lb, _zmpy_ub, _zmpy_lb, _detzmppxb, _detzmppyb;

	// angle range
	double _thetax_max,  _thetax_min,  _thetay_max, _thetay_min, _detthetax,_detthetay;	
	
	// torque range
	double _torquex_max, _torquex_min, _torquey_max, _torquey_min, _dettorquex,_dettorquey; 
	
	// swing foot velocity constraints
	double _footx_vmax, _footx_vmin,_footy_vmax,_footy_vmin,_detfootvx,_detfootvy;


        double _fx, _fy;
	double _fxx_global, _fyy_global;	
	
	// solution preparation
// 	int _Nt;
	Eigen::Matrix<double, _Nt,1> _V_ini;
	Eigen::Vector2d _V_inix;			
	Eigen::Matrix<double,_Nt, _nsum> _V_optimal;
	
	Eigen::Matrix<double, _nsum,1> _flag, _flag_global;
	
        Eigen::Matrix<double, _nstep,1> _Lx_ref, _Ly_ref,_Lz_ref;
	Eigen::Matrix<double, _nh,1> _comx_center_ref, _comy_center_ref,_comz_center_ref, _thetax_center_ref, _thetay_center_ref; 	
	
	
	// weight coefficients
	double _Rx,     _Ry,     _Rz;
	double _alphax, _alphay, _alphaz;
	double _beltax, _beltay, _beltaz;
	double _gamax,  _gamay,  _gamaz;
	double _Rthetax, _Rthetay;
	double _alphathetax, _alphathetay;
	double _beltathetax, _beltathetay;
	
	// time cost consumption
	Eigen::RowVectorXd _tcpu;
// predictive model control_tracking with time_varying height		
	int _loop;
	
	Eigen::Matrix<double,_nh,_nh> A_unit;
	Eigen::Matrix<double,_nstep,_nstep> C_unit;		

	int _bjxx_1, _bjxx;
	Eigen::Matrix<double,_nh,1> _t_f;
	
//	int _bjx1, _bjx2, _mx;
	Eigen::Vector4d  _tnx;
	
	int _n_vis, xxx, xxx1,xxx2;

	// sqp model	
	Eigen::Matrix<double,_nh,_nh> _WX, _WY, _WZ, _WthetaX, _WthetaY;
	Eigen::Matrix<double,_nstep,_nstep> _PHIX, _PHIY, _PHIZ;
	Eigen::Matrix<double,_Nt,_Nt> _Q_goal, _Q_goal1;
	Eigen::Matrix<double,_Nt,1> _q_goal, _q_goal1;	

	Eigen::Matrix<double,_nh,_Nt> _Sjx,_Sjy,_Sjz,_Sjthetax,_Sjthetay;
	Eigen::Matrix<double,_nstep,_Nt>   _Sfx,_Sfy,_Sfz;
	
	// inequality constraints	
	Eigen::Matrix<double,_nh,_Nt> _H_q_upx,_H_q_lowx,_H_q_upy,_H_q_lowy;
	Eigen::Matrix<double,_nh,1> _F_zmp_upx,_F_zmp_lowx,_F_zmp_upy,_F_zmp_lowy;
	
	// zmp constraints
	Eigen::Matrix<double,1,_nh> _Si;
// 	Eigen::Matrix<double,_Nt,_Nt> _phi_i_x_up,_phi_i_x_low,_phi_i_y_up,_phi_i_y_low;
	Eigen::Matrix<double,1,_Nt> _phi_i_x_up,_phi_i_y_up;
// 	Eigen::Matrix<double,_Nt,_nh> _p_i_x_t_up,_p_i_x_t_low,_p_i_y_t_up,_p_i_y_t_low;	
	Eigen::Matrix<double,_nh,_Nt> _p_i_x_t_up,_p_i_x_t_low,_p_i_y_t_up,_p_i_y_t_low;	
	Eigen::Matrix<double,1,_nh> _del_i_x_up,_del_i_x_low,_del_i_y_up,_del_i_y_low, _det_del_i_x,_det_del_i_y;
	
	Eigen::Matrix<double,_nh,_Nt> _q_upx,_q_lowx,_q_upy,_q_lowy;
	Eigen::Matrix<double,_nh,1> _qq_upx,_qq_lowx,_qq_upy,_qq_lowy,_qq1_upx,_qq1_lowx,_qq1_upy,_qq1_lowy;

	
	Eigen::Matrix<double,_nh,_Nt> _t_upx,_t_lowx,_t_upy,_t_lowy;
	Eigen::Matrix<double,_nh,1> _tt_upx,_tt_lowx,_tt_upy,_tt_lowy,_tt1_upx,_tt1_lowx,_tt1_upy,_tt1_lowy;
	
	Eigen::Matrix<double,_nh,_Nt> _H_h_upz,_H_h_lowz;
	Eigen::Matrix<double,_nh,1> _F_h_upz,_F_h_lowz, _delta_footz_up, _delta_footz_low;
	
	Eigen::Matrix<double,_nh,_Nt> _H_hacc_lowz;
	Eigen::Matrix<double,_nh,1> _F_hacc_lowz, _delta_footzacc_up;
	

	Eigen::Matrix<double,1,_Nt> _Footvx_max,_Footvx_min,_Footvy_max,_Footvy_min;
	Eigen::Matrix<double,1,1> _footubxv,_footlbxv,_footubyv,_footlbyv;
	
	Eigen::Matrix<double,_nstep,_Nt> _H_q_footx_up,_H_q_footx_low,_H_q_footy_up,_H_q_footy_low;
	Eigen::Matrix<double,_nstep,1> _F_foot_upx,_F_foot_lowx,_F_foot_upy,_F_foot_lowy;

	
	// equality constraints
	Eigen::Matrix<double, 1, _Nt> _H_q_footz;
	Eigen::Matrix<double, 1, 1> _F_footz;
	
	Eigen::Matrix<double,_nh, _Nt> _h_h;
	Eigen::Matrix<double,_nh, 1> _hhhx;
	
	Eigen::Matrix<double,_nh, _Nt> _a_hx,  _a_hy;
	Eigen::Matrix<double,_nh, 1> _a_hxx, _a_hyy;

	
	//foot location constraints
	Eigen::RowVector2d _Sfoot;

// 	int xxx_vector=30;
	
	
	/// for c++11
// 	vector <Eigen::Matrix<double,_Nt, _Nt> > ZMPx_constraints_offfline;
	//vector <Eigen::Matrix<double,_Nt, _Nt>> ZMPy_constraints_offfline;
	//vector <Eigen::Matrix<double,_Nt, _Nt>> _phi_i_x_up_est, _phi_i_y_up_est;	

	//vector <Eigen::Matrix<double,_nh, _Nt>> ZMPx_constraints_half;
	//vector <Eigen::Matrix<double,_nh, _Nt>> ZMPy_constraints_half;
	//vector <Eigen::Matrix<double,3, _Nt>> _xkZMPx_constraints,_ykZMPy_constraints;
	//vector <Eigen::Matrix<double,3, _Nt>> _zkZMPx_constraints,_zkZMPy_constraints;
    
	//// for c++98
	//vector <Eigen::Matrix<double,_Nt,_Nt>> ZMPy_constraints_offfline;
	std::vector<Eigen::Matrix<double,_Nt,_Nt>,  Eigen::aligned_allocator<Eigen::Matrix<double,_Nt,_Nt>>> ZMPx_constraints_offfline,ZMPy_constraints_offfline;
	std::vector<Eigen::Matrix<double,_Nt,_Nt>,  Eigen::aligned_allocator<Eigen::Matrix<double,_Nt,_Nt>>> _phi_i_x_up_est, _phi_i_y_up_est;	

	std::vector <Eigen::Matrix<double,_nh, _Nt>, Eigen::aligned_allocator<Eigen::Matrix<double,_nh,_Nt>>> ZMPx_constraints_half, ZMPy_constraints_half;
	std::vector <Eigen::Matrix<double,3,_Nt>, Eigen::aligned_allocator<Eigen::Matrix<double,3,_Nt>>> _xkZMPx_constraints,_ykZMPy_constraints;
	std::vector <Eigen::Matrix<double,3,_Nt>, Eigen::aligned_allocator<Eigen::Matrix<double,3,_Nt>>> _zkZMPx_constraints,_zkZMPy_constraints;


        Eigen::Matrix<double,_Nt, _Nt> _x_offline1_va;
        Eigen::Matrix<double,_nh, _Nt> _x_offline2_va;
        Eigen::Matrix<double,3, _Nt> _x_offline3_va;
        Eigen::Matrix<double,3, 3> _x_offline4_va;    
        
        
        Eigen::Matrix<double,1, _Nt> ZMPx_constraints_halfxxx1,ZMPx_constraints_halfxxx2;
        

        
        Eigen::Matrix<double,_Nt,1> ZMPx_constraints_halfyyy1,ZMPx_constraints_halfyyy2,ZMPx_constraints_halfyyy3,ZMPx_constraints_halfyyy4;
	
	
        Eigen::Matrix<double,1, _Nt> _x_offline4_vax;       
        Eigen::Matrix<double,_Nt,1> _x_offline4_vay;        
        
        
        
        
	
	double _Footx_global_relative,_Footy_global_relative;
	
	
	Eigen::Matrix<double,_Nt, _nh> _ZMPx_constraints_half2, _ZMPy_constraints_half2;
	Eigen::Matrix<double,_Nt, _Nt> _phi_i_x_up1, _phi_i_y_up1;
        
        Eigen::Matrix<double,_nh, _Nt> _ZMPx_constraints_half_va,_ZMPy_constraints_half_va;
        
        
	
	Eigen::Matrix<double,28,_nsum> CoMMM_ZMP_foot;
	
// 	Eigen::Matrix4d AAAaaa1;
	double _ZMP_ratio;
	double _t_end_walking;
	
	Eigen::Vector4d  _comy_matrix_inv;
	
	
	
	
///    offline calculated matrices	
	Eigen::Matrix<double, _nh, _nh> _ppu_T;
	Eigen::Matrix<double, _nh, 3> _pvupvs, _ppupps;
	
	
	Eigen::Matrix<double, _nh, _Nt> _ppuSjx, _ppuSjy;
	Eigen::Matrix<double, _nh, _Nt> _pauSjx, _pauSjy;	
	Eigen::Matrix<double, _nh, _Nt> _pauSjz1,_pauSjz2,_pauSjz11,_pauSjz21;	
	Eigen::Matrix<double, _nh, _Nt> _pauSjthetay,_pauSjthetax;
	
	vector <Eigen::Matrix3d> _xkzk_constraints;	
	Eigen::Matrix<double, 1, 1> _gzmpxub,_gzmpxlb,_gzmpyub,_gzmpylb;


////==========================================================================	
	///////////////for polynomial intepolation
	Eigen::Matrix<double, 4,4>  _AAA_inv;
	
		//strange error !!!!!!!!!!!!!!!!!!!!be careful to modified the function: maybe caused by c++98 bug
//	void solve_AAA_inv(Eigen::Matrix<double, 4, 1> t_plana); 
	void solve_AAA_inv(const Eigen::Matrix4d& t_plana); 

	Eigen::Matrix<double, 7, 7> solve_AAA_inv_x(Eigen::Vector3d t_plan); 
	
	
	////for lowel level control: ZMP optimal distribution
// 	Vector3d _F_R, _F_L,_M_R,_M_L;
	Matrix<double,3,3> _Co_L, _Co_R;
		 
	void Force_torque_calculate(Vector3d comxyzx1,Vector3d comaxyzx1,Vector3d thetaaxyx1,Vector3d Lfootxyz,Vector3d Rfootxyz);

	void zmp_interpolation(int t_int,int walktime, double dt_sample);
	
// 	void Zmp_distributor(int walktime, double dt_sample);	
	
	

};
