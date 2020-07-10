/*****************************************************************************
MPCClass.cpp

Description:    source file of MPCClass
*****************************************************************************/
#include "QP\QPBaseClass.h"        //////////////remember to include the h file in this header file
#include "MPCClass.h"
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>
//#include <boost/math/special_functions/round.hpp> //c++98 round function


using namespace Eigen;
using namespace std;


MPCClass::MPCClass()                    ///declaration function
	: QPBaseClass()
	, _robot_name("cogimon")
	, _robot_mass(69)
	, _lift_height(0.06)
	, _GGG(9.8)
	,_HCOM(0.89)
	,_half_hip_width(0.103)	
	,_foot_width(0.1)
	, _method_flag(0)
	,_n_end_walking(1)
	,_j_period(0)
	,_F_R(0,0,0)
	,_F_L(0,0,0)
	,_M_R(0,0,0)
	,_M_L(0,0,0)		
{
  cout<<"MPCClass object!!!!!!!!!!"<<endl;
}

//////////////////step parameters input============================================================
void MPCClass::FootStepInputs(double stepwidth, double steplength, double stepheight)
{	
	_steplength.setConstant(steplength);
	_steplength(0) = 0;
	/// for yanlong huang: icra2020
 //	_steplength(1) = steplength/2;
// 	_steplength(2) = steplength/2;
	
	_stepwidth.setConstant(stepwidth);
	_stepwidth(0) = stepwidth/2;
	
	_stepheight.setConstant(stepheight);     
    _steplength(_footstepsnumber-1) = 0;
    _steplength(_footstepsnumber-2) = 0;
    _steplength(_footstepsnumber-3) = 0;
    _steplength(_footstepsnumber-4) = 0;	
    _steplength(_footstepsnumber-5) = stepwidth/2;

	_lift_height_ref.setConstant(_lift_height);
	_lift_height_ref(0) = 0.00;
	_lift_height_ref(1) = 0.04;
    _lift_height_ref(_footstepsnumber-1) = 0;
    _lift_height_ref(_footstepsnumber-2) = 0.04;	
    _lift_height_ref(_footstepsnumber-3) = 0.04; 
    _lift_height_ref(_footstepsnumber-4) = 0.04; 
	
	cout<<"steplength definition"<<_steplength<<"\n"<<endl;

}

///////////////////////// initialize all the variables============================================================
void MPCClass::Initialize()
{ // //    variance step length
	_steplength(5) = 0.15;
	_steplength(6) = 0.05;
	_steplength(7) = 0;
	_steplength(8) = -0.05;	
	_steplength(9) = -0.05;
	_steplength(10) = 0;
	_steplength(11) = 0.05;	
	_steplength(12) = 0.15;	

//// desired motion///
        std::cout << "MPC_initialization00X...\n";
 	// ==step loctions setup==
    _footx_ref.setZero();
	_footy_ref.setZero();
	_footz_ref.setZero();	
  	for (int i = 1; i < _footstepsnumber; i++) {
 	  _footx_ref(i) = _footx_ref(i-1) + _steplength(i-1);
	  _footy_ref(i) = _footy_ref(i-1) + ((int)pow(-1.0,i-1))*_stepwidth(i-1);   
	  _footz_ref(i) = _footz_ref(i-1) + _stepheight(i-1);
	}	
        //cout <<"_footy_ref:"<<"\n"<<_footy_ref<<endl;

//       sampling time & step cycle	
	_ts.setConstant(_tstep);
	_td = 0.2*_ts;                                // dsp time
	
// 	_tx = Eigen::VectorXd::Zero(_footstepsnumber);   // start time for each step cycle
	_tx.setZero();   // start time for each step cycle	
  	for (int i = 1; i < _footstepsnumber; i++) {
 	  _tx(i) = _tx(i-1) + _ts(i-1);
	  _tx(i) = ((int)boost::math::iround(_tx(i)/_dt))*_dt -0.00001;	  
	}	
 //   cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(0.1/_dt))<<endl;
	//cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(0.149/_dt))<<endl;
	//cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(0.151/_dt))<<endl;
	//cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(0.2/_dt))<<endl;
 //   cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(-0.1/_dt))<<endl;
	//cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(-0.149/_dt))<<endl;
	//cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(-0.151/_dt))<<endl;
	//cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(-0.2/_dt))<<endl;
	//cout<< "boost round function test"<<"\n"<<((int)boost::math::iround(-0.05/_dt))<<endl;

	_t.setLinSpaced(_nsum,_dt,_tx(_footstepsnumber-1));  ///sampling time sequence for entire step period	
        // ==initial parameters & matrix for MPC==
        _hcom = _HCOM-_height_offset;                                                    /// intial comz
	_Hcom.setConstant(_hcom);                             /// for comz reference during the predictive window: matrix operations: height variance between support feet and body center	
	_ggg.setConstant(_GGG);                                            //// acceleration of gravity: 1*1 matrix: matrix operations

//////////////==============================state variable============================================
	_zmpx_real.setZero(); _zmpy_real.setZero();                              
	_comx.setZero(); _comvx.setZero(); _comax.setZero();
	_comy.setZero(); _comvy.setZero(); _comay.setZero();
	_comz.setZero(); _comvz.setZero(); _comaz.setZero();	
	_thetax.setZero(); _thetavx.setZero(); _thetaax.setZero();
	_thetay.setZero(); _thetavy.setZero(); _thetaay.setZero();
	_thetaz.setZero(); _thetavz.setZero(); _thetaaz.setZero();	
	_torquex_real.setZero(); _torquey_real.setZero();  
          //// optimized CoM postion, foot trajectory, torse angle, footstep location;
	_xk.setZero(); 
	_yk.setZero(); 
	_zk.setZero();
	_thetaxk.setZero(); 
	_thetayk.setZero();
	_x_vacc_k.setZero(); 
	_y_vacc_k.setZero();
	_z_vacc_k.setZero(); 
	_thetax_vacc_k.setZero(); 
	_thetay_vacc_k.setZero(); 	
///////================================================================================================			
	_a << 1, _dt, pow(_dt,2)/2,    
	      0,   1,            _dt,
	      0,   0,              1;
	_b << pow(_dt,3)/6,    
	      pow(_dt,2)/2,
	               _dt;		
	
	_cp.setZero();
	_cp(0,0) = 1;
	_cv.setZero();
	_cv(0,1) = 1;
	_ca.setZero();
	_ca(0,2) = 1;		
		
	//predictive model matrixs: just calculated once
	_pps.setZero(); _ppu.setZero();
	_pvs.setZero(); _pvu.setZero();
	_pas.setZero(); _pau.setZero();
	
    _pps = Matrix_ps(_a,_nh,_cp);
	_pvs = Matrix_ps(_a,_nh,_cv);
	_pas = Matrix_ps(_a,_nh,_ca);
	_ppu = Matrix_pu(_a,_b,_nh,_cp);
	_pvu = Matrix_pu(_a,_b,_nh,_cv);
	_pau = Matrix_pu(_a,_b,_nh,_ca);
	
	_pvu_2 = _pvu.transpose()*_pvu;
	_ppu_2 = _ppu.transpose()*_ppu;
	
 //   cout<<"pps:"<<"\n"<<_pps<<endl; 
	//cout<<"pvs:"<<"\n"<<_pvs<<endl; 
	//cout<<"pas:"<<"\n"<<_pas<<endl; 
 //   cout<<"ppu:"<<"\n"<<_ppu<<endl; 
	//cout<<"pvu:"<<"\n"<<_pvu<<endl; 
	//cout<<"pau:"<<"\n"<<_pau<<endl; 
 //   cout<<"ppu2:"<<"\n"<<_ppu_2<<endl; 
	//cout<<"pvu2:"<<"\n"<<_pvu_2<<endl; 



	xyz1 = 0;  //flag for Indexfind function: 
	xyz2 = 1;
	_j_period = 0; // the number for step cycle indefind
	
        //footz refer:
	_Zsc.setZero();		
  	for (int i = 0; i < _nsum-1; i++) {	  
          Indexfind(_t(i),xyz1);	  
	  _Zsc(i) = _footz_ref(_j_period);   
	  _j_period = 0; 
	}		

    _yk.topRows(1).setConstant(_footy_ref(0)); 
	_zk.topRows(1).setConstant(_hcom);


	/// for footstep location reference generation: the foot-ref = _v_i*fx + _VV_i*_Lx_ref	
	_v_i.setZero();                    ///current step cycle
	_VV_i.setZero();	      /// next step cycles: 2 cycles maximal
	
	// optimized footstep location
	_footx_real.setZero();  
	_footy_real.setZero(); 
	_footz_real.setZero();	
	_footxyz_real.setZero();	
	_footx_real_next.setZero(); 
	_footy_real_next.setZero(); 
	_footz_real_next.setZero();
	_footx_real_next1.setZero();  
	_footy_real_next1.setZero(); 
	_footz_real_next1.setZero();
	
	/// for foot trajectory generation
    _Lfootx.setZero();
    _Lfooty.setZero(); _Lfooty.setConstant(_stepwidth(0)); _Lfootz.setZero(); 
	_Lfootvx.setZero(); _Lfootvy.setZero();_Lfootvz.setZero(); 
	_Lfootax.setZero(); _Lfootay.setZero();_Lfootaz.setZero();
	_Rfootx.setZero(); 
    _Rfooty.setZero(); _Rfooty.setConstant(-_stepwidth(0));_Rfootz.setZero(); 
	_Rfootvx.setZero(); _Rfootvy.setZero();_Rfootvz.setZero(); 
	_Rfootax.setZero(); _Rfootay.setZero();_Rfootaz.setZero();	
	_ry_left_right = 0;
///////=========================================constraints initialize========================================
	_ZMP_ratio = 0.9;	
	  //vertical height constraints	
	_z_max=0.1;
	_z_min=-0.1;	
	
	if(_robot_name == "coman"){
	  _rad = 0.1; 	  
	  _footx_max=0.3;
	  _footx_min=-0.2;	  	  
	  /// zmp-constraints	
	  _zmpx_ub=0.07;  
	  _zmpx_lb=-0.03;
	  _zmpy_ub=0.05; 
	  _zmpy_lb=-0.05;		  	  
	}
	else if (_robot_name == "bigman")
	{
	  _rad = 0.2; 	  
	  _footx_max=0.5;
	  _footx_min=-0.2;	  	  
	  /// zmp-constraints	
	  _zmpx_ub=0.07;  
	  _zmpx_lb=-0.04;
	  _zmpy_ub=(_foot_width/2*_ZMP_ratio); 
	  _zmpy_lb=(-_foot_width/2*_ZMP_ratio);		
	}
	else if (_robot_name == "cogimon")
        {	  
	  _rad = 0.2; 	  
	  _footx_max=0.4;
	  _footx_min=-0.2;	  	  
	  /// zmp-constraints	
	  _zmpx_ub=0.07;  
	  _zmpx_lb=-0.04;
	  _zmpy_ub=(_foot_width/2*_ZMP_ratio); 
	  _zmpy_lb=(-_foot_width/2*_ZMP_ratio);	
        } 
        else {
	  cout<<"Errorrrrrrrr for IK\n"<<endl;}		
	  
	_mass = _robot_mass; 		
	_j_ini = _mass* pow(_rad,2);		

	 
	_footy_max=2*_half_hip_width + 0.2; 
	_footy_min=_half_hip_width +0.02;
	
	// angle range
	_thetax_max=10*M_PI/180;  
	_thetax_min=-5*M_PI/180;
	_thetay_max=10*M_PI/180;  
	_thetay_min=-10*M_PI/180;
	
	// torque range
	_torquex_max=80/_j_ini; 
	_torquex_min=-60/_j_ini;
	_torquey_max=80/_j_ini;  
	_torquey_min=-80/_j_ini;	

	// swing foot velocity constraints	
	_footx_vmax=3;
	_footx_vmin=-2;
	_footy_vmax=2; 
	_footy_vmin=-2;	

	_fx= 0;    //current step location in local coordinate
	_fy= 0;
	_fxx_global= 0; //global footstep location in local coordinate
	_fyy_global= 0;		

///===========initiallize: preparation for MPC solution	
	// sulotion preparation	
	_V_ini.setZero();                                /// initialize optimal variable
    _V_inix.setZero();	
	_V_optimal.setZero();	      ///optimization det_V
	_flag.setZero();	        // 	 store n_vis: actual following step numbers under predictive window
	_flag_global.setZero();    // store step cycle sequence	
	
	_Lx_ref.setZero();             ///reference following footstep locations during the predictive window
	_Ly_ref.setZero(); 
	_Lz_ref.setZero();                    
	_comx_center_ref.setZero();
	_comy_center_ref.setZero();
	_comz_center_ref.setZero();	
	_thetax_center_ref.setZero(); 
	_thetay_center_ref.setZero();	
		  
	if(_robot_name == "coman"){
         //////// for methx ==2:reactive step + body inclination + height variance: for flat ground walking and up-down stairs: offline
	 _Rx = 1;           _Ry = 1;            _Rz =1;                     //com acceleration
	_alphax = 1;       _alphay = 1;        _alphaz = 100;               //com velocity
	_beltax = 5000;   _beltay = 10;        _beltaz = 20000000;          //com position
	_gamax =  10000000; _gamay = 10000000;  _gamaz = 200;               //footstep location
	_Rthetax = 1; _Rthetay = 1;                                         // theta acceleration
	_alphathetax =1; _alphathetay = 1;                                  // theta velocity
	_beltathetax = 10; _beltathetay = 10;	                            // theta postion	
	}
	else if(_robot_name  == "bigman"){
         //////// for methx ==2:reactive step + body inclination + height variance: for flat ground walking and up-down stairs: offline
	 _Rx = 1;           _Ry = 1;            _Rz =1;
	_alphax = 1;       _alphay = 10;        _alphaz = 100; 
	_beltax = 100;   _beltay = 1000;        _beltaz = 20000000;
	_gamax =  50000000; _gamay = 1000000000;  _gamaz = 200;
	_Rthetax = 1; _Rthetay = 1;
	_alphathetax =1; _alphathetay = 1;
	_beltathetax = 1000; _beltathetay = 1000;
	}
	else if (_robot_name == "cogimon"){
	  //////// for methx ==2:reactive step + body inclination + height variance: for flat ground walking and up-down stairs: offline !!!! to guarantee the solution: _alphay and _beltay should be large enough
         //desired motion: steplength=0.15
	  _Rx = 10;           _Ry = 10;            _Rz =10;
	  _alphax = 10;     _alphay = 50;        _alphaz = 100; 
	  _beltax = 40000;   _beltay = 20000;     _beltaz = 50000000;
	  _gamax =  50000000; _gamay = 50000000;  _gamaz = 200;
	  _Rthetax = 1; _Rthetay = 1;
	  _alphathetax =10; _alphathetay = 10;
	  _beltathetax = 1000000; _beltathetay = 1000000;	  	  
	  
        } 
	else
	{cout<<"Errorrrrrrrr for IK\n"<<endl;}
	
        cout<<"start mpc initialization3x"<<endl; 	
	
	// time cost consumption======================mark
	_tcpu.setZero(_nsum);

	///////////////////=====================================//////////////////////////
	/// for offline calculation	
	_loop = 2;   /// loop number for SQP
///// next code just run once	
	A_unit.setIdentity();
	C_unit.setIdentity();
	cout<<"C_unit:"<<"\n"<<C_unit<<endl;

	
  /////////// initialize each variable
    _bjxx_1 = 0; 
	_bjxx = 0; 
	_t_f.setZero();     ///predictive window time-period
	_bjx1 = 0;
	_bjx2 = 0;
    _mx = 0;
    _tnx.setZero();	  
	  
	_n_vis =0; 
	xxx = 0; 
	xxx1=0; 
	xxx2=0;	 
	
	// optimization objective function 	
	_WX.setZero();
	_WY.setZero();
	_WZ.setZero();
	_WthetaX.setZero();
	_WthetaY.setZero();
	_PHIX.setZero();
	_PHIY.setZero();
	_PHIZ.setZero();
	_Q_goal.setZero();
	_q_goal.setZero();
	_Q_goal1.setZero();
	_q_goal1.setZero();	
	
	_WX = _Rx*0.5 * A_unit + _alphax*0.5 * _pvu_2 + _beltax*0.5 * _ppu_2;	  
	_WY = _Ry/2 * A_unit + _alphay/2 * _pvu_2 + _beltay/2 * _ppu_2;
	_WZ = _Rz/2 * A_unit + _alphaz/2 * _pvu_2 + _beltaz/2 * _ppu_2;  
	_WthetaX = _Rthetax/2 * A_unit + _alphathetax/2 * _pvu_2 + _beltathetax/2 * _ppu_2;
	_WthetaY = _Rthetay/2 * A_unit + _alphathetay/2 * _pvu_2 + _beltathetay/2 * _ppu_2;
	_PHIX  = _gamax/2 * C_unit;
	_PHIY  = _gamay/2 * C_unit;
	_PHIZ  = _gamaz/2 * C_unit;
		
	_Q_goal.block<_nh, _nh>(0, 0) = _WX;
	_Q_goal.block<_nh, _nh>(_nh, _nh) = _WY;
	_Q_goal.block<_nh, _nh>(2*_nh, 2*_nh) = _WZ;
	_Q_goal.block<_nh, _nh>(3*_nh, 3*_nh) = _WthetaX;
	_Q_goal.block<_nh, _nh>(4*_nh, 4*_nh) = _WthetaY;
	_Q_goal.block<_nstep,_nstep>(5*_nh, 5*_nh) = _PHIX;
	_Q_goal.block<_nstep,_nstep>(5*_nh+_nstep, 5*_nh+_nstep) = _PHIY;
	_Q_goal.block<_nstep,_nstep>(5*_nh+2*_nstep, 5*_nh+2*_nstep) = _PHIZ;	  
      
	_Q_goal1 = 2 * _Q_goal;	
      
        cout<<"start mpc initialization4x"<<endl;       
	// constraints
	_Sjx.setZero();
	_Sjy.setZero();
	_Sjz.setZero();
	_Sjthetax.setZero();
	_Sjthetay.setZero();
	_Sjx.block<_nh, _nh>(0, 0) = A_unit;
	_Sjy.block<_nh, _nh>(0, _nh) = A_unit;
	_Sjz.block<_nh, _nh>(0, 2*_nh) = A_unit;	
	_Sjthetax.block<_nh, _nh>(0, 3*_nh) = A_unit;
	_Sjthetay.block<_nh, _nh>(0, 4*_nh) = A_unit;
	
	_Sfx.setZero();
	_Sfy.setZero();
	_Sfz.setZero();	 	
	

	// ZMP boundary preparation
	_H_q_upx.setZero();
	_F_zmp_upx.setZero();
	_H_q_lowx.setZero();
	_F_zmp_lowx.setZero();
	_H_q_upy.setZero();
	_F_zmp_upy.setZero();
	_H_q_lowy.setZero();
	_F_zmp_lowy.setZero();

	_phi_i_x_up.setZero();
	_p_i_x_t_up.setZero();
	_del_i_x_up.setZero();
	_p_i_x_t_low.setZero();
	_del_i_x_low.setZero();
	_phi_i_y_up.setZero();
	_p_i_y_t_up.setZero();
	_del_i_y_up.setZero();
	_p_i_y_t_low.setZero();
	_del_i_y_low.setZero();	
	_det_del_i_x.setZero();
	_det_del_i_y.setZero();
	

	// angle boundary preparation
	_q_upx.setZero();
	_qq_upx.setZero();
	_q_lowx.setZero();
	_qq_lowx.setZero();
	_q_upy.setZero();
	_qq_upy.setZero();
	_q_lowy.setZero();
	_qq_lowy.setZero();

	_qq1_upx.setZero();
	_qq1_lowx.setZero();
	_qq1_upy.setZero();
	_qq1_lowy.setZero();	  

	// torque bondary preparation
	_t_upx.setZero();
	_tt_upx.setZero();
	_t_lowx.setZero();
	_tt_lowx.setZero();
	_t_upy.setZero();
	_tt_upy.setZero();
	_t_lowy.setZero();
	_tt_lowy.setZero();

	_tt1_upx.setZero();
	_tt1_lowx.setZero();
	_tt1_upy.setZero();
	_tt1_lowy.setZero();

	// CoM height boundary preparation
	_H_h_upz.setZero();
	_F_h_upz.setZero();
	_H_h_lowz.setZero();
	_F_h_lowz.setZero();
	_delta_footz_up.setZero();
	_delta_footz_low.setZero();

	// CoM height acceleration boundary preparation
	_H_hacc_lowz.setZero();
	_F_hacc_lowz.setZero();
	_delta_footzacc_up.setZero();	  


	//swing foot velocity constraints
	_Footvx_max.setZero();
	_Footvx_min.setZero();
	_Footvy_max.setZero();
	_Footvy_min.setZero();
	_footubxv.setZero();
	_footlbxv.setZero();
	_footubyv.setZero();
	_footlbyv.setZero();
	
	// foot location constraints: be careful that the step number is change: so should be intialized in each whole loop
	_H_q_footx_up.setZero();
	_F_foot_upx.setZero();
	_H_q_footx_low.setZero();
	_F_foot_lowx.setZero();
	_H_q_footy_up.setZero();
	_F_foot_upy.setZero();
	_H_q_footy_low.setZero();
	_F_foot_lowy.setZero();	
	
	

	// foot vertical location-equality constraints
	_H_q_footz.setZero();
	_F_footz.setZero();

	// CoMZ height-equality constraints
	_h_h.setZero();
	_hhhx.setZero();	  

	// body inclination-equality constraints
	_a_hx.setZero();
	_a_hxx.setZero();
	_a_hy.setZero();
	_a_hyy.setZero();


	// foot location constraints
	_Sfoot.setZero();
	_Sfoot(0) = -1;
	_Sfoot(1) = 1;
	  
        
    /// global CoMx and CoMy of each loop
    _Footx_global_relative =0;
    _Footy_global_relative =0;	


	_ZMPx_constraints_half2.setZero();	 
	_ZMPy_constraints_half2.setZero();
	_phi_i_x_up1.setZero();
	_phi_i_y_up1.setZero();
	
	//// data saving
	CoMMM_ZMP_foot.setZero();
	
	///QP size initiallize
	if (_method_flag ==0)
	{
	  
	  int nVars = _Nt;
	  int nEqCon = 1+3*_nh;
	  int nIneqCon = 5*_nh + 4*_nstep+4;  
	  resizeQP(nVars, nEqCon, nIneqCon);		   
        }
	else if (_method_flag ==1)
	{
	  int nVars = _Nt;
	  int nEqCon = 1+_nh;
	  int nIneqCon = 13*_nh + 4*_nstep +4;
	  resizeQP(nVars, nEqCon, nIneqCon);	   
	}
	else
	{
	  int nVars = _Nt;
	  int nEqCon = 1;
	  int nIneqCon = 15*_nh + 4*_nstep +4;
	  resizeQP(nVars, nEqCon, nIneqCon);		   
	}

	_t_end_walking = _tx(_footstepsnumber-1)- 9*_tstep/4;  // ending time when normal walking ends
	_n_end_walking = (int)boost::math::iround(_t_end_walking/_dt);	
	_comy_matrix_inv.setZero();	
	
	
    cout<<"start mpc initialization5x"<<endl; 
               
        
        
	//// offline calculated matrices
	
	_ppu_T = _ppu.transpose();
        std::cout << "MPC_initialization1X...\n"; 
        
	_pvupvs = (_pvu.transpose()) * _pvs;
        std::cout << "MPC_initialization2X...\n";         
        _ppupps = (_ppu_T) * _pps;
	
	std::cout << "MPC_initialization3X...\n"; 
        
        
	//// offline calulated the ZMP constraints coefficient==================================
	  ////////initiallize vector
	vector <Eigen::Matrix<double,_Nt, _Nt>,  Eigen::aligned_allocator<Eigen::Matrix<double,_Nt,_Nt>>> x_offline1(_nh)  ;
        _x_offline1_va.setZero();               
	for (int j=0;j<_nh; j++)
	{
	  x_offline1[j]= _x_offline1_va;
	}
	ZMPx_constraints_offfline = x_offline1;
	ZMPy_constraints_offfline = x_offline1;		
	_phi_i_x_up_est = x_offline1;
	_phi_i_y_up_est = x_offline1;

        
	vector <Eigen::Matrix<double,_nh, _Nt>, Eigen::aligned_allocator<Eigen::Matrix<double,_nh,_Nt>>> x_offline2(_nh)  ;
        _x_offline2_va.setZero();              
	for (int j=0;j<_nh; j++)
	{
	  x_offline2[j]= _x_offline2_va;
	}	
	ZMPx_constraints_half = x_offline2;	
	ZMPy_constraints_half = x_offline2;
    
    ZMPx_constraints_halfxxx1.setZero(); 
    ZMPx_constraints_halfxxx2.setZero();  
  
    _x_offline4_vax.setZero();   	
    _x_offline4_vay.setZero();  
      
    ZMPx_constraints_halfyyy1.setZero(); 
    ZMPx_constraints_halfyyy2.setZero();          
    ZMPx_constraints_halfyyy3.setZero(); 
    ZMPx_constraints_halfyyy4.setZero();         
        
    _matrix_large1.setZero();
    _matrix_large2.setZero();
    _ZMPx_constraints_half_va.setZero();
    _ZMPy_constraints_half_va.setZero();
	
    std::cout << "local vector x_offline...\n";         
		
	for(int jxx=1; jxx<=_nh; jxx++)
	{
	  _Si.setZero();
	  _Si(0,jxx-1) = 1;
	  // ZMP constraints	      		 
	 //ZMPx_constraints_half[jxx-1] = - (_Si).transpose() * _Si * _pau * _Sjz;
	 //ZMPy_constraints_half[jxx-1] = - (_Si).transpose() * _Si * _pau * _Sjz;
 
		ZMPx_constraints_half[jxx-1] = - (_Si).transpose() * _Si * _pau * _Sjz;
		ZMPy_constraints_half[jxx-1] = ZMPx_constraints_half[jxx-1];
         
		ZMPx_constraints_halfxxx1 = _Si * _pau * _Sjz;
		ZMPx_constraints_halfxxx2 = _Si * _ppu * _Sjz;
		ZMPx_constraints_halfyyy1 = (_Si * _ppu * _Sjx).transpose();
		ZMPx_constraints_halfyyy2 = (_Si * _pau * _Sjx).transpose();         
		ZMPx_constraints_halfyyy3 = (_Si * _ppu * _Sjy).transpose();
		ZMPx_constraints_halfyyy4 = (_Si * _pau * _Sjy).transpose();          
         
			 /// write a function to deal with the matrix multiply	 
		 Matrix_large(1);        
		 ZMPx_constraints_offfline[jxx-1] = _matrix_large1 - _matrix_large2; 	 
		_matrix_large1.setZero();
		_matrix_large2.setZero();
         
		Matrix_large(2);	  
		 ZMPy_constraints_offfline[jxx-1] = _matrix_large1 - _matrix_large2;
		_matrix_large1.setZero();
		_matrix_large2.setZero();
	//	  ZMPx_constraints_offfline[jxx-1] = ZMPx_constraints_halfyyy1 * ZMPx_constraints_halfxxx1 - ZMPx_constraints_halfyyy2 * ZMPx_constraints_halfxxx2;                  
	// 	  ZMPy_constraints_offfline[jxx-1] = ZMPx_constraints_halfyyy3 * ZMPx_constraints_halfxxx1 - ZMPx_constraints_halfyyy4 * ZMPx_constraints_halfxxx2;

	}


	vector <Eigen::Matrix<double,3, _Nt>, Eigen::aligned_allocator<Eigen::Matrix<double,3,_Nt>>> xx_offline1(_nh);
        _x_offline3_va.setZero();
        
	for (int j=0;j<_nh; j++)
	{
	  xx_offline1[j]= _x_offline3_va;
	}

	_xkZMPx_constraints = xx_offline1;
	_zkZMPx_constraints = xx_offline1;
	_ykZMPy_constraints = xx_offline1;
	_zkZMPy_constraints = xx_offline1;

	
	for(int jxx=1; jxx<=_nh; jxx++)
	{
 	  _xkZMPx_constraints[jxx-1] = (_pps.row(jxx-1)).transpose()*_pau.row(jxx-1)*_Sjz - (_pas.row(jxx-1)).transpose()* _ppu.row(jxx-1)*_Sjz;
 	  _zkZMPx_constraints[jxx-1] = (_pas.row(jxx-1)).transpose()*_ppu.row(jxx-1)*_Sjx - (_pps.row(jxx-1)).transpose() *_pau.row(jxx-1)*_Sjx;
 	  _zkZMPy_constraints[jxx-1] = (_pas.row(jxx-1)).transpose()*_ppu.row(jxx-1)*_Sjy - (_pps.row(jxx-1)).transpose() *_pau.row(jxx-1)*_Sjy;	  
 	  
	}
	_ykZMPy_constraints = _xkZMPx_constraints;
	
	
	
	_ppuSjx.setZero();
	_ppuSjx = _ggg(0,0)*_ppu*_Sjx;
	_ppuSjy.setZero();
	_ppuSjy = _ggg(0,0)*_ppu*_Sjy;	
	
	_pauSjx.setZero();
	_pauSjx = _pau*_Sjx;
	_pauSjy.setZero();
	_pauSjy = _pau*_Sjy;	
	
	
	_pauSjz1.setZero();
	_pauSjz1 = _zmpx_ub*_pau*_Sjz;
	_pauSjz2.setZero();   ////error ( _p_i_x_t_low -_p_i_x_t_up)   
	_pauSjz2 = _mass*_zmpx_lb*_pau*_Sjz -_pauSjz1;	
	_pauSjz11.setZero();
	_pauSjz11 = _zmpy_ub*_pau*_Sjz;
	_pauSjz21.setZero(); /// error  ( _p_i_x_t_low -_p_i_x_t_up)
	_pauSjz21 = _mass*_zmpy_lb*_pau*_Sjz - _pauSjz11;	
	
	
	_pauSjthetay.setZero();
	_pauSjthetay = _j_ini * _pau * _Sjthetay;
	_pauSjthetax.setZero();
	_pauSjthetax = _j_ini * _pau * _Sjthetax;


	vector <Eigen::Matrix3d> xxx_offline1xx(_nh);
    _x_offline4_va.setZero();
	for (int j=0;j<_nh; j++)
	{
	  xxx_offline1xx[j]= _x_offline4_va;
	}
	_xkzk_constraints = xxx_offline1xx;

	for(int jxx=1; jxx<=_nh; jxx++)
	{
	  _xkzk_constraints[jxx-1] = (_pps.row(jxx-1)).transpose() *_pas.row(jxx-1) - (_pas.row(jxx-1)).transpose() *_pps.row(jxx-1);

	}	
	
	_gzmpxub = _ggg *_zmpx_ub;
	_gzmpxlb = _ggg *_zmpx_lb;	
	_gzmpyub = _ggg *_zmpy_ub;
	_gzmpylb = _ggg *_zmpy_lb;		
	
      
	
	//angle range constraints: 1 order coefficient
	_q_upx = _ppu* _Sjthetax;
	_q_lowx = -_q_upx;	         

	_q_upy = _Sjthetay;
	_q_lowy = -_q_upy;	
	//torque range constraints
	_t_upx = _pau* _Sjthetax;
	_t_lowx = -_t_upx;
	
	_t_upy = _pau* _Sjthetay;
	_t_lowy = -_t_upy;	 		
	// body height constraints
	_H_h_upz = _ppu* _Sjz;
	_H_h_lowz = -_H_h_upz;   	
	
	// body height acceleration constraints	      
	_H_hacc_lowz = -_pau* _Sjz;   

	//equality constraints: footz height constraints:   
	_h_h = _ppu * _Sjz;
	
	_a_hx = _ppu * _Sjthetax;
	_a_hy = _ppu * _Sjthetay;	
    
    
        _detzmppxb = _zmpx_ub- _zmpx_lb;
	_detzmppyb = _zmpy_ub- _zmpy_lb;
	_detthetax = _thetax_max-_thetax_min;
	_detthetay = _thetay_max-_thetay_min;
	_detfootx  = _footx_max - _footx_min; 
	_detfooty  = _footy_max - _footy_min; 
	_detfootvx = _footx_vmax - _footx_vmin;
	_detfootvy = _footy_vmax - _footy_vmin;
	_dettorquex = _torquex_max - _torquex_min;
	_dettorquey = _torquey_max - _torquey_min;
	_detzb = _z_max - _z_min;


	
	_nTdx =0;

	//////polynomial intepolation for lower level interpolation
	_AAA_inv.setZero();

	
	
	_j_count = 0;
	
	
	/////////////for ZMP distribution
// 	_F_R.setZero();  _F_L.setZero();   _M_R.setZero(); _M_L.setZero();
	_F_R(2) = _F_L(2) = 0.5 * _mass * _GGG;
	_Co_L.setZero(); _Co_R.setZero();
	
	_comxyzx.setZero(); _comvxyzx.setZero(); _comaxyzx.setZero(); 
	_thetaxyx.setZero(); _thetavxyx.setZero(); _thetaaxyx.setZero();
	_Lfootxyzx.setZero(); _Rfootxyzx.setZero();
	_ZMPxy_realx.setZero();
	
	_comxyzx(2) = _HCOM;			
}



///////////////////////// local coordinate CoM solution---modified---------------------------------
void MPCClass::CoM_foot_trajection_generation_local(int i, const Eigen::Matrix<double,18,1>& estimated_state, Eigen::Vector3d _Rfoot_location_feedback, Eigen::Vector3d _Lfoot_location_feedback,double lamda, bool _stopwalking)
{
 
    if (i<_n_end_walking)   ///////normal walking
    {
	  /// modified the footy_min
	  if (i==((int)boost::math::iround(2*_ts(1)/_dt))+1) ///update the footy_limit
	  {
	      _footy_min = _foot_width+0.01;
	      _detfooty  = _footy_max - _footy_min; 
	  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
  // 	================ iterative calulcation: predictive control_tracking with time-varying height+angular momentum	    	    
	    /// run once
	    ///////////////////////////////////////////////////////	    

	    Indexfind((i-1)*_dt,xyz1);	  	      //// step cycle number when (i-1)*dt fall into      
	    _bjxx_1 = _j_period+1;
	    _j_period = 0;
	    
	    
	    Indexfind(i*_dt,xyz1);                   //// step cycle number when (i)*dt fall into : current sampling time
	    _bjxx = _j_period+1;  //coincidence with matlab 
	    _j_period = 0;	  
	    
	    // com_center_ref = ZMP_center_ref = v_i*f + V_i*L_ref
	    //solve the following steps: 1.5s may couve 2 or three  steps, so  one/two following steps	    
	    _t_f.setLinSpaced(_nh,(i+1)*_dt, (i+_nh)*_dt);
	    
	    Indexfind(_t_f(0),xyz1);                /// step cycle number when (i+1)*dt fall into : current sampling time
	    _bjx1 = _j_period+1;
	    _j_period = 0;
	    
	    Indexfind(_t_f(_nh-1),xyz1);           /// step cycle number when (i+_nh)*dt fall into : current sampling time
	    _bjx2 = _j_period+1;
	    _j_period = 0;	  
	    	    
	    ////================================================================
	    /// judge if stop walking enable: if (from _bjx1 step, reference step length and step height will be set to be zero)
	    if(_stopwalking)  
	    {    
	      for (int i_t = _bjx1; i_t < _footstepsnumber; i_t++) {
		_steplength(i_t) = 0;	
		_footx_ref(i_t) = _footx_ref(i_t-1) + _steplength(i_t-1); 	      
	      }	  
	    }
	    
	    _mx = _bjx2 - _bjx1 +1;
            /// find out the relative postion of start-time of predictive step cycles to the predictive window 	    
	    for (int j=1;j<_mx; j++)
	    {
	      Indexfind(_tx(_bjx1+j-1),xyz2);
	      _tnx(j-1) = _j_period; 
	      _j_period = 0;	      
	    }

		//c++98 from double to int===> not accurate
	    _v_i.setZero();
	    _VV_i.setZero();
	    xxx = (int)_tnx(0);
	    // be careful that the position is from 0;;;;;         
	    if (fabs(_tnx(0) - _nT) <=0.00001)
	    {
	      _n_vis =2;
	      for (int jjj = 1; jjj <= _mx; jjj++)
	      {
		
		if (jjj == 1)
		{
		  _VV_i.block(0, 0, xxx, 1).setOnes();   ///dynamic size 
		}
		else
		{	
		  xxx1 = _nh-(int)_tnx(0);
		  _VV_i.block(xxx, 1, xxx1, 1).setOnes();		
		}
	      }	    
	    }
	    else
	    {	
	      _v_i.head(xxx).setOnes();	      
/*	      _v_i.segment(0, xxx).setOnes();*/	      
	      if (abs(_mx - 2) <=0.00001)
	      {
		_n_vis = 1;
		_VV_i.block(xxx, 0, _nh -xxx, 1).setOnes();
	      }
	      else
	      {
		_n_vis = 2;
		xxx2 = _nh -(int)_tnx(1);
		_VV_i.block<_nT, 1>(xxx, 0).setOnes();
		_VV_i.block((int)_tnx(1), 1, xxx2, 1).setOnes();	      	      
	      }
	    }
	     
	    _flag(i-1,0)= _n_vis;	  ////_n_vis flag: for judging if there is a switch between current support leg for i-1 and i sampling time    
	    _flag_global(i-1,0) = _bjxx;  // walking cycle of current sampling time
	    _flag_global(i,0) = _bjx1;	  // walking cycle of next one sampling time 
  
  //================cuttent state switch to the current support; ====================//	    
	    ///// pass the actual state into the control loop:
	    /// output _Footx_global_relative and _Footy_global_relative are the global xk(0) and yk(0)
	    if (i>1)
	    {
	      _Footx_global_relative = _xk(0,i-1) + _fxx_global;    ///_fxx_global is the old one: 
	      _Footy_global_relative = _yk(0,i-1) + _fyy_global;	    
	    }	     
	    
	    // current foot location
	    _fx =0;
	    _fy = 0;	  	    
	    _fxx_global = _footx_real(_bjxx-1);      ////update the supporting foot
	    _fyy_global = _footy_real(_bjxx-1);	     ////update the supporting foot

  // 	      // relative state switch	: the output xk(0) and yk(0) is relative to the current foot location      
	    if (i>1)
	    {
	      if (_flag_global(i-2,0) < _flag_global(i-1,0) )  //current and last support foot are not the same
	      {
		/// reference relative state switch
		_xk(0,i-1) = _Footx_global_relative - _fxx_global; 
		_yk(0,i-1) = _Footy_global_relative - _fyy_global;		
	      }
	    }
	    
    
            /////adjust the reference step length and step width	    
	    if (_n_vis ==1)
	    {
	      _Lx_ref(0) = _footx_ref(_bjx2-1) - _fxx_global;  //relative to the current location ///here, position tracking,thus, use the reference location postion:_footx_ref
	      _Ly_ref(0) = _footy_ref(_bjx2-1) - _fyy_global;  // relative to the current location//here,position tracking, thus, use the reference location postion:_footx_ref
	      _Lz_ref(0) = _footz_ref(_bjx2-1);
	      _Lx_ref(1) = 0;
	      _Ly_ref(1) = 0;
	      _Lz_ref(1) = 0;	    
	    }
	    else
	    {
	      _Lx_ref(0) = _footx_ref(_bjx2-2) - _fxx_global;
	      _Ly_ref(0) = _footy_ref(_bjx2-2) - _fyy_global;
	      _Lz_ref(0) = _footz_ref(_bjx2-2);
	      _Lx_ref(1) = _footx_ref(_bjx2-1) - _fxx_global;
	      _Ly_ref(1) = _footy_ref(_bjx2-1) - _fyy_global;
	      _Lz_ref(1) = _footz_ref(_bjx2-1);	    
	    }	    
	  // com_center_ref: set to be the center of supporting foot
	    _comx_center_ref = _v_i*_fx + _VV_i*_Lx_ref;
	    _comy_center_ref = _v_i*_fy + _VV_i*_Ly_ref;
	    _comz_center_ref = _Zsc.segment<_nh>(i) + _Hcom;
	    
	    /// hot start:initilize _V_ini with last loop results
	    if (i==1)
	    {
	      _V_ini(5*_nh) = _footx_ref(1);
	      _V_ini(5*_nh+1) = _footy_ref(1);	    
	    }
	    else
	    {
	      _V_ini.topRows(5*_nh) = _V_optimal.block<5*_nh, 1>(0, i-2);
	      if (_n_vis > _flag(i-2))
	      { 
			_V_ini(_Nt -1-1) = _V_optimal(5*_nh+6-1, i-2);
			_V_ini(_Nt -3-1) = _V_optimal(5*_nh+6-2-1, i-2);
			_V_ini(_Nt -5-1) = _V_optimal(5*_nh+6-4-1, i-2); 
	      }
	      else
	      {
			if (_n_vis < _flag(i-2))
			{ 
			  _V_ini(_Nt-1) = _V_optimal(5*_nh+6-1, i-2);
			  _V_ini(_Nt -1-1) = _V_optimal(5*_nh+6-2-1, i-2);
			  _V_ini(_Nt -2-1) = _V_optimal(5*_nh+6-4-1, i-2); 
			}	   
			else
			{
			  if (_n_vis ==1)
			  { 
				_V_ini(_Nt-1) = _V_optimal(5*_nh+6-1, i-2);
				_V_ini(_Nt -1-1) = _V_optimal(5*_nh+6-2-1, i-2);
				_V_ini(_Nt -2-1) = _V_optimal(5*_nh+6-4-1, i-2); 
			  }
			  else
			  {
				_V_ini.bottomRows(6) = _V_optimal.block<6, 1>(5*_nh+6-6, i-2);
			  }
			}
	      }
	    }
	      	  
	    // foot location constraints: be careful that the step number is change: so should be setZero in each whole loop
	    _H_q_footx_up.setZero();
	    _F_foot_upx.setZero();
	    _H_q_footx_low.setZero();
	    _F_foot_lowx.setZero();
	    _H_q_footy_up.setZero();
	    _F_foot_upy.setZero();
	    _H_q_footy_low.setZero();
	    _F_foot_lowy.setZero();	    

	    // boundary initialzation: footx, footy,footz select matrix: 2*(5*_nh+3*_nstep)
	    _Sfx.setZero();
	    _Sfy.setZero();
	    _Sfz.setZero();	    
	    if (_n_vis ==1)
	    {
	      _Sfx(0,5*_nh) = 1;
	      _Sfy(0,5*_nh+_nstep) = 1;
	      _Sfz(0,5*_nh+2*_nstep) = 1;
	    }  
	    else
	    {
	      _Sfx(0,5*_nh) = 1;
	      _Sfx(1,5*_nh+1) = 1;
	      _Sfy(0,5*_nh+_nstep) = 1;
	      _Sfy(1,5*_nh+_nstep+1) = 1;
	      _Sfz(0,5*_nh+2*_nstep) = 1;	 
	      _Sfz(1,5*_nh+2*_nstep+1) = 1;		      
	    }

	    //SQP MOdels	 
	      _q_goal.block<_nh, 1>(0, 0) = (_alphax * _pvupvs + _beltax * _ppupps )* _xk.col(i-1) - _beltax * _ppu_T * _comx_center_ref;
	      _q_goal.block<_nh, 1>(_nh, 0) = (_alphay * _pvupvs  + _beltay * _ppupps )* _yk.col(i-1) - _beltay * _ppu_T * _comy_center_ref;
	      _q_goal.block<_nh, 1>(2*_nh, 0) = (_alphaz * _pvupvs + _beltaz * _ppupps )* _zk.col(i-1) - _beltaz * _ppu_T * _comz_center_ref;
	      _q_goal.block<_nh, 1>(3*_nh, 0) = (_alphathetax * _pvupvs + _beltathetax * _ppupps )* _thetaxk.col(i-1) - _beltathetax * _ppu_T * _thetax_center_ref;
	      _q_goal.block<_nh, 1>(4*_nh, 0) = (_alphathetay * _pvupvs + _beltathetay * _ppupps )* _thetayk.col(i-1) - _beltathetay * _ppu_T * _thetay_center_ref;
	      _q_goal.block<_nstep, 1>(5*_nh, 0) = -_gamax * _Lx_ref;
	      _q_goal.block<_nstep, 1>(5*_nh+_nstep, 0) = -_gamay * _Ly_ref;
	      _q_goal.block<_nstep, 1>(5*_nh+2*_nstep, 0) = -_gamaz * _Lz_ref;
	
	    ///// the following code only run once in each loop:   
	      for(int jxx=1; jxx<=_nh; jxx++)
	      {
		// ZMP constraints: merger the simliar item and calculated offline
		// x-ZMP upper boundary                                      
	        _p_i_x_t_up.row(jxx-1) = _mass * (  (_xk.col(i-1)).transpose() * _xkZMPx_constraints[jxx-1] + (_zk.col(i-1)).transpose()*_zkZMPx_constraints[jxx-1] + _ppuSjx.row(jxx-1) + _Zsc.row(i+jxx-1)*_pauSjx.row(jxx-1) - ( (_pas.row(jxx-1) * _zk.col(i-1)).transpose() *_VV_i.row(jxx-1)*_Sfx ) - _ggg*_VV_i.row(jxx-1)*_Sfx - _pauSjz1.row(jxx-1) ) - _pauSjthetay.row(jxx-1);		
		_del_i_x_up.col(jxx-1) = _mass * ((_xk.col(i-1)).transpose()*_xkzk_constraints[jxx-1]*_zk.col(i-1) + _ggg*_pps.row(jxx-1) * _xk.col(i-1) + (_pas.row(jxx-1) * _xk.col(i-1)).transpose() *_Zsc.row(i+jxx-1) - _zmpx_ub*_pas.row(jxx-1)*_zk.col(i-1) - _gzmpxub) - _j_ini * _pas.row(jxx-1) * _thetayk.col(i-1);

		// x-ZMP low boundary
		_p_i_x_t_low.row(jxx-1) = _p_i_x_t_up.row(jxx-1) + _pauSjz2.row(jxx-1);	      
		_del_i_x_low.col(jxx-1) = _del_i_x_up.col(jxx-1) +_mass*(_detzmppxb*_pas.row(jxx-1)*_zk.col(i-1)+  _gzmpxub- _gzmpxlb);
		
		// y-ZMP upper boundary
		_p_i_y_t_up.row(jxx-1) = _mass * ( (_yk.col(i-1)).transpose() *_ykZMPy_constraints[jxx-1] + (_zk.col(i-1)).transpose()*_zkZMPy_constraints[jxx-1] + _ppuSjy.row(jxx-1) + _Zsc.row(i+jxx-1)*_pauSjy.row(jxx-1) - ((_pas.row(jxx-1) * _zk.col(i-1)).transpose() *_VV_i.row(jxx-1)*_Sfy ) - _ggg*_VV_i.row(jxx-1)*_Sfy - _pauSjz11.row(jxx-1)) + _pauSjthetax.row(jxx-1);
		_del_i_y_up.col(jxx-1) = _mass * ( (_yk.col(i-1)).transpose() *_xkzk_constraints[jxx-1]*_zk.col(i-1) + _ggg*_pps.row(jxx-1) * _yk.col(i-1) + (_pas.row(jxx-1) * _yk.col(i-1)).transpose() *_Zsc.row(i+jxx-1) - _zmpy_ub*_pas.row(jxx-1)*_zk.col(i-1) - _gzmpyub) + _j_ini * _pas.row(jxx-1) * _thetaxk.col(i-1); 	      
	
		// y-ZMP low boundary
		_p_i_y_t_low.row(jxx-1) = _p_i_y_t_up.row(jxx-1) + _pauSjz21.row(jxx-1);	      
		_del_i_y_low.col(jxx-1) = _del_i_y_up.col(jxx-1) +_mass*(_detzmppyb*_pas.row(jxx-1)*_zk.col(i-1)+  _gzmpyub -_gzmpylb);	      	      	     

		
		
		//angle range constraints
	         
                ////// attention : don't forger to test on the robot pC: if  _qq1_upx = _pps* _thetaxk.col(i-1) can run in the real time!!!!!!
		_qq1_upx.row(jxx-1) = _pps.row(jxx-1)* _thetaxk.col(i-1);
		_qq1_upx(jxx-1,0) = _qq1_upx(jxx-1,0)-_thetax_max;
			
		
		_qq1_upy.row(jxx-1) = _pps.row(jxx-1)* _thetayk.col(i-1);
		_qq1_upy(jxx-1,0) = _qq1_upy(jxx-1,0)-_thetay_max;    

		//torque range constraints
		
		_tt1_upx.row(jxx-1) = _pas.row(jxx-1)* _thetaxk.col(i-1);
		_tt1_upx(jxx-1,0) = _tt1_upx(jxx-1,0)-_torquex_max;
			 
		
		_tt1_upy.row(jxx-1) = _pas.row(jxx-1)* _thetayk.col(i-1);
		_tt1_upy(jxx-1,0) = _tt1_upy(jxx-1,0)-_torquey_max;		
		
		// body height constraints 	
		_delta_footz_up.row(jxx-1) = _pps.row(jxx-1)*_zk.col(i-1) - _comz_center_ref.row(jxx-1) ;
		_delta_footz_up(jxx-1,0) = _delta_footz_up(jxx-1,0) - _z_max;
		
		// body height acceleration constraints	      
		_delta_footzacc_up.row(jxx-1) = _pas.row(jxx-1)*_zk.col(i-1) + _ggg;
	      }

  ///       only one-time caluation	  
	    _ZMPx_constraints_half2 = (_VV_i* _Sfx).transpose();
	    _ZMPy_constraints_half2 = (_VV_i* _Sfy).transpose();	    

	    for(int jxx=1; jxx<=_nh; jxx++)
	    {
	      // ZMP constraints
	      // x-ZMP upper boundary
//	      _phi_i_x_up1 = ZMPx_constraints_offfline[jxx-1] + _ZMPx_constraints_half2 * ZMPx_constraints_half[jxx-1];	      
//	      _phi_i_x_up_est[jxx-1] = _mass * (_phi_i_x_up1 + _phi_i_x_up1.transpose())/2;
    
	      // y-ZMP upper boundary
//	      _phi_i_y_up1 = ZMPy_constraints_offfline[jxx-1] + _ZMPy_constraints_half2 * ZMPy_constraints_half[jxx-1];
//	      _phi_i_y_up_est[jxx-1] = _mass * (_phi_i_y_up1 + _phi_i_y_up1.transpose())/2;   
              
 //           _phi_i_x_up1 = ZMPx_constraints_offfline[jxx-1] + _ZMPx_constraints_half2 * ZMPx_constraints_half[jxx-1];   
//            _phi_i_y_up1 = ZMPy_constraints_offfline[jxx-1] + _ZMPy_constraints_half2 * ZMPy_constraints_half[jxx-1];   
                
                _ZMPx_constraints_half_va = ZMPx_constraints_half[jxx-1];
                _ZMPy_constraints_half_va = ZMPy_constraints_half[jxx-1];

                Matrix_large(3);
                _phi_i_x_up1 = ZMPx_constraints_offfline[jxx-1] +_matrix_large1;
                _phi_i_x_up1 = ZMPx_constraints_offfline[jxx-1] +_matrix_large1; 

                _phi_i_x_up_est[jxx-1] = _mass * (_phi_i_x_up1 + _phi_i_x_up1.transpose())/2;
                _phi_i_y_up_est[jxx-1] = _mass * (_phi_i_y_up1 + _phi_i_y_up1.transpose())/2; 
                
                _matrix_large1.setZero();
                _matrix_large2.setZero();
	    }
    
	    // constraints:Swing+Foot velocity: 
	    _Footvx_max = _Sfx.row(0);
	    _Footvx_min = -_Sfx.row(0);
	    _Footvy_max = _Sfy.row(0);
	    _Footvy_min = -_Sfy.row(0);	  		    
	    
	    //equality constraints: footz height constraints:
	    _H_q_footz = _Sfz.row(0);	    
// 	    _h_h = _ppu * _Sjz;
// 	    
// 	    _a_hx = _ppu * _Sjthetax;
// 	    _a_hy = _ppu * _Sjthetay;
	   
	    // SEQUENCE QUADARTIC PROGRAMMING: lOOP_until the maximal loops reaches		
	  // calculated the control loop

            _det_del_i_x = _del_i_x_low -  _del_i_x_up;
	    _det_del_i_y = _del_i_y_low -  _del_i_y_up;
	    
	    for (int xxxx=1; xxxx <= _loop; xxxx++)
	    {	
              /// attention: V_ini would be updated in each loop.
	      _q_goal1 = _Q_goal1 * _V_ini + _q_goal;	  
		      
	      
	      /// time consuming process	    
	      
	      for(int jxx=1; jxx<=_nh; jxx++)
	      {
			// x-ZMP upper constraints
			_phi_i_x_up = _V_ini.transpose() *_phi_i_x_up_est[jxx-1];
			_H_q_upx.row(jxx-1) = 2*_phi_i_x_up + _p_i_x_t_up.row(jxx-1); 
			_F_zmp_upx.row(jxx-1) = -((_phi_i_x_up + _p_i_x_t_up.row(jxx-1)) * _V_ini + _del_i_x_up.col(jxx-1)); 

			// x-ZMP low boundary   	      
			_H_q_lowx.row(jxx-1) = - _pauSjz2.row(jxx-1) - _H_q_upx.row(jxx-1);  
			_F_zmp_lowx.row(jxx-1) = _pauSjz2.row(jxx-1) * _V_ini + _det_del_i_x.col(jxx-1)-_F_zmp_upx.row(jxx-1); 
		
		
			// y-ZMP upper boundary	      
			_phi_i_y_up = _V_ini.transpose()*_phi_i_y_up_est[jxx-1];	      
			_H_q_upy.row(jxx-1) = 2*_phi_i_y_up + _p_i_y_t_up.row(jxx-1); 
			_F_zmp_upy.row(jxx-1) = -((_phi_i_y_up + _p_i_y_t_up.row(jxx-1)) * _V_ini + _del_i_y_up.col(jxx-1)); 
		
			// y-ZMP low boundary	      	      
			_H_q_lowy.row(jxx-1) = - _pauSjz21.row(jxx-1) - _H_q_upy.row(jxx-1); 
			_F_zmp_lowy.row(jxx-1) = _pauSjz21.row(jxx-1) * _V_ini + _det_del_i_y.col(jxx-1)-_F_zmp_upy.row(jxx-1);      	      
	    

		
			//angle range constraints
			_qq_upx.row(jxx-1) = -(_q_upx.row(jxx-1)* _V_ini + _qq1_upx.row(jxx-1));	 
			_qq_lowx(jxx-1,0) = _detthetax - _qq_upx(jxx-1,0);	 
				
			_qq_upy.row(jxx-1) = -(_q_upy.row(jxx-1)* _V_ini + _qq1_upy.row(jxx-1));
			_qq_lowy(jxx-1,0) = _detthetay - _qq_upy(jxx-1,0);

			//torque range constraints	      
			_tt_upx.row(jxx-1) = -(_t_upx.row(jxx-1)* _V_ini +  _tt1_upx.row(jxx-1));
			_tt_lowx(jxx-1,0) = _dettorquex - _tt_upx(jxx-1,0);	
		
			_tt_upy.row(jxx-1) = -(_t_upy.row(jxx-1)* _V_ini +  _tt1_upy.row(jxx-1));	      
			_tt_lowy(jxx-1,0) = _dettorquey - _tt_upy(jxx-1,0);	
		
			// body height constraints	      
			_F_h_upz.row(jxx-1) = -(_H_h_upz.row(jxx-1)*_V_ini + _delta_footz_up.row(jxx-1));      
			_F_h_lowz(jxx-1,0) = _detzb - _F_h_upz(jxx-1,0);
		
			// body height acceleration constraints	      
			_F_hacc_lowz.row(jxx-1) = (-_H_hacc_lowz.row(jxx-1)*_V_ini + _delta_footzacc_up.row(jxx-1));	      	      
	      }
	      
	      // foot location constraints
	      if (_n_vis == 1)  //one next steo
	      {
			_H_q_footx_up.row(0) = _Sfx.row(0);
			_F_foot_upx.row(0) = -(_H_q_footx_up.row(0) * _V_ini); 
			_F_foot_upx(0,0) = _F_foot_upx(0,0) +_fx + _footx_max; 
		
			_H_q_footx_low.row(0) = -_Sfx.row(0);
	// 		_F_foot_lowx.row(0) = (_H_q_footx_up.row(0) * _V_ini);
			_F_foot_lowx(0,0) = -_F_foot_upx(0,0)+_detfootx;
		
		
			// footy location constraints
			if (_bjxx % 2 == 0) //odd
			{
			  _H_q_footy_up.row(0) = _Sfy.row(0);
			  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini); 
			  _F_foot_upy(0,0) = _F_foot_upy(0,0)+_fy - _footy_min; 
		  
			  _H_q_footy_low.row(0) = -_Sfy.row(0);	
			  _F_foot_lowy(0,0) = -_F_foot_upy(0,0) +_detfooty;	
			}
			else
			{
			  _H_q_footy_up.row(0) = _Sfy.row(0);
			  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini); 
			  _F_foot_upy(0,0) = _F_foot_upy(0,0)+_fy + _footy_max; 
			  _H_q_footy_low.row(0) = -_Sfy.row(0);
	// 		  _F_foot_lowy.row(0) = (_H_q_footy_up.row(0) * _V_ini );
			  _F_foot_lowy(0,0) = -_F_foot_upy(0,0)+ _detfooty;
			}	 
	      }
	      else   //two next steps
	      {
			_H_q_footx_up.row(0) = _Sfx.row(0);
			_F_foot_upx.row(0) = -(_H_q_footx_up.row(0) * _V_ini); 
			_F_foot_upx(0,0) = _F_foot_upx(0,0) +_fx + _footx_max;
			_H_q_footx_low.row(0) = -_Sfx.row(0);
	// 		_F_foot_lowx.row(0) = (_H_q_footx_up.row(0) * _V_ini);
			_F_foot_lowx(0,0) = -_F_foot_upx(0,0)+ _detfootx;
		
			// footy location constraints
			if (_bjxx % 2 == 0) //odd
			{
			  _H_q_footy_up.row(0) = _Sfy.row(0);
			  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini); 
			  _F_foot_upy(0,0) = _F_foot_upy(0,0)+_fy - _footy_min;
			  _H_q_footy_low.row(0) = -_Sfy.row(0);
	// 		  _F_foot_lowy.row(0) = (_H_q_footy_up.row(0) * _V_ini);
			  _F_foot_lowy(0,0) = -_F_foot_upy(0,0)+_detfooty;
			}
			else
			{
			  _H_q_footy_up.row(0) = _Sfy.row(0);
			  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini); 
			  _F_foot_upy(0,0) = _F_foot_upy(0,0)+_fy + _footy_max;
			  _H_q_footy_low.row(0) = -_Sfy.row(0);
	// 		  _F_foot_lowy.row(0) = (_H_q_footy_up.row(0) * _V_ini);
			  _F_foot_lowy(0,0) = -_F_foot_upy(0,0)+ _detfooty;
			}
		
			// the next two steps 
			_H_q_footx_up.row(1) = _Sfoot* _Sfx;
			_F_foot_upx.row(1) = -(_H_q_footx_up.row(1) * _V_ini); 	 
			_F_foot_upx(1,0) = _F_foot_upx(1,0) +_footx_max; 
		
			_H_q_footx_low.row(1) = -_Sfoot* _Sfx;
			_F_foot_lowx(1,0) = -_F_foot_upx(1,0)+_detfootx;
		
			// footy location constraints
			if (_bjxx % 2 == 0) //odd
			{
			  _H_q_footy_up.row(1) = _Sfoot* _Sfy;
			  _F_foot_upy.row(1) = -(_H_q_footy_up.row(1) * _V_ini); 
			  _F_foot_upy(1,0) = _F_foot_upy(1,0) + _footy_max; 
		  
			  _H_q_footy_low.row(1) = -_H_q_footy_up.row(1);
	// 		  _F_foot_lowy.row(1) = (_H_q_footy_up.row(1) * _V_ini);
			  _F_foot_lowy(1,0) = -_F_foot_upy(1,0)+_detfooty;		  
			}
			else
			{
			  _H_q_footy_up.row(1) = _Sfoot* _Sfy;
			  _F_foot_upy.row(1) = -(_H_q_footy_up.row(1) * _V_ini); 
			  _F_foot_upy(1,0) = _F_foot_upy(1,0) - _footy_min; 
		  
			  _H_q_footy_low.row(1) = -_H_q_footy_up.row(1);
	// 		  _F_foot_lowy.row(1) = (_H_q_footy_up.row(1) * _V_ini);
			  _F_foot_lowy(1,0) = -_F_foot_upy(1,0)+_detfooty;		  
			}	      	     	      
	      }	      
	      
	      //swing foot veloctiy boundary
	      if (i ==1)
	      {
			_footubxv = -(_Sfx.row(0) * _V_ini - _footx_real_next.row(i+_nT-2));
			_footubxv(0,0) = _footubxv(0,0)  + _footx_max;
		
	// 		_footlbxv = (_Sfx.row(0) * _V_ini - _footx_real_next.row(i+_nT-2));
			_footlbxv(0,0) = -_footubxv(0,0)+ _detfootx;		
		
			_footubyv = -(_Sfy.row(0) * _V_ini - _footy_real_next.row(i+_nT-2));
			_footubyv(0,0) = _footubyv(0,0) +_footy_max;		
		
	// 		_footlbyv = (_Sfy.row(0) * _V_ini - _footy_real_next.row(i+_nT-2));
			_footlbyv(0,0) = -_footubyv(0,0)+_detfooty;		
		
	      }
	      else
	      {
			if (fabs(i*_dt - _tx(_bjxx-1))<=0.01)
			{		
			  _Footvx_max.setZero();  _Footvx_min.setZero();  _Footvy_max.setZero();  _Footvy_min.setZero();		  
			  _footubxv.setZero(); _footlbxv.setZero(); _footubyv.setZero(); _footlbyv.setZero();			
			}
			else
			{
			  _footubxv = -(_Sfx.row(0) * _V_ini - _footx_real_next.row(i+_nT-2));
			  _footubxv(0,0) = _footubxv(0,0) + _footx_vmax*_dt;
	// 		  _footlbxv = (_Sfx.row(0) * _V_ini - _footx_real_next.row(i+_nT-2));
			  _footlbxv(0,0) = -_footubxv(0,0) +_detfootvx*_dt;		  
			  _footubyv = -(_Sfy.row(0) * _V_ini - _footy_real_next.row(i+_nT-2));
			  _footubyv(0,0) = _footubyv(0,0) + _footy_vmax*_dt;		  
	// 		  _footlbyv = (_Sfy.row(0) * _V_ini - _footy_real_next.row(i+_nT-2));	
			  _footlbyv(0,0) = -_footubyv(0,0)+_detfootvy*_dt;		  
			}
	      }

	    ///////////// equality equation	    
	    //equality constraints
	    _F_footz = _Sfz.row(0)*_V_ini - _Lz_ref.row(0);	    
	    _hhhx = _h_h*_V_ini + _pps * _zk.col(i-1) - _Hcom;	

	    _a_hxx = _a_hx * _V_ini + _pps * _thetaxk.col(i-1);
	    _a_hyy = _a_hy * _V_ini + _pps * _thetayk.col(i-1);		    	      
  //////////////////////////////===========////////////////////////////////////////////////////	    
  	    // quadratic program GetSolution
	      
	      if (_method_flag ==0)
	      {
			solve_reactive_step(); ///merely the reactive steps
	      }
	      else
	      {
			if (_method_flag ==1)
			{
			  solve_reactive_step_body_inclination(); /// the reactive steps + body inclination
			}
			else
			{
			  solve_reactive_step_body_inclination_CoMz(); /// the reactive steps + body inclination + height variance
			}
	      }

	    }
    
    
//  ////////////////////////////////results postprocessed////////////
//  /////////===================================================%%%%%	  
//	  // results postprocessed:	  
//	    if (_n_vis == 1)
//	    {
//	      _V_optimal.block< 5*_nh, 1>(0, i-1)  = _V_ini.topRows(5*_nh);
//	      _V_optimal.block< 2, 1>(5*_nh, i-1)  = _V_inix.setConstant(_V_ini(5*_nh,0));
//	      _V_optimal.block< 2, 1>(5*_nh+2, i-1)  = _V_inix.setConstant(_V_ini(5*_nh+1,0));
//	      _V_optimal.block< 2, 1>(5*_nh+4, i-1)  = _V_inix.setConstant(_V_ini(5*_nh+2,0));	      
//	    }
//	    else
//	    {
//	      _V_optimal.col(i-1) = _V_ini;
//	    }
//	    
//	    //next step location
//	    _x_vacc_k.col(i-1) = _V_ini.row(0);
//	    _footx_real.row(_bjxx) = _V_ini.row(5*_nh);
//	    _footx_real_next.row(i+_nT -1) = _V_ini.row(5*_nh);
//	    _xk.col(i) = _a * _xk.col(i-1)  + _b* _x_vacc_k.col(i-1);
//	    _comx(0,i)=_xk(0,i); _comvx(0,i) = _xk(1,i); _comax(0,i)=_xk(2,i); 	  
//	    
//	    _y_vacc_k.col(i-1) = _V_ini.row(0+_nh);
//	    _footy_real.row(_bjxx) = _V_ini.row(5*_nh + _nstep);
//	    _footy_real_next.row(i+_nT -1) = _V_ini.row(5*_nh + _nstep);
//	    _yk.col(i) = _a * _yk.col(i-1)  + _b* _y_vacc_k.col(i-1);
//	    _comy(0,i)=_yk(0,i); _comvy(0,i) = _yk(1,i); _comay(0,i)=_yk(2,i); 
//	    
//	    
//	    
//	    _comx(0,i) +=  _fxx_global;
//	    _comy(0,i) +=  _fyy_global;	  
//	    _footx_real(_bjxx) = _footx_real(_bjxx) + _fxx_global;
//	    _footy_real(_bjxx) = _footy_real(_bjxx) + _fyy_global;	  	  
//	    _footx_real_next1(i+_nT -1) = _footx_real_next(i+_nT -1) + _fxx_global;	  
//	    _footy_real_next1(i+_nT -1) = _footy_real_next(i+_nT -1) + _fyy_global;
//
//	    
//	    
//	    _z_vacc_k.col(i-1) = _V_ini.row(0+2*_nh);
//	    _footz_real.row(_bjxx) = _V_ini.row(5*_nh + 2*_nstep);
//	    _footz_real_next.row(i+_nT -1) = _V_ini.row(5*_nh + 2*_nstep);	  
//	    _zk.col(i) = _a * _zk.col(i-1)  + _b* _z_vacc_k.col(i-1);
//	    _comz(0,i)=_zk(0,i); _comvz(0,i) = _zk(1,i); _comaz(0,i)=_zk(2,i); 
//
//	    _thetax_vacc_k.col(i-1) = _V_ini.row(0+3*_nh);
//	    _thetaxk.col(i) = _a * _thetaxk.col(i-1)  + _b* _thetax_vacc_k.col(i-1);
//	    _thetax(0,i) = _thetaxk(0,i); _thetavx(0,i) = _thetaxk(1,i); _thetaax(0,i)=_thetaxk(2,i); 	
//
//
//	    _thetay_vacc_k.col(i-1) = _V_ini.row(0+4*_nh);
//	    _thetayk.col(i) = _a * _thetayk.col(i-1)  + _b* _thetay_vacc_k.col(i-1);
//	    _thetay(0,i) = _thetayk(0,i); _thetavy(0,i) = _thetayk(1,i); _thetaay(0,i)=_thetayk(2,i); 	
//	    
//	    
//	    // reference relative state    
//	    /// /// relative state to the actual foot lcoation: very good
//	    if (_bjxx % 2 == 0)  // odd : left support
//	    {
//	      estimated_state(0,0) =  estimated_state(0,0) - _Lfoot_location_feedback(0);
//	      estimated_state(3,0) =  estimated_state(3,0) - _Lfoot_location_feedback(1);	 	    
//	    }
//	    else
//	    {
//	      estimated_state(0,0) =  estimated_state(0,0) - _Rfoot_location_feedback(0);
//	      estimated_state(3,0) =  estimated_state(3,0) - _Rfoot_location_feedback(1);	  
//	    }
//  //////============================================================================================================================	      
//  ////////////////////////////// state modified:====================================================================================
//
//	    if (_method_flag <2)
//	    {
//	      estimated_state(6,0) =  _comz(0,i-1);  
//	      estimated_state(7,0) =  _comvz(0,i-1);  
//	      estimated_state(8,0) =  _comaz(0,i-1);
//	    }	  
//	  	    
//	    if (_method_flag <=0)
//	    {
//	      estimated_state(9,0) =  _thetax(0,i-1);  
//	      estimated_state(10,0) =  _thetavx(0,i-1);  
//	      estimated_state(11,0) =  _thetaax(0,i-1);	  
//	      estimated_state(12,0) =  _thetay(0,i-1);  
//	      estimated_state(13,0) =  _thetavy(0,i-1);  
//	      estimated_state(14,0) =  _thetaay(0,i-1);	    
//	    }
//	
//  /*	  ///////////////================state feedback: determined by ratio parameter: lamda==============================////
//	    /// model0================COMX+COMY feedback
//
//	    _xk(0,i) = (estimated_state(0,0)+2*_xk(0,i))/3;             
//	    _xk(1,i) = (estimated_state(1,0)+2*_xk(1,i))/3; 
//	    _xk(2,i) = (estimated_state(2,0)+2*_xk(2,i))/3;
//	    
//
//	    if (i<(int)boost::math::iround(2*_ts(1)/_dt))
//	    {
//		
//	    _yk(0,i) = (estimated_state(3,0)+24*_yk(0,i))/25; _yk(1,i) = (estimated_state(4,0)+24*_yk(1,i))/25;
//	    _yk(2,i) = (estimated_state(5,0)+24*_yk(2,i))/25;		    
//	      
//	    }
//	    else
//	    {
//		
//	    _yk(0,i) = (estimated_state(3,0)+2*_yk(0,i))/3; _yk(1,i) = (estimated_state(4,0)+2*_yk(1,i))/3;
//	    _yk(2,i) = (estimated_state(5,0)+2*_yk(2,i))/3;		    
//	    }
//
//	    ///model1===================COMX+COMY + thetax+ thetay feedback
//	    _thetaxk(0,i) = (estimated_state(9,0)+1*_thetaxk(0,i))/2;
//	    _thetaxk(1,i) = (estimated_state(10,0)+1*_thetaxk(1,i))/2; 
//	    _thetaxk(2,i) = (estimated_state(11,0)+1*_thetaxk(2,i))/2;	  
//	    _thetayk(0,i) = (estimated_state(12,0)+1*_thetayk(0,i))/2; 
//	    _thetayk(1,i) = (estimated_state(13,0)+1*_thetayk(1,i))/2; 
//	    _thetayk(2,i) = (estimated_state(14,0)+1*_thetayk(2,i))/2;	
//	    
//	    //model2====================COMX+COMY + thetax+ thetay + CoMz feedback
//	    _zk(0,i) = (estimated_state(6,0)+2*_zk(0,i))/3; 
//	    _zk(1,i) = (estimated_state(7,0)+2*_zk(1,i))/3; 
//	    _zk(2,i) = (estimated_state(8,0)+2*_zk(2,i))/3;*/	
//	    
//	    
//  ////////////////===============================================================================================	  
//	  /// next two sample time:	actually the preictive value is not reliable  
//// 	    _comx(0,i+1) = _comx(0,i) + _dt * _comvx(0,i); 	  	  
//// 	    _comy(0,i+1) = _comy(0,i) + _dt * _comvy(0,i); 	 
//// 	    _comz(0,i+1) = _comz(0,i) + _dt * _comvz(0,i); 	 
//// 	    _thetax(0,i+1) = _thetax(0,i)+ _dt * _thetavx(0,i); 	
//// 	    _thetay(0,i+1) = _thetay(0,i)+ _dt * _thetavy(0,i); 	
//// 	    
//	    
//	    _torquex_real.col(i) = _j_ini * _thetaax.col(i);
//	    _torquey_real.col(i) = _j_ini * _thetaay.col(i);
//	    
//	    _zmpx_real(0,i) = _comx(0,i) - (_comz(0,i) - _Zsc(i))/(_comaz(0,i)+_ggg(0))*_comax(0,i) - _j_ini * _thetavy(0,i)/(_mass * (_ggg(0) + _comaz(0,i)));
//	    _zmpy_real(0,i) = _comy(0,i) - (_comz(0,i) - _Zsc(i))/(_comaz(0,i)+_ggg(0))*_comay(0,i) + _j_ini * _thetavx(0,i)/(_mass * (_ggg(0) + _comaz(0,i)));
//	    
//
//	    
//	    _footxyz_real.row(0) = _footx_real.transpose();
//	    _footxyz_real.row(1) = _footy_real.transpose();	  
//	    _footxyz_real.row(2) = _footz_real.transpose();
//
//	    
//	    
//	    
//	    /////  generate the trajectory during the double support phase'
//	    _nTdx = (int)boost::math::iround(_td(1)/_dt)+2;
//      // 	    cout <<"nTdx:"<<_nTdx<<endl;
//	    for (int jxx=2; jxx <=_nTdx; jxx++)
//	    {
//	      _x_vacc_k.col(i+jxx-2) = _V_ini.row(jxx-1);
//	      _xk.col(i+jxx-1) = _a * _xk.col(i+jxx-2)  + _b* _x_vacc_k.col(i+jxx-2);
//	      _comx(0,i+jxx-1)=_xk(0,i+jxx-1); _comvx(0,i+jxx-1) = _xk(1,i+jxx-1); _comax(0,i+jxx-1)=_xk(2,i+jxx-1); 	  
//	      
//	      _y_vacc_k.col(i+jxx-2) = _V_ini.row(jxx-1+_nh);
//	      _yk.col(i+jxx-1) = _a * _yk.col(i+jxx-2)  + _b* _y_vacc_k.col(i+jxx-2);
//	      _comy(0,i+jxx-1)=_yk(0,i+jxx-1); _comvy(0,i+jxx-1) = _yk(1,i+jxx-1); _comay(0,i+jxx-1)=_yk(2,i+jxx-1); 	     	      
//	      
//	      _comx(0,i+jxx-1) +=  _fxx_global;
//	      _comy(0,i+jxx-1) +=  _fyy_global;	  	      
//	      
//	      _z_vacc_k.col(i+jxx-2) = _V_ini.row(jxx-1+2*_nh);	  
//	      _zk.col(i+jxx-1) = _a * _zk.col(i+jxx-2)  + _b* _z_vacc_k.col(i+jxx-2);
//	      _comz(0,i+jxx-1)=_zk(0,i+jxx-1); _comvz(0,i+jxx-1) = _zk(1,i+jxx-1); _comaz(0,i+jxx-1)=_zk(2,i+jxx-1); 
//
//	      _thetax_vacc_k.col(i+jxx-2) = _V_ini.row(jxx-1+3*_nh);
//	      _thetaxk.col(i+jxx-1) = _a * _thetaxk.col(i+jxx-2)  + _b* _thetax_vacc_k.col(i+jxx-2);
//	      _thetax(0,i+jxx-1) = _thetaxk(0,i+jxx-1); _thetavx(0,i+jxx-1) = _thetaxk(1,i+jxx-1); _thetaax(0,i+jxx-1)=_thetaxk(2,i+jxx-1); 	
//
//
//	      _thetay_vacc_k.col(i+jxx-2) = _V_ini.row(jxx-1+4*_nh);
//	      _thetayk.col(i+jxx-1) = _a * _thetayk.col(i+jxx-2)  + _b* _thetay_vacc_k.col(i+jxx-2);
//	      _thetay(0,i+jxx-1) = _thetayk(0,i+jxx-1); _thetavy(0,i+jxx-1) = _thetayk(1,i+jxx-1); _thetaay(0,i+jxx-1)=_thetayk(2,i+jxx-1); 		      
//	    
//	      _torquex_real.col(i+jxx-1) = _j_ini * _thetaax.col(i+jxx-1);
//	      _torquey_real.col(i+jxx-1) = _j_ini * _thetaay.col(i+jxx-1);
//	      
//	      _zmpx_real(0,i+jxx-1) = _comx(0,i+jxx-1) - (_comz(0,i+jxx-1) - _Zsc(i+jxx-1))/(_comaz(0,i+jxx-1)+_ggg(0))*_comax(0,i+jxx-1) - _j_ini * _thetaay(0,i+jxx-1)/(_mass * (_ggg(0) + _comaz(0,i+jxx-1)));
//	      _zmpy_real(0,i+jxx-1) = _comy(0,i+jxx-1) - (_comz(0,i+jxx-1) - _Zsc(i+jxx-1))/(_comaz(0,i+jxx-1)+_ggg(0))*_comay(0,i+jxx-1) + _j_ini * _thetaax(0,i+jxx-1)/(_mass * (_ggg(0) + _comaz(0,i+jxx-1)));
//	      
//	    }	    
//	    
//	    _comx(0,i+1) = _comx(0,i) + _dt * _comvx(0,i); 	  	  
//	    _comy(0,i+1) = _comy(0,i) + _dt * _comvy(0,i); 	 
//	    _comz(0,i+1) = _comz(0,i) + _dt * _comvz(0,i); 	 
//	    _thetax(0,i+1) = _thetax(0,i)+ _dt * _thetavx(0,i); 	
//	    _thetay(0,i+1) = _thetay(0,i)+ _dt * _thetavy(0,i); 	
//	
//	  
//	  if (i>=1)
//	  {	  
//	    _Rfootx(0) = _Rfootx(1);
//	    _Lfootx(0) = _Lfootx(1);
//	    _Rfooty(0) = _Rfooty(1);
//	    _Lfooty(0) = _Lfooty(1);
//	    _comx(0) = _comx(1);	
//	    _comy(0) = _comy(1);
//	    _comz(0) = _comz(1);	  
//	  }	 
// 
//	    _tcpu(0,i-1) = 0; 
// 
// 
      }
    else
    {
	 if(i <= _n_end_walking+(int)boost::math::iround(_tstep/_dt/2)){
	   
 	    Indexfind(i*_dt,xyz1);                   //// step cycle number when (i)*dt fall into : current sampling time
 	    _bjxx = _j_period+1;  //coincidence with matlab 
 	    _j_period = 0;		   
 	   
 	    _t_f.setLinSpaced(_nh,(i+1)*_dt, (i+_nh)*_dt);
 	    
 	    Indexfind(_t_f(0),xyz1);                /// step cycle number when (i+1)*dt fall into : current sampling time
 	    _bjx1 = _j_period+1;
 	    _j_period = 0;
 	    
 	   Eigen::Matrix<double, 4, 1> _comy_temp;
 	   _comy_temp.setZero();
 	   _comy_temp(0) = _comy(0,_n_end_walking-2);
 	   _comy_temp(1) = _comy(0,_n_end_walking-1);
 	   _comy_temp(2) = 0;
 	   _comy_temp(3) = 0;
 
 	   
 	   Eigen::Matrix<double, 4, 4> _comy_matrix;
 	   
 	   double ix_temp1 = -1.0;
 	   double ix_temp2 = 0.0;	   
 	   double ix_temp3 = (int)boost::math::iround(_tstep/_dt)/2+1.0;
 	   
 	  
 	  Eigen::Matrix4d AAA_inv;
 	  
 	  double abx1, abx2, abx3, abx4;
 	  abx1 = ((ix_temp1 - ix_temp2)*pow(ix_temp1 - ix_temp3, 2));
 	  abx2 = ((ix_temp1 - ix_temp2)*pow(ix_temp2 - ix_temp3, 2));
 	  abx3 =(pow(ix_temp1 - ix_temp3, 2)*pow(ix_temp2 - ix_temp3, 2));
 	  abx4 = ((ix_temp1 - ix_temp3)*(ix_temp2 - ix_temp3));
 	  
 
 	  AAA_inv(0,0) = 1/ abx1;
 	  AAA_inv(0,1) =  -1/ abx2;
 	  AAA_inv(0,2) = (ix_temp1 + ix_temp2 - 2*ix_temp3)/ abx3;
 	  AAA_inv(0,3) = 1/ abx4;
 	  
 	  AAA_inv(1,0) = -(ix_temp2 + 2*ix_temp3)/ abx1;
 	  AAA_inv(1,1) = (ix_temp1 + 2*ix_temp3)/ abx2;
 	  AAA_inv(1,2) = -(pow(ix_temp1, 2) + ix_temp1*ix_temp2 + pow(ix_temp2, 2) - 3*pow(ix_temp3, 2))/ abx3;
 	  AAA_inv(1,3) = -(ix_temp1 + ix_temp2 + ix_temp3)/ abx4;
 	  
 	  AAA_inv(2,0) = (ix_temp3*(2*ix_temp2 + ix_temp3))/ abx1;
 	  AAA_inv(2,1) = -(ix_temp3*(2*ix_temp1 + ix_temp3))/ abx2;
 	  AAA_inv(2,2) = (ix_temp3*(2*pow(ix_temp1, 2) + 2*ix_temp1*ix_temp2 - 3*ix_temp3*ix_temp1 + 2*pow(ix_temp2, 2) - 3*ix_temp3*ix_temp2))/ abx3;
 	  AAA_inv(2,3) = (ix_temp1*ix_temp2 + ix_temp1*ix_temp3 + ix_temp2*ix_temp3)/ abx4;
 	  
 	  AAA_inv(3,0) = -(ix_temp2*pow(ix_temp3, 2))/ abx1;
 	  AAA_inv(3,1) = (ix_temp1*pow(ix_temp3, 2))/ abx2;
 	  AAA_inv(3,2) = (ix_temp1*ix_temp2*(ix_temp1*ix_temp2 - 2*ix_temp1*ix_temp3 - 2*ix_temp2*ix_temp3 + 3*pow(ix_temp3, 2)))/ abx3;
 	  AAA_inv(3,3) = -(ix_temp1*ix_temp2*ix_temp3)/ abx4;
 	   
 	   _comy_matrix_inv = AAA_inv * _comy_temp;
 
 
 	    /////  generate the trajectory during the double support phase'
 	    _nTdx = (int)boost::math::iround(_ts(1)/_dt);
 // 	    cout <<"nTdx:"<<_nTdx<<endl;
 	    for (int jxx=1; jxx <=_nTdx; jxx++)
 	    {
 	      
 	      if (i+jxx-1 <= _n_end_walking+(int)boost::math::iround(_tstep/_dt/2))
 	      {
 		Eigen::Matrix<double, 1, 4> _t_temp;
 		
 		_t_temp(0) = pow((double)(i+jxx-1-_n_end_walking+1), 3);
 		_t_temp(1) = pow((double)(i+jxx-1-_n_end_walking+1), 2);
 		_t_temp(2) = pow((double)(i+jxx-1-_n_end_walking+1), 1);
 		_t_temp(3) = pow((double)(i+jxx-1-_n_end_walking+1), 0);	
 		
 		_comy.col(i+jxx-1) = _t_temp* _comy_matrix_inv;
 		
 		Eigen::Matrix<double, 1, 4> _t_tempv;
 		
 		_t_tempv(0) = 3*pow((double)(i+jxx-1-_n_end_walking+1), 2);
 		_t_tempv(1) = 2*pow((double)(i+jxx-1-_n_end_walking+1), 1);
 		_t_tempv(2) = 1;
 		_t_tempv(3) = 0;	
 		
 		_comvy.col(i+jxx-1) = _t_tempv* _comy_matrix_inv;
 
 		Eigen::Matrix<double, 1, 4> _t_tempa;
 		
 		_t_tempa(0) = 6*pow((double)(i+jxx-1-_n_end_walking+1), 1);
 		_t_tempa(1) = 2;
 		_t_tempa(2) = 0;
 		_t_tempa(3) = 0;	
 		
 		_comay.col(i+jxx-1) = _t_tempa* _comy_matrix_inv;				
 	      }
 	      else
 	      {
                _comy(0,i+jxx-1)=_comy(0,i+jxx-2); _comvy(0,i+jxx-1) = 0; _comay(0,i+jxx-1)= 0;		
 	      }
 	      
 	      _comx(0,i+jxx-1)=_comx(0,_n_end_walking-1); _comvx(0,i+jxx-1) = 0; _comax(0,i+jxx-1)= 0; 	  
 	       	     	      
 	      _comz(0,i+jxx-1)=_comz(0,_n_end_walking-1); _comvz(0,i+jxx-1) = 0; _comaz(0,i+jxx-1)= 0; 
 
 	      _thetax(0,i+jxx-1) = _thetax(0,_n_end_walking-1); _thetavx(0,i+jxx-1) = 0; _thetaax(0,i+jxx-1)= 0; 	
 
 	      _thetay(0,i+jxx-1) = _thetay(0,_n_end_walking-1); _thetavy(0,i+jxx-1) = 0; _thetaay(0,i+jxx-1)= 0; 		      
 	   
 	      _torquex_real.col(i+jxx-1) = _j_ini * _thetaax.col(i+jxx-1);
 	      _torquey_real.col(i+jxx-1) = _j_ini * _thetaay.col(i+jxx-1);
 	      
 	      _zmpx_real(0,i+jxx-1) = _comx(0,i+jxx-1) - (_comz(0,i+jxx-1) - _Zsc(i+jxx-1))/(_comaz(0,i+jxx-1)+_ggg(0))*_comax(0,i+jxx-1) - _j_ini * _thetaay(0,i+jxx-1)/(_mass * (_ggg(0) + _comaz(0,i+jxx-1)));
 	      _zmpy_real(0,i+jxx-1) = _comy(0,i+jxx-1) - (_comz(0,i+jxx-1) - _Zsc(i+jxx-1))/(_comaz(0,i+jxx-1)+_ggg(0))*_comay(0,i+jxx-1) + _j_ini * _thetaax(0,i+jxx-1)/(_mass * (_ggg(0) + _comaz(0,i+jxx-1)));
 	      
 	    }	    
 	   
 	  _footx_real_next.row(i+_nT -1)=_footx_real_next.row(_n_end_walking-1+_nT -1)	;
 	  _footy_real_next.row(i+_nT -1)=_footy_real_next.row(_n_end_walking-1-_nT -1)	;	   
 	   
 	   for (int jxxx = _bjxx+1; jxxx<_footstepsnumber; jxxx++){
 	   	    _footx_real(jxxx) = _footx_real(_bjxx) ;
 	   	    _footy_real(jxxx) = _footy_real(_bjxx-1);		     
 	     
 	  }
 	    _footxyz_real.row(0) = _footx_real.transpose();
 	    _footxyz_real.row(1) = _footy_real.transpose();	  
 	    _footxyz_real.row(2) = _footz_real.transpose();
 // 	    cout<<"_footy_real:"<<_footy_real<<endl;
  
	   
	}
	 else
	{
	   
	    Indexfind(i*_dt,xyz1);                   //// step cycle number when (i)*dt fall into : current sampling time
	    _bjxx = _j_period+1;  //coincidence with matlab 
	    _j_period = 0;	
	    
	    _t_f.setLinSpaced(_nh,(i+1)*_dt, (i+_nh)*_dt);
	    
	    Indexfind(_t_f(0),xyz1);                /// step cycle number when (i+1)*dt fall into : current sampling time
	    _bjx1 = _j_period+1;
	    _j_period = 0;	
	    	      
	   _comy(0,i) = _comy(0,i-1);
	   _comy(0,i+1) = _comy(0,i);	    
	   
	   _comx(0,i) = _comx(0,_n_end_walking-1);
	   _comx(0,i+1) = _comx(0,_n_end_walking-1);
	   _comz(0,i) = _comz(0,_n_end_walking-1);
	   _comz(0,i+1) = _comz(0,_n_end_walking-1);	   
	   _thetax(0,i) = _thetax(0,_n_end_walking-1);
	   _thetax(0,i+1) = _thetax(0,_n_end_walking-1);	   
	   _thetay(0,i) = _thetay(0,_n_end_walking-1);
	   _thetay(0,i+1) = _thetay(0,_n_end_walking-1);
	   	   
	    /////  generate the trajectory during the double support phase'
	    _nTdx = (int)boost::math::iround(_ts(1)/_dt)+2;
	    for (int jxx=2; jxx <=_nTdx; jxx++)
	    {
	      _comx(0,i+jxx-1)=_comx(0,i+jxx-2); _comvx(0,i+jxx-1) = 0; _comax(0,i+jxx-1)= 0; 	  

	      _comy(0,i+jxx-1)=_comy(0,i+jxx-2); _comvy(0,i+jxx-1) = 0; _comay(0,i+jxx-1)= 0; 	     	      

	      _comz(0,i+jxx-1)=_comz(0,i+jxx-2); _comvz(0,i+jxx-1) = 0; _comaz(0,i+jxx-1)= 0; 

	      _thetax(0,i+jxx-1) = _thetax(0,i+jxx-2); _thetavx(0,i+jxx-1) = 0; _thetaax(0,i+jxx-1)= 0; 	

	      _thetay(0,i+jxx-1) = _thetay(0,i+jxx-2); _thetavy(0,i+jxx-1) = 0; _thetaay(0,i+jxx-1)= 0; 		      
	   
	      _torquex_real.col(i+jxx-1) = _j_ini * _thetaax.col(i+jxx-1);
	      _torquey_real.col(i+jxx-1) = _j_ini * _thetaay.col(i+jxx-1);
	      
	      _zmpx_real(0,i+jxx-1) = _comx(0,i+jxx-1) - (_comz(0,i+jxx-1) - _Zsc(i+jxx-1))/(_comaz(0,i+jxx-1)+_ggg(0))*_comax(0,i+jxx-1) - _j_ini * _thetaay(0,i+jxx-1)/(_mass * (_ggg(0) + _comaz(0,i+jxx-1)));
	      _zmpy_real(0,i+jxx-1) = _comy(0,i+jxx-1) - (_comz(0,i+jxx-1) - _Zsc(i+jxx-1))/(_comaz(0,i+jxx-1)+_ggg(0))*_comay(0,i+jxx-1) + _j_ini * _thetaax(0,i+jxx-1)/(_mass * (_ggg(0) + _comaz(0,i+jxx-1)));
	      
	    }	    
	    	   	   
		  _footx_real_next.row(i+_nT -1)=_footx_real_next.row(_n_end_walking-1+_nT -1)	;
		  _footy_real_next.row(i+_nT -1)=_footy_real_next.row(_n_end_walking-1-_nT -1)	;	   
	   
		  for (int jxxx = _bjxx+1; jxxx<_footstepsnumber; jxxx++){
			  _footx_real(jxxx) = _footx_real(_bjxx) ;
			  _footy_real(jxxx) = _footy_real(_bjxx-1);		     
	     
		  }
	    _footxyz_real.row(0) = _footx_real.transpose();
	    _footxyz_real.row(1) = _footy_real.transpose();	  
	    _footxyz_real.row(2) = _footz_real.transpose(); 
	} 
	}
}


//////////////////////////// modified
void MPCClass::Indexfind(double goalvari, int xyz)
{
        _j_period = 0;
	if (xyz<0.05)
	{
	  while (goalvari >= _tx(_j_period))
	  {
	    _j_period++;
	  }
	  
	  _j_period = _j_period-1;	  
	}
	else
	{
	  while ( fabs(goalvari - _t_f(_j_period)) >0.0001 )
	  {
	    _j_period++;
	  }
	  	  
	}	  
}

////////////////////////////modified
int MPCClass::is_sigular(int num)
{
	if ( num & 1)
	{	  
	  return 1; //sigular 
	}
	else
	{
	  return 0;	// odd  	 
	}	   
}

///// only works once when initialize
Eigen::MatrixXd  MPCClass::Matrix_ps(Eigen::Matrix<double,3,3> a, int nh,Eigen::RowVector3d cxps)
{
//   Eigen::MatrixXd matrixps(nh,3);
  Eigen::MatrixXd matrixps;
  matrixps.setZero(nh,3);  
  
  
  
  Eigen::MatrixXd A;
//   A.setIdentity(a.rows(),a.cols());
  
  for (int i = 0; i < nh; i++) {
    A.setIdentity(a.rows(),a.cols());
    for (int j = 1; j < i+2; j++)
    {
      A = A*a;
    }  
    
     matrixps.middleRows(i, 1)= cxps * A;      
  }
    
  return matrixps;
}

//
Eigen::MatrixXd MPCClass::Matrix_pu(Eigen::Matrix<double,3,3> a, Eigen::Matrix<double,3,1> b, int nh, Eigen::RowVector3d cxpu)
{
  Eigen::MatrixXd matrixpu;
  matrixpu.setZero(nh,nh);
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd Tempxx;
  
  
  for (int i = 1; i < nh+1; i++) {
    for (int j = 1; j < i+1; j++)
    { 
      A.setIdentity(a.rows(),a.cols());      
      if (j==i)
      {
	Tempxx = cxpu * A * b;
	matrixpu(i-1,j-1) = Tempxx(0,0);
      }
      else
      {	
	for (int k = 1; k < i-j+1; k++)
	{
	  A = A*a;
	}
	Tempxx = cxpu * A * b;
	matrixpu(i-1,j-1) = Tempxx(0,0);
      }          
    }       
  }
    
  return matrixpu;  
}


void MPCClass::Matrix_large(int i_f)
{
  if (i_f == 1)
  {
    for(int jxx=1; jxx<=_Nt; jxx++)
    {
      _matrix_large1.row(jxx-1) = ZMPx_constraints_halfyyy1.row(jxx-1) * ZMPx_constraints_halfxxx1;
      
      _matrix_large2.row(jxx-1) = ZMPx_constraints_halfyyy2.row(jxx-1) * ZMPx_constraints_halfxxx2;
    }
         
  }
  else
  {
    if (i_f ==2)
    {    
      for(int jxx=1; jxx<=_Nt; jxx++)
      {      
	_matrix_large1.row(jxx-1) = ZMPx_constraints_halfyyy3.row(jxx-1) * ZMPx_constraints_halfxxx1;
	
	_matrix_large2.row(jxx-1) = ZMPx_constraints_halfyyy4.row(jxx-1) * ZMPx_constraints_halfxxx2;  
      }
    }
    else
    {             
        if (i_f ==3)
        {    
            for(int jxx=1; jxx<=_Nt; jxx++)
            {      
                
 //           _phi_i_x_up1 = ZMPx_constraints_offfline[jxx-1] + _ZMPx_constraints_half2 * ZMPx_constraints_half[jxx-1];   
//            _phi_i_y_up1 = ZMPy_constraints_offfline[jxx-1] + _ZMPy_constraints_half2 * ZMPy_constraints_half[jxx-1];                 
                 _matrix_large1.row(jxx-1) = _ZMPx_constraints_half2.row(jxx-1)*_ZMPx_constraints_half_va;
                 _matrix_large2.row(jxx-1) = _ZMPy_constraints_half2.row(jxx-1)*_ZMPy_constraints_half_va;
 
            }
        }
    }  
  }
}






///DATA SAVING:modified=========================================================
void MPCClass::File_wl()
{
        
// 	CoMMM_ZMP_foot.setZero();
	CoMMM_ZMP_foot.block<1,_nsum>(0, 0) = _comx;
	CoMMM_ZMP_foot.block<1,_nsum>(1, 0) = _comy;	
	CoMMM_ZMP_foot.block<1,_nsum>(2, 0) = _comz;	
	CoMMM_ZMP_foot.block<1,_nsum>(3, 0) = _zmpx_real;	
	CoMMM_ZMP_foot.block<1,_nsum>(4, 0) = _zmpy_real;
	CoMMM_ZMP_foot.block<1,_nsum>(5, 0) = _thetax;	
	CoMMM_ZMP_foot.block<1,_nsum>(6, 0) = _thetay;	
	CoMMM_ZMP_foot.block<1,_nsum>(7, 0) = _torquex_real;
	CoMMM_ZMP_foot.block<1,_nsum>(8, 0) = _torquey_real;
	CoMMM_ZMP_foot.block<1,_nsum>(9, 0) = _footx_real_next1.transpose();	
	CoMMM_ZMP_foot.block<1,_nsum>(10, 0) = _footy_real_next1.transpose();	
	CoMMM_ZMP_foot.block<1,_nsum>(11, 0) = _footz_real_next.transpose();
	
	CoMMM_ZMP_foot.block<1,_nsum>(12, 0) = _Lfootx;	
	CoMMM_ZMP_foot.block<1,_nsum>(13, 0) = _Lfooty;	
	CoMMM_ZMP_foot.block<1,_nsum>(14, 0) = _Lfootz;
	CoMMM_ZMP_foot.block<1,_nsum>(15, 0) = _Rfootx;	
	CoMMM_ZMP_foot.block<1,_nsum>(16, 0) = _Rfooty;	
	CoMMM_ZMP_foot.block<1,_nsum>(17, 0) = _Rfootz;

	CoMMM_ZMP_foot.block<1,_nsum>(18, 0) = _comvx;
	CoMMM_ZMP_foot.block<1,_nsum>(19, 0) = _comax;
	
	CoMMM_ZMP_foot.block<1,_nsum>(20, 0) = _comvy;	
	CoMMM_ZMP_foot.block<1,_nsum>(21, 0) = _comay;
	
	CoMMM_ZMP_foot.block<1,_nsum>(22, 0) = _comvz;	
	CoMMM_ZMP_foot.block<1,_nsum>(23, 0) = _comaz;	

	CoMMM_ZMP_foot.block<1,_nsum>(24, 0) = _thetavx;
	CoMMM_ZMP_foot.block<1,_nsum>(25, 0) = _thetaax;
	
	CoMMM_ZMP_foot.block<1,_nsum>(26, 0) = _thetavy;	
	CoMMM_ZMP_foot.block<1,_nsum>(27, 0) = _thetaay;

	
	
	
	
  
	std::string fileName = "NMPC_runtime.txt" ;
	std::ofstream outfile( fileName.c_str() ) ; // file name and the operation type. 	
        for(int i=0; i<_tcpu.rows(); i++){
           for(int j=0; j<_tcpu.cols(); j++){
                 outfile << (double) _tcpu(i,j) << " " ; 
           }
           outfile << std::endl;       // a   newline
        }
        outfile.close();	


	std::string fileName1 = "NMPC_optimal_trajectory.txt" ;
	std::ofstream outfile1( fileName1.c_str() ) ; // file name and the operation type.        
	
        for(int i=0; i<CoMMM_ZMP_foot.rows(); i++){
           for(int j=0; j<CoMMM_ZMP_foot.cols(); j++){
                 outfile1 << (double) CoMMM_ZMP_foot(i,j) << " " ; 
           }
           outfile1 << std::endl;       // a   newline
        }
        outfile1.close();	
	
	
	
}

/// three model MPC solution :modified================================================================
void MPCClass::solve_reactive_step()
{ 

  _G = _Q_goal1;
  _g0 = _q_goal1;
  _X = _V_ini;
	    
  

  _CI.block<_Nt,_nh>(0,0) = _H_q_upx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,_nh) = _H_q_lowx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,2*_nh) = _H_q_upy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,3*_nh) = _H_q_lowy.transpose() * (-1); 
  _CI.block<_Nt,_nh>(0,4*_nh) = _H_hacc_lowz.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh) = _H_q_footx_up.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh+_nstep) = _H_q_footx_low.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh+2*_nstep) = _H_q_footy_up.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh+3*_nstep) = _H_q_footy_low.transpose() * (-1);  

  
  _CI.block<_Nt,1>(0,5*_nh+4*_nstep) = _Footvx_max.transpose() * (-1);
  _CI.block<_Nt,1>(0,5*_nh+4*_nstep+1) = _Footvx_min.transpose() * (-1);
  _CI.block<_Nt,1>(0,5*_nh+4*_nstep+2) = _Footvy_max.transpose() * (-1);
  _CI.block<_Nt,1>(0,5*_nh+4*_nstep+3) = _Footvy_min.transpose() * (-1);

  

  
  _ci0.block(0, 0,_nh,1) = _F_zmp_upx;
  _ci0.block(_nh, 0,_nh,1) = _F_zmp_lowx;
  _ci0.block(2*_nh, 0,_nh,1) = _F_zmp_upy;
  _ci0.block(3*_nh, 0,_nh,1) = _F_zmp_lowy;
  _ci0.block(4*_nh, 0,_nh,1) = _F_hacc_lowz;
   
  _ci0.block(5*_nh, 0,_nstep,1) = _F_foot_upx;
  _ci0.block(5*_nh+_nstep, 0,_nstep,1) = _F_foot_lowx;
  _ci0.block(5*_nh+2*_nstep, 0,_nstep,1) = _F_foot_upy;
  _ci0.block(5*_nh+3*_nstep, 0,_nstep,1) = _F_foot_lowy;
   
  _ci0.block(5*_nh+4*_nstep, 0,1,1) = _footubxv;
  _ci0.block(5*_nh+4*_nstep+1, 0,1,1) = _footlbxv;
  _ci0.block(5*_nh+4*_nstep+2, 0,1,1) = _footubyv;
  _ci0.block(5*_nh+4*_nstep+3, 0,1,1) = _footlbyv;

  
  _CE.block(0,0, _Nt,1) = _H_q_footz.transpose();
  _CE.block(0,1, _Nt,_nh) = _h_h.transpose();
  _CE.block(0,_nh+1, _Nt,_nh) = _a_hx.transpose();
  _CE.block(0,2*_nh+1, _Nt,_nh) = _a_hy.transpose();
  
  _ce0.block(0,0, 1,1) = _F_footz;
  _ce0.block(1,0, _nh,1) = _hhhx;  
  _ce0.block(1+_nh,0, _nh,1) = _a_hxx;  
  _ce0.block(1+2*_nh,0, _nh,1) = _a_hyy;    
  
  
  Solve();  

}

void MPCClass::solve_reactive_step_body_inclination()
{    

  _G = _Q_goal1;
  _g0 = _q_goal1;
  _X = _V_ini;

  
  _CI.block<_Nt,_nh>(0,0) = _H_q_upx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,_nh) = _H_q_lowx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,2*_nh) = _H_q_upy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,3*_nh) = _H_q_lowy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,4*_nh) = _H_hacc_lowz.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh) = _H_q_footx_up.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh+_nstep) = _H_q_footx_low.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh+2*_nstep) = _H_q_footy_up.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,5*_nh+3*_nstep) = _H_q_footy_low.transpose() * (-1);	    

  _CI.block<_Nt,_nh>(0,5*_nh+4*_nstep) = _q_upx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,6*_nh+4*_nstep) = _q_lowx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,7*_nh+4*_nstep) = _q_upy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,8*_nh+4*_nstep) = _q_lowy.transpose() * (-1);
  
  _CI.block<_Nt,_nh>(0,9*_nh+4*_nstep) = _t_upx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,10*_nh+4*_nstep) = _t_lowx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,11*_nh+4*_nstep) = _t_upy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,12*_nh+4*_nstep) = _t_lowy.transpose() * (-1);
  
  _CI.block<_Nt,1>(0,13*_nh+4*_nstep) = _Footvx_max.transpose() * (-1);
  _CI.block<_Nt,1>(0,13*_nh+4*_nstep+1) = _Footvx_min.transpose() * (-1);
  _CI.block<_Nt,1>(0,13*_nh+4*_nstep+2) = _Footvy_max.transpose() * (-1);
  _CI.block<_Nt,1>(0,13*_nh+4*_nstep+3) = _Footvy_min.transpose() * (-1);

  

  
  _ci0.block(0, 0,_nh,1) = _F_zmp_upx;
  _ci0.block(_nh, 0,_nh,1) = _F_zmp_lowx;
  _ci0.block(2*_nh, 0,_nh,1) = _F_zmp_upy;
  _ci0.block(3*_nh, 0,_nh,1) = _F_zmp_lowy; 
  _ci0.block(4*_nh, 0,_nh,1) = _F_hacc_lowz;
  _ci0.block(5*_nh, 0,_nstep,1) = _F_foot_upx;
  _ci0.block(5*_nh+_nstep, 0,_nstep,1) = _F_foot_lowx;
  _ci0.block(5*_nh+2*_nstep, 0,_nstep,1) = _F_foot_upy;
  _ci0.block(5*_nh+3*_nstep, 0,_nstep,1) = _F_foot_lowy;	    
  _ci0.block(5*_nh+4*_nstep, 0,_nh,1) = _qq_upx;
  _ci0.block(6*_nh+4*_nstep, 0,_nh,1) = _qq_lowx;
  _ci0.block(7*_nh+4*_nstep, 0,_nh,1) = _qq_upy;
  _ci0.block(8*_nh+4*_nstep, 0,_nh,1) = _qq_lowy;
  
  _ci0.block(9*_nh+4*_nstep, 0,_nh,1) = _tt_upx;
  _ci0.block(10*_nh+4*_nstep, 0,_nh,1) = _tt_lowx;
  _ci0.block(11*_nh+4*_nstep, 0,_nh,1) = _tt_upy;
  _ci0.block(12*_nh+4*_nstep, 0,_nh,1) = _tt_lowy;
  
  _ci0.block(13*_nh+4*_nstep, 0,1,1) = _footubxv;
  _ci0.block(13*_nh+4*_nstep+1, 0,1,1) = _footlbxv;
  _ci0.block(13*_nh+4*_nstep+2, 0,1,1) = _footubyv;
  _ci0.block(13*_nh+4*_nstep+3, 0,1,1) = _footlbyv;

  
  _CE.block(0,0, _Nt,1) = _H_q_footz.transpose();
  _CE.block(0,1, _Nt,_nh) = _h_h.transpose();
  
  _ce0.block(0,0, 1,1) = _F_footz;
  _ce0.block(1,0, _nh,1) = _hhhx;  

  Solve();  

}

void MPCClass::solve_reactive_step_body_inclination_CoMz()
{  

  _G = _Q_goal1;
  _g0 = _q_goal1;
  _X = _V_ini;

  
  _CI.block<_Nt,_nh>(0,0) = _H_q_upx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,_nh) = _H_q_lowx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,2*_nh) = _H_q_upy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,3*_nh) = _H_q_lowy.transpose() * (-1);
  
  _CI.block<_Nt,_nh>(0,4*_nh) = _H_h_upz.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,5*_nh) = _H_h_lowz.transpose() * (-1);
  
  _CI.block<_Nt,_nh>(0,6*_nh) = _H_hacc_lowz.transpose() * (-1);
  
  _CI.block<_Nt,_nstep>(0,7*_nh) = _H_q_footx_up.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,7*_nh+_nstep) = _H_q_footx_low.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,7*_nh+2*_nstep) = _H_q_footy_up.transpose() * (-1);
  _CI.block<_Nt,_nstep>(0,7*_nh+3*_nstep) = _H_q_footy_low.transpose() * (-1);	    

  _CI.block<_Nt,_nh>(0,7*_nh+4*_nstep) = _q_upx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,8*_nh+4*_nstep) = _q_lowx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,9*_nh+4*_nstep) = _q_upy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,10*_nh+4*_nstep,_Nt,_nh) = _q_lowy.transpose() * (-1);
  
  _CI.block<_Nt,_nh>(0,11*_nh+4*_nstep) = _t_upx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,12*_nh+4*_nstep) = _t_lowx.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,13*_nh+4*_nstep) = _t_upy.transpose() * (-1);
  _CI.block<_Nt,_nh>(0,14*_nh+4*_nstep) = _t_lowy.transpose() * (-1);
  
  _CI.block<_Nt,1>(0,15*_nh+4*_nstep) = _Footvx_max.transpose() * (-1);
  _CI.block<_Nt,1>(0,15*_nh+4*_nstep+1) = _Footvx_min.transpose() * (-1);
  _CI.block<_Nt,1>(0,15*_nh+4*_nstep+2) = _Footvy_max.transpose() * (-1);
  _CI.block<_Nt,1>(0,15*_nh+4*_nstep+3) = _Footvy_min.transpose() * (-1);

  

  
  _ci0.block(0, 0,_nh,1) = _F_zmp_upx;
  _ci0.block(_nh, 0,_nh,1) = _F_zmp_lowx;
  _ci0.block(2*_nh, 0,_nh,1) = _F_zmp_upy;
  _ci0.block(3*_nh, 0,_nh,1) = _F_zmp_lowy;
  
  _ci0.block(4*_nh, 0,_nh,1) = _F_h_upz;
  _ci0.block(5*_nh, 0,_nh,1) = _F_h_lowz;
  
  _ci0.block(6*_nh, 0,_nh,1) = _F_hacc_lowz;
  
  _ci0.block(7*_nh, 0,_nstep,1) = _F_foot_upx;
  _ci0.block(7*_nh+_nstep, 0,_nstep,1) = _F_foot_lowx;
  _ci0.block(7*_nh+2*_nstep, 0,_nstep,1) = _F_foot_upy;
  _ci0.block(7*_nh+3*_nstep, 0,_nstep,1) = _F_foot_lowy;	    

  _ci0.block(7*_nh+4*_nstep, 0,_nh,1) = _qq_upx;
  _ci0.block(8*_nh+4*_nstep, 0,_nh,1) = _qq_lowx;
  _ci0.block(9*_nh+4*_nstep, 0,_nh,1) = _qq_upy;
  _ci0.block(10*_nh+4*_nstep, 0,_nh,1) = _qq_lowy;
  
  _ci0.block(11*_nh+4*_nstep, 0,_nh,1) = _tt_upx;
  _ci0.block(12*_nh+4*_nstep, 0,_nh,1) = _tt_lowx;
  _ci0.block(13*_nh+4*_nstep, 0,_nh,1) = _tt_upy;
  _ci0.block(14*_nh+4*_nstep, 0,_nh,1) = _tt_lowy;
  
  _ci0.block(15*_nh+4*_nstep, 0,1,1) = _footubxv;
  _ci0.block(15*_nh+4*_nstep+1, 0,1,1) = _footlbxv;
  _ci0.block(15*_nh+4*_nstep+2, 0,1,1) = _footubyv;
  _ci0.block(15*_nh+4*_nstep+3, 0,1,1) = _footlbyv;

  
  
  
  
  
  _CE = _H_q_footz.transpose();
  _ce0 = _F_footz;
  
  
  
   Solve();  

}

void MPCClass::Solve()
{
// min 0.5 * x G x + g0 x
// _s.t.
// 		CE^T x + ce0 = 0
// 		CI^T x + ci0 >= 0
		solveQP();
		if (_X.rows() == _Nt)
		{
		  _V_ini += _X;
		}

}

/////////////////////////////////////////////////////////////////////////////////=================swing foot trajectory==============================================////////////////////////
/////////////============================================================================================================================////////////////////////////////////////////////
//// foot trajectory solve--------polynomial ================================================
void MPCClass::Foot_trajectory_solve(int j_index,bool _stopwalking)
{
  // judge if stop  
  if(_stopwalking)  
  {
    for (int i_t = _bjx1+1; i_t < _footstepsnumber; i_t++) {	  
      _lift_height_ref(i_t) = 0;  
    }	  
  }  

  _footxyz_real(1,0) = -_stepwidth(0);
  
  // foot trajectory generation:
  if (_bjx1 >= 2)
  {
    
    if (_bjx1 % 2 == 0)           //odd:left support
    {
  //     no change on the left support location
	_Lfootx(j_index) = _Lfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Lfooty(j_index) = _Lfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Lfootz(j_index) = _Lfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	
	_Lfootx(j_index+1) = _Lfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Lfooty(j_index+1) = _Lfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Lfootz(j_index+1) = _Lfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);    
    
	/// right swing
	if ((j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double support
	{
	  _Rfootx(j_index) = _Rfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Rfooty(j_index) = _Rfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Rfootz(j_index) = _Rfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);	
	  
	  _Rfootx(j_index+1) = _Rfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Rfooty(j_index+1) = _Rfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Rfootz(j_index+1) = _Rfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);		  
	}
	else
	{
	  double t_des = (j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt) +1)*_dt;
	  Eigen::Vector3d t_plan(0,0,0);
// 	  t_plan(0) = t_des - _dt;
// 	  t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
// 	  t_plan(2) = _ts(_bjx1-1) - 0.0045;
	  t_plan(0) = t_des - _dt;
	  t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
// 	  t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 - 0.01;
	  t_plan(2) = _ts(_bjx1-1) - 0.001;	  
	  
	   
	  if (abs(t_des - _ts(_bjx1-1)) <= ( + 0.0005))
	  {
	    _Rfootx(j_index) = _footxyz_real(0,_bjxx); 
	    _Rfooty(j_index) = _footxyz_real(1,_bjxx);
	    _Rfootz(j_index) = _footxyz_real(2,_bjxx); 	
	    
	    _Rfootx(j_index+1) = _footxyz_real(0,_bjxx); 
	    _Rfooty(j_index+1) = _footxyz_real(1,_bjxx);
	    _Rfootz(j_index+1) = _footxyz_real(2,_bjxx); 	  	    
	  }
	  else
	  {
	    Eigen::Matrix<double,7,7> AAA_inv = solve_AAA_inv_x(t_plan);
		    
	    Eigen::Matrix<double, 1, 7> t_a_plan;
	    t_a_plan.setZero();
	    t_a_plan(0) = pow(t_des, 6);   t_a_plan(1) = pow(t_des, 5);   t_a_plan(2) = pow(t_des, 4);  t_a_plan(3) = pow(t_des, 3);
	    t_a_plan(4) = pow(t_des, 2);   t_a_plan(5) = pow(t_des, 1);   t_a_plan(6) = 1;
	    

	    Eigen::Matrix<double, 1, 7> t_a_planv;
	    t_a_planv.setZero();
	    t_a_planv(0) = 6*pow(t_des, 5);   t_a_planv(1) = 5*pow(t_des, 4);   t_a_planv(2) = 4*pow(t_des, 3);  t_a_planv(3) = 3*pow(t_des, 2);
	    t_a_planv(4) = 2*pow(t_des, 1);   t_a_planv(5) = 1;                 t_a_planv(6) = 0;
	    
	    
	    Eigen::Matrix<double, 1, 7> t_a_plana;
	    t_a_plana.setZero();
	    t_a_plana(0) = 30*pow(t_des, 4);   t_a_plana(1) = 20*pow(t_des, 3);   t_a_plana(2) = 12*pow(t_des, 2);  t_a_plana(3) = 6*pow(t_des, 1);
	    t_a_plana(4) = 2;                  t_a_plana(5) = 0;                  t_a_plana(6) = 0;
	    
  // 	  cout<<"AAA="<<endl<<AAA<<endl;
  // 	  cout <<"AAA_inverse="<<endl<<AAA.inverse()<<endl;	  
  // 	  cout <<"t_des="<<endl<<t_des<<endl;
  // 	  cout <<"t_plan="<<endl<<t_plan<<endl;
	    
	    ////////////////////////////////////////////////////////////////////////////
	    Eigen::Matrix<double, 7, 1> Rfootx_plan;
	    Rfootx_plan.setZero();	
	    Rfootx_plan(0) = _Rfootvx(j_index-1);     Rfootx_plan(1) = _Rfootax(j_index-1); Rfootx_plan(2) = _Rfootx(j_index-1); Rfootx_plan(3) = _Lfootx(j_index);
	    Rfootx_plan(4) = _footxyz_real(0,_bjxx);  Rfootx_plan(5) = 0;                   Rfootx_plan(6) = 0;
	    
	    
	    Eigen::Matrix<double, 7, 1> Rfootx_co;
	    Rfootx_co.setZero();
	    Rfootx_co = AAA_inv * Rfootx_plan;
	    
	    _Rfootx(j_index) = t_a_plan * Rfootx_co;
	    _Rfootvx(j_index) = t_a_planv * Rfootx_co;
	    _Rfootax(j_index) = t_a_plana * Rfootx_co;
	    
	    /////////////////////////////////////////////////////////////////////////////
	    if ((j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1)+_dt)
	    {
	      _ry_left_right = (_footxyz_real(1,_bjxx) + _footxyz_real(1,_bjxx-2))/2;
	    }
	    
	    Eigen::Matrix<double, 7, 1> Rfooty_plan;
	    Rfooty_plan.setZero();	
	    Rfooty_plan(0) = _Rfootvy(j_index-1);     Rfooty_plan(1) = _Rfootay(j_index-1); Rfooty_plan(2) = _Rfooty(j_index-1); Rfooty_plan(3) = _ry_left_right;
	    Rfooty_plan(4) = _footxyz_real(1,_bjxx);  Rfooty_plan(5) = 0;                   Rfooty_plan(6) = 0;	    
	    
	    Eigen::Matrix<double, 7, 1> Rfooty_co;
	    Rfooty_co.setZero();
	    Rfooty_co = AAA_inv * Rfooty_plan;
	    
	    _Rfooty(j_index) = t_a_plan * Rfooty_co;
	    _Rfootvy(j_index) = t_a_planv * Rfooty_co;
	    _Rfootay(j_index) = t_a_plana * Rfooty_co;	
	    
	    
	    //////////////////////////////////////////////////////////
	    Eigen::Matrix<double, 7, 1> Rfootz_plan;
	    Rfootz_plan.setZero();	
	    Rfootz_plan(0) = _Rfootvz(j_index-1);     Rfootz_plan(1) = _Rfootaz(j_index-1); Rfootz_plan(2) = _Rfootz(j_index-1); Rfootz_plan(3) = _Lfootz(j_index)+_lift_height_ref(_bjx1-1);
	    Rfootz_plan(4) = _footxyz_real(2,_bjxx);  Rfootz_plan(5) = 0;                   Rfootz_plan(6) = 0.000;	
	    
	    Eigen::Matrix<double, 7, 1> Rfootz_co;
	    Rfootz_co.setZero();
	    Rfootz_co = AAA_inv * Rfootz_plan;
	    
	    _Rfootz(j_index) = t_a_plan * Rfootz_co;
	    _Rfootvz(j_index) = t_a_planv * Rfootz_co;
	    _Rfootaz(j_index) = t_a_plana * Rfootz_co;	
	  
  	  _Rfootx(j_index+1) = _Rfootx(j_index)+_dt * _Rfootvx(j_index);
  	  _Rfooty(j_index+1) = _Rfooty(j_index)+_dt * _Rfootvy(j_index);
  	  _Rfootz(j_index+1) = _Rfootz(j_index)+_dt * _Rfootvz(j_index);
	  }
    }
	}
    else                       //right support
    {
  //       no change on right support
	_Rfootx(j_index) = _Rfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Rfooty(j_index) = _Rfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Rfootz(j_index) = _Rfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	
	_Rfootx(j_index+1) = _Rfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Rfooty(j_index+1) = _Rfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	_Rfootz(j_index+1) = _Rfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);      
    
	/// left swing
	if ((j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double suppot
	{
  // 	cout << "dsp"<<endl;
	  _Lfootx(j_index) = _Lfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Lfooty(j_index) = _Lfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Lfootz(j_index) = _Lfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);

	  _Lfootx(j_index+1) = _Lfootx((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Lfooty(j_index+1) = _Lfooty((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  _Lfootz(j_index+1) = _Lfootz((int)boost::math::iround(_tx(_bjx1-1)/_dt) -1-1);
	  
	}
	else
	{
  // 	cout << "ssp"<<endl;
	  //initial state and final state and the middle state
	  double t_des = (j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt) +1)*_dt;
	  Eigen::Vector3d t_plan(0,0,0);
// 	  t_plan(0) = t_des - _dt;
// 	  t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
// 	  t_plan(2) = _ts(_bjx1-1) - 0.0045;
	  t_plan(0) = t_des - _dt;
	  t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
// 	  t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 - 0.01;
	  t_plan(2) = _ts(_bjx1-1) - 0.001;	  
	  
	  
	  if (abs(t_des - _ts(_bjx1-1)) <= ( + 0.0005))
	  {
	    
	    _Lfootx(j_index) = _footxyz_real(0,_bjxx); 
	    _Lfooty(j_index) = _footxyz_real(1,_bjxx);
	    _Lfootz(j_index) = _footxyz_real(2,_bjxx); 

	    _Lfootx(j_index+1) = _footxyz_real(0,_bjxx); 
	    _Lfooty(j_index+1) = _footxyz_real(1,_bjxx);
	    _Lfootz(j_index+1) = _footxyz_real(2,_bjxx); 
	    
	  }
	  else
	  {
	    Eigen::Matrix<double,7,7> AAA_inv = solve_AAA_inv_x(t_plan);
	    
	    Eigen::Matrix<double, 1, 7> t_a_plan;
	    t_a_plan.setZero();
	    t_a_plan(0) = pow(t_des, 6);   t_a_plan(1) = pow(t_des, 5);   t_a_plan(2) = pow(t_des, 4);  t_a_plan(3) = pow(t_des, 3);
	    t_a_plan(4) = pow(t_des, 2);   t_a_plan(5) = pow(t_des, 1);   t_a_plan(6) = 1;
	    

	    Eigen::Matrix<double, 1, 7> t_a_planv;
	    t_a_planv.setZero();
	    t_a_planv(0) = 6*pow(t_des, 5);   t_a_planv(1) = 5*pow(t_des, 4);   t_a_planv(2) = 4*pow(t_des, 3);  t_a_planv(3) = 3*pow(t_des, 2);
	    t_a_planv(4) = 2*pow(t_des, 1);   t_a_planv(5) = 1;                 t_a_planv(6) = 0;
	    
	    
	    Eigen::Matrix<double, 1, 7> t_a_plana;
	    t_a_plana.setZero();
	    t_a_plana(0) = 30*pow(t_des, 4);   t_a_plana(1) = 20*pow(t_des, 3);   t_a_plana(2) = 12*pow(t_des, 2);  t_a_plana(3) = 6*pow(t_des, 1);
	    t_a_plana(4) = 2;                  t_a_plana(5) = 0;                  t_a_plana(6) = 0;
	    
	    
	    ////////////////////////////////////////////////////////////////////////////
	    Eigen::Matrix<double, 7, 1> Lfootx_plan;
	    Lfootx_plan.setZero();	
	    Lfootx_plan(0) = _Lfootvx(j_index-1);     Lfootx_plan(1) = _Lfootax(j_index-1); Lfootx_plan(2) = _Lfootx(j_index-1); Lfootx_plan(3) = _Rfootx(j_index);
	    Lfootx_plan(4) = _footxyz_real(0,_bjxx);  Lfootx_plan(5) = 0;                   Lfootx_plan(6) = 0;	  
	    
	    
	    Eigen::Matrix<double, 7, 1> Lfootx_co;
	    Lfootx_co.setZero();
	    Lfootx_co = AAA_inv * Lfootx_plan;
	    
	    _Lfootx(j_index) = t_a_plan * Lfootx_co;
	    _Lfootvx(j_index) = t_a_planv * Lfootx_co;
	    _Lfootax(j_index) = t_a_plana * Lfootx_co;
	    
	    /////////////////////////////////////////////////////////////////////////////
	    if ((j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1)+_dt)
	    {
	      _ry_left_right = (_footxyz_real(1,_bjxx) + _footxyz_real(1,_bjxx-2))/2;
	    }
	    
	    Eigen::Matrix<double, 7, 1> Lfooty_plan;
	    Lfooty_plan.setZero();	
	    Lfooty_plan(0) = _Lfootvy(j_index-1);     Lfooty_plan(1) = _Lfootay(j_index-1); Lfooty_plan(2) = _Lfooty(j_index-1); Lfooty_plan(3) = _ry_left_right;
	    Lfooty_plan(4) = _footxyz_real(1,_bjxx);  Lfooty_plan(5) = 0;                   Lfooty_plan(6) = 0;		  
	    
	    
	    Eigen::Matrix<double, 7, 1> Lfooty_co;
	    Lfooty_co.setZero();
	    Lfooty_co = AAA_inv * Lfooty_plan;
	    
	    _Lfooty(j_index) = t_a_plan * Lfooty_co;
	    _Lfootvy(j_index) = t_a_planv * Lfooty_co;
	    _Lfootay(j_index) = t_a_plana * Lfooty_co;	
	    
	    
	    //////////////////////////////////////////////////////////
	    Eigen::Matrix<double, 7, 1> Lfootz_plan;
	    Lfootz_plan.setZero();		
	    Lfootz_plan(0) = _Lfootvz(j_index-1);     Lfootz_plan(1) = _Lfootaz(j_index-1); Lfootz_plan(2) = _Lfootz(j_index-1); Lfootz_plan(3) = _Rfootz(j_index)+_lift_height_ref(_bjx1-1);
	    Lfootz_plan(4) = _footxyz_real(2,_bjxx);  Lfootz_plan(5) = 0;                   Lfootz_plan(6) = 0.000;		  
	    
	    
	    Eigen::Matrix<double, 7, 1> Lfootz_co;
	    Lfootz_co.setZero();
	    Lfootz_co = AAA_inv * Lfootz_plan;
	    
	    _Lfootz(j_index) = t_a_plan * Lfootz_co;
	    _Lfootvz(j_index) = t_a_planv * Lfootz_co;
	    _Lfootaz(j_index) = t_a_plana * Lfootz_co;
	    
  	  
	    _Lfootx(j_index+1) = _Lfootx(j_index)+_dt * _Lfootvx(j_index);
	    _Lfooty(j_index+1) = _Lfooty(j_index)+_dt * _Lfootvy(j_index);
	    _Lfootz(j_index+1) = _Lfootz(j_index)+_dt * _Lfootvz(j_index);
	    
	    

// 	    if (j_index>=73)
// 	    {
// 	      cout<<"_Lfootx_plan:"<<Lfootx_plan<<endl;
// 	      cout<<"_Lfooty_plan:"<<Lfooty_plan<<endl;
// 	      cout<<"t_a_plan:"<<t_a_plan<<endl;
// 	      cout<<"_tx(_bjx1-1)"<<_tx(_bjx1-1)<<endl;
// 	      cout<<"(_bjx1-1)"<<_bjx1-1<<endl;
// 	      cout<<"t-des"<<(j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt) +1)*_dt<<endl;
// 	    cout<<"j_index:"<<j_index<<endl;	  
// 	    cout<<"t-des"<<t_des<<endl;
// 	    cout<<"_bjx1:"<<_bjx1<<endl;
// 	    cout<<"_Lfootx(j_index):"<<_Lfootx(j_index)<<endl;
// 	    cout<<"_Lfooty(j_index):"<<_Lfooty(j_index)<<endl;	  	      
// 	      
// 	    }
	  }

	}

    }

  }
  else
  {
    _Rfooty(j_index) = -_stepwidth(0);
    _Lfooty(j_index) = _stepwidth(0);
  }
}



Vector3d MPCClass::X_CoM_position_squat(int walktime, double dt_sample)
{
  Vector3d com_inte;
  com_inte.setZero();


  double t_des;
  t_des = walktime * dt_sample;

  Eigen::Vector3d t_plan;
  t_plan(0) = 0.00001;
  t_plan(1) = _height_squat_time/2+0.0001;
  t_plan(2) = _height_squat_time+0.0001;


  if (t_des<=_height_squat_time)
  {
    Eigen::Matrix<double,7,7> AAA_inv = solve_AAA_inv_x(t_plan);
	    
    Eigen::Matrix<double, 1, 7> t_a_plan;
    t_a_plan.setZero();
    t_a_plan(0) = pow(t_des, 6);   t_a_plan(1) = pow(t_des, 5);   t_a_plan(2) = pow(t_des, 4);  t_a_plan(3) = pow(t_des, 3);
    t_a_plan(4) = pow(t_des, 2);   t_a_plan(5) = pow(t_des, 1);   t_a_plan(6) = 1;

    ////////////////////////////////////////////////////////////////////////////
    Eigen::Matrix<double, 7, 1> com_squat_plan;
    com_squat_plan.setZero();	
    com_squat_plan(0) = 0;     com_squat_plan(1) = 0; 
    com_squat_plan(2) = _HCOM; 
    com_squat_plan(3) = _HCOM-_height_offset/2;
    com_squat_plan(4) = _HCOM-_height_offset;  
    com_squat_plan(5) = 0;                   
    com_squat_plan(6) = 0;


    Eigen::Matrix<double, 7, 1> com_co;
    com_co.setZero();
    com_co = AAA_inv * com_squat_plan;

    com_inte(2) = t_a_plan * com_co;
  }
  else
  {
    com_inte(2) = _HCOM-_height_offset;
  }
  
  return com_inte;
}

///////////////////////// ODE - sampling time maximal ========================================
int MPCClass::Get_maximal_number_reference()
{
  return (_nsum -_nh-1);
}

int MPCClass::Get_maximal_number(double dtx)
{
  
  return (_nsum -_nh-1)*(int)(floor(_dt/dtx));
}

////====================================================================================================================
/////////////////////////// using the lower-level control-loop  sampling time as the reference: every 5ms;  at the same time: just using the next one position + next one velocity

Vector3d MPCClass::XGetSolution_CoM_position(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
  //reference com position
	_comz(0) = _HCOM-_height_offset;
	_comz(1) = _HCOM-_height_offset;
	_comz(2) = _HCOM-_height_offset;
	_comz(3) = _HCOM-_height_offset;
	_comz(4) = _HCOM-_height_offset;
	
	Vector3d com_inte(0,0,0);	
	
	if (walktime>=2)
	{
	  int t_int= (int) floor(walktime / (_dt / dt_sample) );

	  ///// chage to be relative time
	  double t_cur = walktime * dt_sample ;
	  
	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = 0;
	  t_plan(1) = dt_sample;
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);
	  
	  solve_AAA_inv(t_plan);
	  	  
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(2*dt_sample, 3);   t_a_plan(1) = pow(2*dt_sample, 2);   
	  t_a_plan(2) = pow(2*dt_sample, 1);   t_a_plan(3) = pow(2*dt_sample, 0); 
	  
	  Eigen::Matrix<double, 1, 4> t_a_planv;
	  t_a_planv.setZero();
	  t_a_planv(0) = 3*pow(2*dt_sample, 2);   t_a_planv(1) = 2*pow(2*dt_sample, 1);   
	  t_a_planv(2) = 1;   t_a_planv(3) = 0; 	  
	  
/*	  Eigen::Matrix<double, 1, 4> t_a_plana;
	  t_a_plana.setZero();
	  t_a_plana(0) = 6*pow(2*dt_sample, 1);   t_a_plana(1) = 2;   
	  t_a_plana(2) = 0;   t_a_plana(3) = pow(2*dt_sample, 0); */	  
	  	  
	  
	  // COM&&foot trajectory interpolation	   	  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = body_in1(0); temp(1) = body_in2(0); temp(2) = _comx(t_int); temp(3) = _comvx(t_int);	  
	  com_inte(0) = t_a_plan * (_AAA_inv)*temp;
	  _comxyzx(0) = com_inte(0);
	  _comvxyzx(0) = t_a_planv * (_AAA_inv)*temp;
//	  _comaxyzx(0) = t_a_plana * (_AAA_inv)*temp;
	  
	  
	  temp(0) = body_in1(1); temp(1) = body_in2(1); temp(2) = _comy(t_int); temp(3) = _comvy(t_int);	  
	  com_inte(1) = t_a_plan * (_AAA_inv)*temp;
	  _comxyzx(1) = com_inte(1);
	  _comvxyzx(1) = t_a_planv * (_AAA_inv)*temp;
//	  _comaxyzx(1) = t_a_plana * (_AAA_inv)*temp;	  
	  
	  temp(0) = body_in1(2); temp(1) = body_in2(2); temp(2) = _comz(t_int); temp(3) = _comvz(t_int);	  
	  com_inte(2) = t_a_plan *(_AAA_inv)*temp;
	  _comxyzx(2) = com_inte(2);
	  _comvxyzx(2) = t_a_planv * (_AAA_inv)*temp;
//	  _comaxyzx(2) = t_a_plana * (_AAA_inv)*temp;	 
	  
	  /////be careful, the polynomial may cause overfitting
	  double t_des = t_cur-t_int*_dt;
	  if (t_des<=0){
	    t_des =0.00001;
	  }
	    
	  if (t_int>=1)
	  {
	    _comaxyzx(0) = (_comax(t_int)-_comax(t_int-1))/_dt*t_des+_comax(t_int-1);
	    _comaxyzx(1) = (_comay(t_int)-_comay(t_int-1))/_dt*t_des+_comay(t_int-1);
	    _comaxyzx(2) = (_comaz(t_int)-_comaz(t_int-1))/_dt*t_des+_comaz(t_int-1);
	  }
	  else
	  {
	    _comaxyzx(0) = (_comax(t_int)-0)/_dt*t_des+0;
	    _comaxyzx(1) = (_comay(t_int)-0)/_dt*t_des+0;
	    _comaxyzx(2) = (_comaz(t_int)-0)/_dt*t_des+0;	    
	  } 
	}
	else
	{
	  com_inte = body_in3;	
	  _comxyzx = com_inte;
	}

 	return com_inte;
	
}

Vector3d MPCClass::XGetSolution_body_inclination(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
  //reference com position
		
	Vector3d com_inte(0,0,0);		
	if (walktime>=2)
	{
	  int t_int = (int) floor(walktime / (_dt / dt_sample) );
	 
	  
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(2*dt_sample, 3);   t_a_plan(1) = pow(2*dt_sample, 2);   
	  t_a_plan(2) = pow(2*dt_sample, 1);   t_a_plan(3) = pow(2*dt_sample, 0); 
	  
	  Eigen::Matrix<double, 1, 4> t_a_planv;
	  t_a_planv.setZero();
	  t_a_planv(0) = 3*pow(2*dt_sample, 2);   t_a_planv(1) = 2*pow(2*dt_sample, 1);   
	  t_a_planv(2) = 1;   t_a_planv(3) = 0; 	  
	  
/*	  Eigen::Matrix<double, 1, 4> t_a_plana;
	  t_a_plana.setZero();
	  t_a_plana(0) = 6*pow(2*dt_sample, 1);   t_a_plana(1) = 2;   
	  t_a_plana(2) = 0;   t_a_plana(3) = pow(2*dt_sample, 0); */	  
	  
	  //body inclination interpolation	  	  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = body_in1(0); temp(1) = body_in2(0); temp(2) = _thetax(t_int); temp(3) = _thetavx(t_int);	  
	  com_inte(0) = t_a_plan * (_AAA_inv)*temp;
	  _thetaxyx(0) = com_inte(0);
	  _thetavxyx(0) = t_a_planv * (_AAA_inv)*temp;
//	  _thetaaxyx(0) = t_a_plana * (_AAA_inv)*temp;	  
	  
	  temp(0) = body_in1(1); temp(1) = body_in2(1); temp(2) = _thetay(t_int); temp(3) = _thetavy(t_int);	  
	  com_inte(1) = t_a_plan * (_AAA_inv)*temp;
	  _thetaxyx(1) = com_inte(1);
	  _thetavxyx(1) = t_a_planv * (_AAA_inv)*temp;
//	  _thetaaxyx(1) = t_a_plana * (_AAA_inv)*temp;		  
	  
// 	  temp(0) = body_in1(2); temp(1) = body_in2(2); temp(2) = _thetaz(t_int); temp(3) = _thetavz(t_int);	  
// 	  com_inte(2) = t_a_plan *(_AAA_inv)*temp;
// 	  _thetaxyx(2) = com_inte(2);
// 	  _thetavxyx(2) = t_a_planv * (_AAA_inv)*temp;
// 	  _thetaaxyx(2) = t_a_plana * (_AAA_inv)*temp;
	  
	  /////be careful, the polynomial may cause overfitting
	  double t_des = walktime * dt_sample-t_int*_dt;
	  if (t_des<=0){
	    t_des =0.0001;
	  }	  

	  if (t_int>=1)
	  {
	    _thetaxyx(0) = (_thetaax(t_int)-_thetaax(t_int-1))/_dt*t_des+_thetaax(t_int-1);
	    _thetaxyx(1) = (_thetaay(t_int)-_thetaay(t_int-1))/_dt*t_des+_thetaay(t_int-1);
	    _thetaxyx(2) = (_thetaaz(t_int)-_thetaaz(t_int-1))/_dt*t_des+_thetaaz(t_int-1);
	  }
	  else
	  {
	    _thetaxyx(0) = (_thetaax(t_int)-0)/_dt*t_des+0;
	    _thetaxyx(1) = (_thetaay(t_int)-0)/_dt*t_des+0;
	    _thetaxyx(2) = (_thetaaz(t_int)-0)/_dt*t_des+0;	    
	  } 	  
	}
	else
	{
	  com_inte = body_in3;
	  _thetaxyx = com_inte;
	}
 	return com_inte;
}

Vector3d MPCClass::XGetSolution_Foot_positionR(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
	
	Vector3d com_inte(0,0,0);	
	
	if (walktime>=2)
	{
	  int t_int = (int) floor(walktime / (_dt / dt_sample) );

	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(2*dt_sample, 3);   t_a_plan(1) = pow(2*dt_sample, 2);   
	  t_a_plan(2) = pow(2*dt_sample, 1);   t_a_plan(3) = pow(2*dt_sample, 0);  
  
	//  foot trajectory interpolation	  
	  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = body_in1(0); temp(1) = body_in2(0); temp(2) = _Rfootx(t_int); temp(3) = _Rfootvx(t_int);	  
	  com_inte(0) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(1); temp(1) = body_in2(1); temp(2) = _Rfooty(t_int); temp(3) = _Rfootvy(t_int);
	  com_inte(1) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(2); temp(1) = body_in2(2); temp(2) = _Rfootz(t_int); temp(3) = _Rfootvz(t_int);	  
	  com_inte(2) = t_a_plan *(_AAA_inv)*temp;	    
	}
	else
	{
          com_inte = body_in3;	    
	}
	
	_Rfootxyzx = com_inte;
	
 	return com_inte;
         	
}

Vector3d MPCClass::XGetSolution_Foot_positionL(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
	
	Vector3d com_inte(0,0,0);	
	
	if (walktime>=2)
	{
	  int t_int = (int) floor(walktime / (_dt / dt_sample) );
	  	  
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(2*dt_sample, 3);   t_a_plan(1) = pow(2*dt_sample, 2);   
	  t_a_plan(2) = pow(2*dt_sample, 1);   t_a_plan(3) = pow(2*dt_sample, 0); 
	  
	//  foot trajectory interpolation
	  
	  Eigen::Vector4d  temp;
	  temp.setZero();
	  temp(0) = body_in1(0); temp(1) = body_in2(0); temp(2) = _Lfootx(t_int); temp(3) = _Lfootvx(t_int);
	  com_inte(0) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(1); temp(1) = body_in2(1); temp(2) = _Lfooty(t_int); temp(3) = _Lfootvy(t_int);	  
	  com_inte(1) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(2); temp(1) = body_in2(2); temp(2) = _Lfootz(t_int); temp(3) = _Lfootvz(t_int);		  
	  com_inte(2) = t_a_plan * (_AAA_inv)*temp;
	  
	}
	else
	{
	  com_inte = body_in3;	  
	}

	_Lfootxyzx = com_inte;
	
 	return com_inte;
	
}

////solve the inverse matrix of 4*4 coefficient matrices
void MPCClass::solve_AAA_inv(const Eigen::Matrix<double, 4, 1>& t_plan)
{
// 	  Eigen::MatrixXd AAA1;	
// 
// 	  AAA1.setZero(4,4);	
// 	  AAA1(0,0) = pow(t_plan(0), 3); AAA1(0,1) = pow(t_plan(0), 2); AAA1(0,2) = pow(t_plan(0), 1); AAA1(0,3) = pow(t_plan(0), 0); 
// 	  AAA1(1,0) = pow(t_plan(1), 3); AAA1(1,1) = pow(t_plan(1), 2); AAA1(1,2) = pow(t_plan(1), 1); AAA1(1,3) = pow(t_plan(0), 0); 
// 	  AAA1(2,0) = pow(t_plan(2), 3); AAA1(2,1) = pow(t_plan(2), 2); AAA1(2,2) = pow(t_plan(2), 1); AAA1(2,3) = pow(t_plan(0), 0); 
// 	  AAA1(3,0) = 3*pow(t_plan(2), 2); AAA1(3,1) = 2*pow(t_plan(2), 1); AAA1(3,2) = pow(t_plan(2), 0); AAA1(3,3) = 0;  

  
  
  double abx1 = ((t_plan(0) - t_plan(1))*pow(t_plan(0) - t_plan(2), 2));
  double abx2 = ((t_plan(0) - t_plan(1))*pow(t_plan(1) - t_plan(2), 2));
  double abx3 =(pow(t_plan(0) - t_plan(2), 2)*pow(t_plan(1) - t_plan(2), 2));
  double abx4 = ((t_plan(0) - t_plan(2))*(t_plan(1) - t_plan(2)));
  

  _AAA_inv(0,0) = 1/ abx1;
  _AAA_inv(0,1) =  -1/ abx2;
  _AAA_inv(0,2) = (t_plan(0) + t_plan(1) - 2*t_plan(2))/ abx3;
  _AAA_inv(0,3) = 1/ abx4;
  
  _AAA_inv(1,0) = -(t_plan(1) + 2*t_plan(2))/ abx1;
  _AAA_inv(1,1) = (t_plan(0) + 2*t_plan(2))/ abx2;
  _AAA_inv(1,2) = -(pow(t_plan(0), 2) + t_plan(0)*t_plan(1) + pow(t_plan(1), 2) - 3*pow(t_plan(2), 2))/ abx3;
  _AAA_inv(1,3) = -(t_plan(0) + t_plan(1) + t_plan(2))/ abx4;
  
  _AAA_inv(2,0) = (t_plan(2)*(2*t_plan(1) + t_plan(2)))/ abx1;
  _AAA_inv(2,1) = -(t_plan(2)*(2*t_plan(0) + t_plan(2)))/ abx2;
  _AAA_inv(2,2) = (t_plan(2)*(2*pow(t_plan(0), 2) + 2*t_plan(0)*t_plan(1) - 3*t_plan(2)*t_plan(0) + 2*pow(t_plan(1), 2) - 3*t_plan(2)*t_plan(1)))/ abx3;
  _AAA_inv(2,3) = (t_plan(0)*t_plan(1) + t_plan(0)*t_plan(2) + t_plan(1)*t_plan(2))/ abx4;
  
  _AAA_inv(3,0) = -(t_plan(1)*pow(t_plan(2), 2))/ abx1;
  _AAA_inv(3,1) = (t_plan(0)*pow(t_plan(2), 2))/ abx2;
  _AAA_inv(3,2) = (t_plan(0)*t_plan(1)*(t_plan(0)*t_plan(1) - 2*t_plan(0)*t_plan(2) - 2*t_plan(1)*t_plan(2) + 3*pow(t_plan(2), 2)))/ abx3;
  _AAA_inv(3,3) = -(t_plan(0)*t_plan(1)*t_plan(2))/ abx4;   
}

////solve the inverse matrix of 7*7 coefficient matrices
Eigen::Matrix<double, 7, 7> MPCClass::solve_AAA_inv_x(Eigen::Vector3d t_plan)
{
  Eigen::Matrix<double,7,7> AAA;
  AAA.setZero();
  Eigen::Matrix<double, 1, 7> aaaa;
  aaaa.setZero();

  aaaa(0) = 6*pow(t_plan(0), 5);   aaaa(1) =  5*pow(t_plan(0), 4);  aaaa(2) =  4*pow(t_plan(0), 3);   aaaa(3) =  3*pow(t_plan(0), 2);
  aaaa(4) = 2*pow(t_plan(0), 1);   aaaa(5) = 1;                     aaaa(6) =  0;
  AAA.row(0) = aaaa;

  aaaa(0) = 30*pow(t_plan(0), 4);  aaaa(1) =  20*pow(t_plan(0), 3);  aaaa(2) =  12*pow(t_plan(0), 2);   aaaa(3) =  6*pow(t_plan(0), 1);
  aaaa(4) = 2;                     aaaa(5) = 0;                      aaaa(6) =  0;
  AAA.row(1) = aaaa;

  aaaa(0) = pow(t_plan(0), 6);     aaaa(1) =  pow(t_plan(0), 5);     aaaa(2) =  pow(t_plan(0), 4);   aaaa(3) =  pow(t_plan(0), 3);
  aaaa(4) = pow(t_plan(0), 2);     aaaa(5) = pow(t_plan(0), 1);      aaaa(6) =  1;	  
  AAA.row(2) = aaaa;

  aaaa(0) = pow(t_plan(1), 6);     aaaa(1) =  pow(t_plan(1), 5);     aaaa(2) =  pow(t_plan(1), 4);   aaaa(3) =  pow(t_plan(1), 3);
  aaaa(4) = pow(t_plan(1), 2);     aaaa(5) = pow(t_plan(1), 1);      aaaa(6) =  1;
  AAA.row(3) = aaaa;

  aaaa(0) = pow(t_plan(2), 6);     aaaa(1) =  pow(t_plan(2), 5);     aaaa(2) =  pow(t_plan(2), 4);   aaaa(3) =  pow(t_plan(2), 3);
  aaaa(4) = pow(t_plan(2), 2);     aaaa(5) = pow(t_plan(2), 1);      aaaa(6) =  1;
  AAA.row(4) = aaaa;	  

  aaaa(0) = 6*pow(t_plan(2), 5);   aaaa(1) =  5*pow(t_plan(2), 4);  aaaa(2) =  4*pow(t_plan(2), 3);   aaaa(3) =  3*pow(t_plan(2), 2);
  aaaa(4) = 2*pow(t_plan(2), 1);   aaaa(5) = 1;                     aaaa(6) =  0;
  AAA.row(5) = aaaa;

  aaaa(0) = 30*pow(t_plan(2), 4);  aaaa(1) =  20*pow(t_plan(2), 3);  aaaa(2) =  12*pow(t_plan(2), 2);   aaaa(3) =  6*pow(t_plan(2), 1);
  aaaa(4) = 2;                     aaaa(5) = 0;                      aaaa(6) =  0;
  AAA.row(6) = aaaa;  
  
  Eigen::Matrix<double,7,7> AAA_inv = AAA.inverse(); 
  
  return AAA_inv;
  
}



//////////////////////////////////////////=============================ZMP optimal distribution for lower level adimittance control=================================
////reference_force_torque_distribution========================================

void MPCClass::Zmp_distributor(int walktime, double dt_sample)
{
// //// judge if stop  
//   if(_stopwalking)  
//   {  
//     for (int i_t = _bjx1+1; i_t < _footstepsnumber; i_t++) {	  
//       _lift_height_ref(i_t) = 0;  
//     }	  
// 
//   }  
  
  int j_index = (int)floor(walktime / (_dt / dt_sample));
  
  zmp_interpolation(j_index,walktime,dt_sample);  

// reference_force_torque_distribution 
  if (_bjx1 >= 2)
  {
      if (_bjx1 % 2 == 0)           //odd:left support
      {  
	/// right swing
	if ((j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))
	{
	  int nTx_n = (int)boost::math::iround(_tx(_bjx1-1)/_dt);
	  int nTx_n_dsp = (int)boost::math::iround((_tx(_bjx1-1)+_td(_bjx1-1))/_dt);
	  Vector2d ZMP_init(0,0);
	  ZMP_init(0) = _zmpx_real(0,nTx_n-2);
	  ZMP_init(1) = _zmpy_real(0,nTx_n-2);
	  
	  Vector2d ZMP_end(0,0); 
	  ZMP_end(0) = _zmpx_real(0,nTx_n_dsp-1);
	  ZMP_end(1) = _zmpy_real(0,nTx_n_dsp-1);
	  
// 	  if (abs(ZMP_end(0)-ZMP_init(0))<=0.001) 
// 	  {
// 	      _Co_L(1,1) = abs((_ZMPxy_realx(1)-ZMP_init(1))/(ZMP_end(1)-ZMP_init(1)));
// 	      if (_Co_L(1,1) >1)
// 	      {
// 		_Co_L(1,1)=1;
// 	      }               
// 	      _Co_L(0,0) = _Co_L(1,1);                
// 	      _Co_L(2,2) = sqrt((pow(_Co_L(0,0),2)+pow(_Co_L(0,0),2))/2);  
// 	  }
// 	  else
// 	  {
	    _Co_L(0,0) = abs(((ZMP_end(1)-ZMP_init(1))*(_ZMPxy_realx(1)-ZMP_init(1))+(ZMP_end(0)-ZMP_init(0))*(_ZMPxy_realx(0)-ZMP_init(0)))/(pow(ZMP_end(1)-ZMP_init(1),2)+pow(ZMP_end(0)-ZMP_init(0),2)));
	    _Co_L(1,1) = _Co_L(0,0);
	      
	    if (_Co_L(0,0) >1)
	    {
	      _Co_L(0,0)=1;
	    }
	    if (_Co_L(1,1) >1)
	    {
	      _Co_L(1,1)=1;
	    }
	    _Co_L(2,2) = sqrt((pow(_Co_L(0,0),2)+pow(_Co_L(0,0),2))/2); 
/*	  }  */ 
	  
	  Matrix3d II;
	  II.setIdentity(3,3);
	  _Co_R = II-_Co_L;     

	  Force_torque_calculate(_comxyzx,_comaxyzx,_thetaaxyx,_Lfootxyzx,_Rfootxyzx);	  
	  
	}
	else
	{
	 //Swing leg with left support
	  _Co_L.setIdentity(3,3);
	  _Co_R.setZero();
	  Force_torque_calculate(_comxyzx,_comaxyzx,_thetaaxyx,_Lfootxyzx,_Rfootxyzx);  	  
	}	
      }      
      else                       //right support
      {	/// left swing
	if ((j_index +1 - (int)boost::math::iround(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double suppot
	{
	  int nTx_n = (int)boost::math::iround(_tx(_bjx1-1)/_dt);
	  int nTx_n_dsp = (int)boost::math::iround((_tx(_bjx1-1)+_td(_bjx1-1))/_dt);
	  Vector2d ZMP_init(0,0);
	  ZMP_init(0) = _zmpx_real(0,nTx_n-2);
	  ZMP_init(1) = _zmpy_real(0,nTx_n-2);
	  
	  Vector2d ZMP_end(0,0); 
	  ZMP_end(0) = _zmpx_real(0,nTx_n_dsp-1);
	  ZMP_end(1) = _zmpy_real(0,nTx_n_dsp-1);
	  
// 	  if (abs(ZMP_end(0)-ZMP_init(0))<=0.001) 
// 	  {
// 	      _Co_R(1,1) = abs((_ZMPxy_realx(1)-ZMP_init(1))/(ZMP_end(1)-ZMP_init(1)));
// 	      if (_Co_R(1,1) >1)
// 	      {
// 		_Co_R(1,1)=1;
// 	      }               
// 	      _Co_R(0,0) = _Co_R(1,1);                
// 	      _Co_R(2,2) = sqrt((pow(_Co_R(0,0),2)+pow(_Co_R(0,0),2))/2);  
// 	  }
// 	  else
// 	  {
	    _Co_R(0,0) = abs(((ZMP_end(1)-ZMP_init(1))*(_ZMPxy_realx(1)-ZMP_init(1))+(ZMP_end(0)-ZMP_init(0))*(_ZMPxy_realx(0)-ZMP_init(0)))/(pow(ZMP_end(1)-ZMP_init(1),2)+pow(ZMP_end(0)-ZMP_init(0),2)));
	    _Co_R(1,1) = _Co_R(0,0);
	      
	    if (_Co_R(0,0) >1)
	    {
	      _Co_R(0,0)=1;
	    }
	    if (_Co_R(1,1) >1)
	    {
	      _Co_R(1,1)=1;
	    }
	    _Co_R(2,2) = sqrt((pow(_Co_R(0,0),2)+pow(_Co_R(0,0),2))/2); 
/*	  }*/   
	  
	  Matrix3d II;
	  II.setIdentity(3,3);
	  _Co_L = II-_Co_R;     

	  Force_torque_calculate(_comxyzx,_comaxyzx,_thetaaxyx,_Lfootxyzx,_Rfootxyzx);	  	 	  
	}
	else
	{
	 //Swing leg with right support
	  _Co_R.setIdentity(3,3);
	  _Co_L.setZero();
	  Force_torque_calculate(_comxyzx,_comaxyzx,_thetaaxyx,_Lfootxyzx,_Rfootxyzx);  	  
	}
      }
  }
  else
  {
    if (_bjx1==0)   ///stand still:
    {
//       _Co_R = 0.5*Eigen::Vector3d::Zero().setIdentity(3,3);
//       _Co_L.setZero();
      _Co_R(0,0) = _Co_R(1,1) = _Co_R(2,2) = _Co_L(0,0) = _Co_L(1,1) = _Co_L(2,2) = 0.5;
      Force_torque_calculate(_comxyzx,_comaxyzx,_thetaaxyx,_Lfootxyzx,_Rfootxyzx);       
    }
    else   ///right support: double support in the whole walking pattern: assuming that the initial and final ZMP is located at the foot locations 
    {
		/// left swing
	Vector2d ZMP_init;
	ZMP_init(0) = _footxyz_real(0,_bjx1-1);
	ZMP_init(1) = _footxyz_real(1,_bjx1-1);
	
	Vector2d ZMP_end; 
	ZMP_end(0) = _footxyz_real(0,_bjx1);
	ZMP_end(1) = _footxyz_real(1,_bjx1);
	
// 	if (abs(ZMP_end(0)-ZMP_init(0))<=0.001) 
// 	{
// 	    _Co_R(1,1) = abs((_ZMPxy_realx(1)-ZMP_init(1))/(ZMP_end(1)-ZMP_init(1)));
// 	    if (_Co_R(1,1) >1)
// 	    {
// 	      _Co_R(1,1)=1;
// 	    }               
// 	    _Co_R(0,0) = _Co_R(1,1);                
// 	    _Co_R(2,2) = sqrt((pow(_Co_R(0,0),2)+pow(_Co_R(0,0),2))/2);  
// 	}
// 	else
// 	{
	  _Co_R(0,0) = abs(((ZMP_end(1)-ZMP_init(1))*(_ZMPxy_realx(1)-ZMP_init(1))+(ZMP_end(0)-ZMP_init(0))*(_ZMPxy_realx(0)-ZMP_init(0)))/(pow(ZMP_end(1)-ZMP_init(1),2)+pow(ZMP_end(0)-ZMP_init(0),2)));
	  _Co_R(1,1) = _Co_R(0,0);
	    
	  if (_Co_R(0,0) >1)
	  {
	    _Co_R(0,0)=1;
	  }
	  if (_Co_R(1,1) >1)
	  {
	    _Co_R(1,1)=1;
	  }
	  _Co_R(2,2) = sqrt((pow(_Co_R(0,0),2)+pow(_Co_R(0,0),2))/2); 
/*	}*/   
	
	Matrix3d II;
	II.setIdentity(3,3);
	_Co_L = II-_Co_R;     

	Force_torque_calculate(_comxyzx,_comaxyzx,_thetaaxyx,_Lfootxyzx,_Rfootxyzx);	  	 	       
    }
  }   
}


void MPCClass::zmp_interpolation(int t_int,int walktime, double dt_sample)
{
  //// calculate by the nonlinear model:
//   if (t_int>=1)
//   {
//     _ZMPxy_realx(0) = _comxyzx(0) - (_comxyzx(2) - _Zsc(t_int-1))/(_comaxyzx(2)+_ggg(0))*_comaxyzx(0) - _j_ini * _thetaaxyx(1)/(_mass * (_ggg(0) + _comaxyzx(2)));
//     _ZMPxy_realx(1) = _comxyzx(1) - (_comxyzx(2) - _Zsc(t_int-1))/(_comaxyzx(2)+_ggg(0))*_comaxyzx(1) + _j_ini * _thetaaxyx(0)/(_mass * (_ggg(0) + _comaxyzx(2)));     
//   }
//   else
//   {
//     _ZMPxy_realx(0) = _comxyzx(0) - (_comxyzx(2))/(_comaxyzx(2)+_ggg(0))*_comaxyzx(0) - _j_ini * _thetaaxyx(1)/(_mass * (_ggg(0) + _comaxyzx(2)));
//     _ZMPxy_realx(1) = _comxyzx(1) - (_comxyzx(2))/(_comaxyzx(2)+_ggg(0))*_comaxyzx(1) + _j_ini * _thetaaxyx(0)/(_mass * (_ggg(0) + _comaxyzx(2)));     
//   }
  
  //// linear interpolation of ZMP reference:  
  double t_des = walktime * dt_sample-t_int*_dt;
  if (t_des<=0){
    t_des =0.0001;
  }	  
//   cout <<"//////////////////////////////////////////"<<endl;
//   cout <<"t_des:"<<t_des<<endl;
//   cout <<"//////////////////////////////////////////"<<endl;
  if (t_int>=1)
  {
    _ZMPxy_realx(0) = (_zmpx_real(0,t_int)-_zmpx_real(0,t_int-1))/_dt*t_des+_zmpx_real(0,t_int-1);
    _ZMPxy_realx(1) = (_zmpy_real(0,t_int)-_zmpy_real(0,t_int-1))/_dt*t_des+_zmpy_real(0,t_int-1);

  }
  else
  {
    _ZMPxy_realx(0) = _zmpx_real(0,t_int)/_dt*t_des+0;
    _ZMPxy_realx(1) = _zmpy_real(0,t_int)/_dt*t_des+0;	    
  }  
  
  
}


void MPCClass::Force_torque_calculate(Vector3d comxyzx1,Vector3d comaxyzx1,Vector3d thetaaxyx1,Vector3d Lfootxyz1,Vector3d Rfootxyz1)
{
  Vector3d gra;
  gra << 0,0, -_ggg;
  
  Vector3d F_total = _mass * (comaxyzx1 - gra);
  
  Vector3d the3a;
  the3a << thetaaxyx1(0),thetaaxyx1(1),0;
  Vector3d L_total = _j_ini * the3a;
  
  _F_R = _Co_R * F_total;
  _F_L = _Co_L * F_total;
  
  Vector3d R_det_foot_com = Rfootxyz1 -  comxyzx1;
  
  Vector3d L_det_foot_com = Lfootxyz1 -  comxyzx1;
  
  Vector3d M_total = L_total - _F_R.cross(R_det_foot_com) - _F_L.cross(L_det_foot_com);
  
  _M_R = _Co_R*M_total;  
  _M_L = _Co_L*M_total;
  
}

