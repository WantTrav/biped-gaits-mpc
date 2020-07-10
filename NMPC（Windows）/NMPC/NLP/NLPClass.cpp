/*****************************************************************************
NLPClass.cpp
*****************************************************************************/
#include "NLP/NLPClass.h"
#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>

// #include "KMP/kmp.h"
// #include <armadillo>



using namespace Eigen;
using namespace std;
// using namespace arma;


NLPClass::NLPClass()                    ///declaration function
	: QPBaseClass()
	, _robot_name("")
	, _robot_mass(69)
	, _lift_height(0.06)
	, _GGG(9.8)
	,_HCOM(0.89)
	,_half_hip_width(0.103)
	,_foot_width(0.1)
	, _method_flag_nlp(0)
	,_n_end_walking(1)
	,_j_period(0)
	,_F_R(0,0,0)
	,_F_L(0,0,0)
	,_M_R(0,0,0)
	,_M_L(0,0,0)		
{  
//     cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx MPC READY"<<endl;
  
}


void NLPClass::FootStepInputs( double stepwidth, double steplength, double stepheight)
{	
	_steplength.setConstant(steplength);
	_steplength(0) = 0;
	_steplength(1) = steplength/2;
        _steplength(_footstepsnumber-1) = 0;		
	_steplength(_footstepsnumber-2) = steplength/2;
	
	_steplength(5) = 0.2;
	
	_steplength(8) = 0.2;
	
	
	_stepwidth.setConstant(stepwidth);
	_stepwidth(0) = _stepwidth(0)/2;	
	
	_stepheight.setConstant(stepheight);
	_stepheight(_footstepsnumber)=0;
	_stepheight(_footstepsnumber-1)=0;

	_stepheight(4) = 0.1;
	_stepheight(5) = 0.1;	
	_stepheight(6) = 0.1;
	
	_stepheight(8) = -0.1;
	_stepheight(9) = -0.1;	
	_stepheight(10) = -0.1;	
	
	
	_lift_height_ref.setConstant(_lift_height);
// 	_lift_height_ref(0) = 0.00;
// 	_lift_height_ref(1) = 0.04;
        _lift_height_ref(_footstepsnumber-1) = 0;	
        _lift_height_ref(_footstepsnumber-2) = 0.04; 	



}


void NLPClass::Initialize()
{
//      ///initialize	
// 	_query_kmp = zeros<vec>(1);
// 	_mean_kmp = zeros<vec>(6);
// 	
//   	_inDim_kmp  = 1; 	      		    //input dimension
// 	_outDim_kmp = 3+3; 	      		    //output dimension
// 	_pvFlag_kmp = 1;			    // output: pos (and vel)
// 	///////////// adjust KMP parameters
// 	_lamda_kmp  = 5, _kh_kmp = 0.75;	    	    //set kmp parameters 
// 
// 
//  	static char fileName1[]="/home/jiatao/Dropbox/cat_multi_cogimon/src/multi_strategy_walking/src/walking/KMP/referdata_swing.txt";  
// //        static char fileName1[]="/home/embedded/jiatao_catkin_ws/src/multi_strategy_walking/src/walking/KMP/referdata_swing.txt";  
// 
// 	_data_kmp.load( fileName1 );         	    // load original data
// 	///modify the data sample:
// //        cout<<"load mat finish"<<endl;
//         
//  	 int data_kmp_m = _data_kmp.n_rows;
// 	 
// 	 for (int xj=0; xj<data_kmp_m; xj++)
// 	 {
// 	   _data_kmp(xj,2) = _data_kmp(xj,2) -(- 0.0726)+(-_half_hip_width);
// 	 }
// 	
// 	
// 	kmp_leg_L.kmp_initialize(_data_kmp, _inDim_kmp, _outDim_kmp, _pvFlag_kmp,_lamda_kmp, _kh_kmp); // initialize kmp
// 	kmp_leg_R.kmp_initialize(_data_kmp, _inDim_kmp, _outDim_kmp, _pvFlag_kmp,_lamda_kmp, _kh_kmp); // initialize kmp
// 	
// 	
	
	
	

//        cout<<"kmp initial finish"<<endl;
        
        // ==step parameters initialize==: given by the inputs
        ////NP initialization for step timing optimization
       //// reference steplength and stepwidth for timing optimization  
 	// ==step loctions setup==
        _Lxx_ref = _steplength;
        _Lyy_ref = _stepwidth;	   
        // local coordinate
	_Lyy_ref(0) = 0;
	for (int j =0; j<_footstepsnumber;j++)
	{
	  _Lyy_ref(j) = (int)pow(-1,j)*_stepwidth(j);
	}
  
// 	reference footstep locations setup
        _footx_ref.setZero();
	_footy_ref.setZero();
	_footz_ref.setZero();		
  	for (int i = 1; i < _footstepsnumber; i++) {
 	  _footx_ref(i) = _footx_ref(i-1) + _steplength(i-1);
	  _footy_ref(i) = _footy_ref(i-1) + (int)pow(-1,i-1)*_stepwidth(i-1);   
	  _footz_ref(i) = _footz_ref(i-1) + _stepheight(i-1);
	}
	cout <<"_footz_ref"<<_footz_ref<<endl;
	_footx_offline = _footx_ref;
	_footy_offline = _footy_ref;
	_footz_offline = _footz_ref;
	
	/// support foot location feedbacke for step time optimization
	_footx_real_feed = _footx_ref;  
	_footy_real_feed = _footy_ref;	
	
	
 //	cout<<"_footy_ref:"<<_footy_ref<<endl;
	
	// == step cycle setup
	_ts.setConstant(_tstep);
	_ts(5) = 0.9;
	_ts(10) = 1;	
	
/*	_ts(10) = 1.2;*/	
	
	_td = 0.2*_ts;
	_tx.setZero();
  	for (int i = 1; i < _footstepsnumber; i++) {
 	  _tx(i) = _tx(i-1) + _ts(i-1);
	  _tx(i) = round(_tx(i)/_dt)*_dt -0.000001;	  
	}	

	//whole sampling time sequnece       
	_t.setLinSpaced(round(_tx(_footstepsnumber-1)/_dt),_dt,_tx(_footstepsnumber-1));

	//parameters
        _hcom = _HCOM-_height_offset;	
//         _hcom = 0.4668;
	_g = _GGG;	
	_Wn = sqrt(_g/_hcom);
	
	_Wndt = _Wn*_dt;	
	
        // COM state
	_comx.setZero(); _comvx.setZero(); _comax.setZero();
	_comy.setZero(); _comvy.setZero(); _comay.setZero();
	_comz.setConstant(_hcom); 
	_comvz.setZero(); 
	_comaz.setZero();
        ///actual steplengh,stepwidth and walking period
	_Lxx_ref_real.setZero(); 
	_Lyy_ref_real.setZero(); 
	_Ts_ref_real.setZero();
	

	cout<<"finish step parameters initialization:"<<endl;
         
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %% parameters for first nlp     -step timing adjustment and next one step location
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	_px.setZero(); _py.setZero(); _pz.setZero();
// 	_zmpvx.setZero(); _zmpvy.setZero(); 
	_COMx_is.setZero(); _COMx_es.setZero(); _COMvx_is.setZero(); 
	_COMy_is.setZero(); _COMy_es.setZero(); _COMvy_is.setZero();
	_comx_feed.setZero(); _comvx_feed.setZero(); _comax_feed.setZero();
	_comy_feed.setZero(); _comvy_feed.setZero(); _comay_feed.setZero();
	
	_Vari_ini.setZero(); //Lxx,Lyy,Tr1,Tr2,Lxx1,Lyy1,Tr11,Tr21;
	_vari_ini.setZero();


	
	  
	  _rad = 0.2; 	 
	  
	//step timing constraints:
	_t_min = 0.4; _t_max = 2;
	
	// swing foot velocity constraints	
	_footx_vmax = 3;
	_footx_vmin = -9.875;
	_footy_vmax = 2;
	_footy_vmin = -2;		
	
	// CoM acceleration velocity:
	_comax_max = 5;
	_comax_min = -5;
	_comay_max = 5;
	_comay_min = -5;	  
	//weight coefficient
// 	_aax = 500000000;        _aay  =700000000;
// 	_aaxv= 800000;           _aayv =1000000;         
// 	_bbx = 500000000;        _bby  =500000000;      
// 	_rr1 = 100000000;        _rr2  =100000000;        
// 	_aax1 =300;              _aay1 =300;            
// 	_bbx1 =1000000;          _bby1 =1000000;
// 	_rr11 =1000000;          _rr21 =1000000; 
	
        // for slow walking
	_aax = 1000000;          _aay  =1000000;
	_aaxv= 100000;           _aayv =50000;         
	_bbx = 500000000;        _bby  =1500000000;      
	_rr1 = 100000000;        _rr2  =100000000;        
	_aax1 =300;              _aay1 =100;            
	_bbx1 =1000000;          _bby1 =1000000;
	_rr11 =1000000;          _rr21 =1000000; 


	
	
/*	  /// zmp-constraints	
	  _zmpx_ub=(RobotParaClass::FOOT_LENGTH()/2+RobotParaClass::HIP_TO_ANKLE_X_OFFSET())*_ZMP_ratio;  
	  _zmpx_lb=(-(RobotParaClass::FOOT_LENGTH()/2-RobotParaClass::HIP_TO_ANKLE_X_OFFSET())*_ZMP_ratio);
	  _zmpy_ub=(RobotParaClass::FOOT_WIDTH()/2*_ZMP_ratio); 
	  _zmpy_lb=(-RobotParaClass::FOOT_WIDTH()/2*_ZMP_ratio);*/

        //  foot location constraints 
	_footx_max=0.35;
	_footx_min=-0.15;
// 	_footy_max=2*RobotParaClass::HALF_HIP_WIDTH() + 0.2; 
// 	_footy_min=RobotParaClass::HALF_HIP_WIDTH() - 0.03; 
	  	
	
	
	_footy_max=0.4; 
	_footy_min=0.11; 	
	
//	_mass = _robot_mass; 
	_j_ini = _robot_mass* pow(_rad,2);	
	
	///external force
	_FX =160;  _FY =120;
	_t_last = 0.5;
	_det_xa = _FX/_robot_mass;  _det_ya = _FY/_robot_mass; 
	_det_xv = _det_xa*_t_last; _det_yv = _det_ya*_t_last;
	_det_xp = pow(_det_xv,2)/(2*_det_xa); _det_yp = pow(_det_yv,2)/(2*_det_ya);
	
	
	
	_n_loop_omit = 2*round(_tstep/_dt);
	
	xyz0 = -1; //flag for find function 
	xyz1 = 0;  
	xyz2 = 1;
	_j_period = 0; // the number for step cycle indefind
		
	
	_periond_i = 0;  /// period number falls into which cycle 
	_ki = 0;
	_k_yu = 0;
	_Tk = 0; 
	
	_Lxx_refx=0;_Lyy_refy=0;
	_Lxx_refx1=0;_Lyy_refy1=0;
	
	_tr1_ref = 0;        _tr2_ref = 0;
	_tr1_ref1 = 0;  _tr2_ref1 = 0;  
	
	/// remaining time boundaries
	_tr1_min=0;   _tr2_min=0;  _tr1_max=0;  _tr2_max=0;
	_tr11_min=0;  _tr21_min=0; _tr11_max=0; _tr21_max=0;	
	
	/// selection matrix for variables
        _SS1.setZero(); _SS1(0) =1;
	_SS2.setZero(); _SS2(1) =1;
	_SS3.setZero(); _SS3(2) =1;
	_SS4.setZero(); _SS4(3) =1;
	_SS5.setZero(); _SS5(4) =1;
	_SS6.setZero(); _SS6(5) =1;
	_SS7.setZero(); _SS7(6) =1;
	_SS8.setZero(); _SS8(7) =1;
	
	
	
	_comvx_endref.setZero();
	_comvy_endref.setZero();
	
	_AxO.setZero();_BxO.setZero();_Cx.setZero(); 
	_Axv.setZero();_Bxv.setZero();_Cxv.setZero();
	_AyO.setZero();_ByO.setZero();_Cy.setZero(); 
	_Ayv.setZero();_Byv.setZero();_Cyv.setZero();	

	_SQ_goal0.setZero();	
	_SQ_goal.setZero();
	_SQ_goal1.setZero();
	_SQ_goal20.setZero();
	_SQ_goal2.setZero();	
	_SQ_goal3.setZero();	
	_Sq_goal.setZero();
	_Sq_goal1.setZero();	
	_Sq_goal2.setZero();
	_Sq_goal3.setZero();	
	_Ax.setZero(); 	_Ay.setZero();
	_Bx.setZero();  _By.setZero();
	_ixi.setZero(); _iyi.setZero();
	
	/// remaining time constraints: inequality constraints
	_trx1_up.setZero();  _trx1_lp.setZero();
	_trx2_up.setZero();  _trx2_lp.setZero();	
	_trx3_up.setZero();  _trx3_lp.setZero();	
	_trx4_up.setZero();  _trx4_lp.setZero();	
	_det_trx1_up.setZero();  _det_trx1_lp.setZero();
	_det_trx2_up.setZero();  _det_trx2_lp.setZero();	
	_det_trx3_up.setZero();  _det_trx3_lp.setZero();	
	_det_trx4_up.setZero();  _det_trx4_lp.setZero();	
	
	_trx.setZero();  _det_trx.setZero();
	
	//// tr1&tr2:equation constraints
	_trx12.setZero();     _trx121.setZero();	
	_det_trx12.setZero(); _det_trx121.setZero();		
	_trxx.setZero();      _det_trxx.setZero();

	 _Lxx_ref1_e.setZero();  _Lyy_ref1_e.setZero();  _tr1_ref1_e.setZero();  _tr2_ref1_e.setZero();  
	 _Lxx_ref_e.setZero();   _Lyy_ref_e.setZero();   _tr1_ref_e.setZero();   _tr2_ref_e.setZero();	
	 _Lxx_ref1_eq.setZero(); _Lyy_ref1_eq.setZero(); _tr1_ref1_eq.setZero(); _tr2_ref1_eq.setZero();
	 _Lxx_ref_eq.setZero(); _Lyy_ref_eq.setZero();   _tr1_ref_eq.setZero();  _tr2_ref_eq.setZero();



	
	

	
	_h_lx_up.setZero();  _h_lx_lp.setZero(); _h_ly_up.setZero(); _h_ly_lp.setZero();
	_h_lx_up1.setZero(); _h_lx_lp1.setZero(); _h_ly_up1.setZero(); _h_ly_lp1.setZero();
	_det_h_lx_up.setZero();_det_h_lx_lp.setZero();_det_h_ly_up.setZero();_det_h_ly_lp.setZero();
	_det_h_lx_up1.setZero();_det_h_lx_lp1.setZero();_det_h_ly_up1.setZero();_det_h_ly_lp1.setZero();
	_h_lx_upx.setZero(); _det_h_lx_upx.setZero();
	
	// swing foot velocity constraints
	_h_lvx_up.setZero();  _h_lvx_lp.setZero(); _h_lvy_up.setZero(); _h_lvy_lp.setZero();
	_h_lvx_up1.setZero(); _h_lvx_lp1.setZero(); _h_lvy_up1.setZero(); _h_lvy_lp1.setZero();
	_det_h_lvx_up.setZero();_det_h_lvx_lp.setZero();_det_h_lvy_up.setZero();_det_h_lvy_lp.setZero();
	_det_h_lvx_up1.setZero();_det_h_lvx_lp1.setZero();_det_h_lvy_up1.setZero();_det_h_lvy_lp1.setZero();	
	_h_lvx_upx.setZero(); _det_h_lvx_upx.setZero();	
	
	// CoM acceleration boundary
	_AA = 0; _CCx=0; _BBx=0; _CCy=0; _BBy=0;_AA1x=0;_AA2x=0;_AA3x=0;_AA1y=0;_AA2y=0;_AA3y=0;
	_CoM_lax_up.setZero();  _CoM_lax_lp.setZero();  _CoM_lay_up.setZero();  _CoM_lay_lp.setZero();
	_det_CoM_lax_up.setZero();  _det_CoM_lax_lp.setZero();  _det_CoM_lay_up.setZero();  _det_CoM_lay_lp.setZero();
	_CoM_lax_upx.setZero();
	_det_CoM_lax_upx.setZero();		
	
	
	
	/// CoM velocity_inremental boundary
	_VAA=0; _VCCx=0; _VBBx=0; _VCCy=0; _VBBy=0;_VAA1x=0;_VAA2x=0;_VAA3x=0;_VAA1y=0;_VAA2y=0;_VAA3y=0;
	_CoM_lvx_up.setZero();  _CoM_lvx_lp.setZero();  _CoM_lvy_up.setZero();  _CoM_lvy_lp.setZero();
	_det_CoM_lvx_up.setZero();  _det_CoM_lvx_lp.setZero();  _det_CoM_lvy_up.setZero();  _det_CoM_lvy_lp.setZero();
	_CoM_lvx_upx.setZero();
	_det_CoM_lvx_upx.setZero();	
	
	
	
	/// CoM initial velocity_ boundary
	_VAA1x1=0;_VAA2x1=0;_VAA3x1=0;_VAA1y1=0;_VAA2y1=0;_VAA3y1=0;
	_CoM_lvx_up1.setZero();  _CoM_lvx_lp1.setZero();  _CoM_lvy_up1.setZero();  _CoM_lvy_lp1.setZero();
	_det_CoM_lvx_up1.setZero();  _det_CoM_lvx_lp1.setZero();  _det_CoM_lvy_up1.setZero();  _det_CoM_lvy_lp1.setZero();
	_CoM_lvx_upx1.setZero();
	_det_CoM_lvx_upx1.setZero();	
	
	cout << "finish!!!!!!!!!!! initial for NLP parameters"<<endl;
	///////////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	///%%%%%%%%%%%%% foot trajectory geneartion
	//////////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        _t_f.setZero();	
	_bjx1 = 0;
	_bjx2 = 0;
        _mx = 0;	
	
	_bjxx = 0;
	
	
	_footxyz_real.setZero();
	
	
	_Lfootx.setZero(); _Lfooty.setConstant(_stepwidth(0));_Lfootz.setZero(); _Lfootvx.setZero(); _Lfootvy.setZero();_Lfootvz.setZero(); 
	_Lfootax.setZero(); _Lfootay.setZero();_Lfootaz.setZero();
	_Rfootx.setZero(); _Rfooty.setConstant(-_stepwidth(0));_Rfootz.setZero(); _Rfootvx.setZero(); _Rfootvy.setZero();_Rfootvz.setZero(); 
	_Rfootax.setZero(); _Rfootay.setZero();_Rfootaz.setZero();	
	_ry_left_right = 0;	
	
	
	
	
	
	/////// swing leg trajectory generation using kmp_initialize
	_Lfootx_kmp.setZero();  _Lfooty_kmp.setConstant(_stepwidth(0)); _Lfootz_kmp.setZero(); 
	_Lfootvx_kmp.setZero(); _Lfootvy_kmp.setZero();                 _Lfootvz_kmp.setZero(); 

	_Rfootx_kmp.setZero();  _Rfooty_kmp.setConstant(-_stepwidth(0));_Rfootz_kmp.setZero(); 
	_Rfootvx_kmp.setZero(); _Rfootvy_kmp.setZero();                 _Rfootvz_kmp.setZero(); 
	
	_Lfootx_kmp_old.setZero();  _Lfooty_kmp_old.setConstant(_stepwidth(0)); _Lfootz_kmp_old.setZero(); 
	_Lfootvx_kmp_old.setZero(); _Lfootvy_kmp_old.setZero();                 _Lfootvz_kmp_old.setZero(); 

	_Rfootx_kmp_old.setZero();  _Rfooty_kmp_old.setConstant(-_stepwidth(0));_Rfootz_kmp_old.setZero(); 
	_Rfootvx_kmp_old.setZero(); _Rfootvy_kmp_old.setZero();                 _Rfootvz_kmp_old.setZero(); 


	///// generated CoM final position
	if (_method_flag_nlp==2)
	{
	  ////number of inequality constraints: 8+8+8+4+4+4
	  int nVars = 8;
	  int nEqCon = 2;
	  int nIneqCon = 24 + 12;
	  resizeQP(nVars, nEqCon, nIneqCon);  
	  
	}
	else
	{
	  if (_method_flag_nlp==1)
	  {	    
	    int nVars = 8;
	    int nEqCon = 6;
	    int nIneqCon = 24 + 12;
	    resizeQP(nVars, nEqCon, nIneqCon);  	    
	    
	  }
	  else
	  {
	    int nVars = 8;
	    int nEqCon = 10;
	    int nIneqCon = 24 + 12;
	    resizeQP(nVars, nEqCon, nIneqCon);  	    
	    
	  }
	  
	}


	

	cout << "finish!!!!!!!!!!! initial for nlp_KMP parameters"<<endl;



	///////////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//////////////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %% pamameters for second MPC- ZMP movement, height variance and body inclination
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %% robot parameters		
	_zmpx_real.setZero(_nsum); _zmpy_real.setZero(_nsum);	
	
	
	_thetax.setZero(_nsum); _thetavx.setZero(_nsum); _thetaax.setZero(_nsum);
	_thetay.setZero(_nsum); _thetavy.setZero(_nsum); _thetaay.setZero(_nsum);
	_thetaz.setZero(_nsum); _thetavz.setZero(_nsum); _thetaaz.setZero(_nsum);
	
//	_torquex_real.setZero(_nsum); _torquey_real.setZero(_nsum);
	
//	_xk.setZero(3,_nsum); _yk.setZero(3,_nsum); _zk.setZero(3,_nsum);
//	_thetaxk.setZero(3,_nsum); _thetayk.setZero(3,_nsum);
//	_x_vacc_k.setZero(_nsum); _y_vacc_k.setZero(_nsum); _z_vacc_k.setZero(_nsum); 
//	_thetax_vacc_k.setZero(_nsum); _thetay_vacc_k.setZero(_nsum); 
	
        // ==initial parameters for MPC==
	_ggg.setConstant(1, _GGG);

//	_Hcom1.setConstant(_nh,1,_hcom);
        	
/*	_a.setZero(3,3);
	_a << 1, _dt, pow(_dt,2)/2,    
	      0,   1,            _dt,
	      0,   0,              1;
	_b.setZero(3,1);
	_b << pow(_dt,3)/6,    
	      pow(_dt,2)/2,
	               _dt;	
	_c.setZero(1,3);
	_c << 1,0,(-1)*_hcom /_g;	
	
	_cp.setZero(1,3);
	_cp(0,0) = 1;
	_cv.setZero(1,3);
	_cv(0,1) = 1;
	_ca.setZero(1,3);
	_ca(0,2) = 1;*/	
		

	//vertical height constraints
//	_z_max.setConstant(_nsum,0.1);
//	_z_min.setConstant(_nsum,-0.1);	
	
        //footz refer: height of step
// 	_Zsc.setZero(_nsum,1);			
//   	for (int i = 0; i < _nsum-1; i++) {
//           Indexfind(_t(i),xyz1);	  
// 	  _Zsc(i) = _footz_ref(_j_period);   
// 	  _j_period = 0;   
// 	}		

//         _yk.topRows(1).setConstant(_footy_ref(0)); 
// 	_zk.topRows(1).setConstant(_hcom);
	

	//predictive model
// 	_pps.setZero(_nh,3); _ppu.setZero(_nh,_nh);
// 	_pvs.setZero(_nh,3); _pvu.setZero(_nh,_nh);
// 	_pas.setZero(_nh,3); _pau.setZero(_nh,_nh);
// 
//         _pps = Matrix_ps(_a,_nh,_cp);
// 	_pvs = Matrix_ps(_a,_nh,_cv);
// 	_pas = Matrix_ps(_a,_nh,_ca);
// 
// 	_ppu = Matrix_pu(_a,_b,_nh,_cp);
// 	_pvu = Matrix_pu(_a,_b,_nh,_cv);
// 	_pau = Matrix_pu(_a,_b,_nh,_ca);	
	
//	_footx_real.setZero(_footstepsnumber);  _footy_real.setZero(_footstepsnumber); _footz_real.setZero(_footstepsnumber);
//	_footx_real_next.setZero(_nsum);  _footy_real_next.setZero(_nsum); _footz_real_next.setZero(_nsum);
//	_footx_real_next1.setZero(_nsum);  _footy_real_next1.setZero(_nsum); _footz_real_next1.setZero(_nsum);
	

// 	_fx.setZero(1);
// 	_fy.setZero(1);

// 	_fxx_global.setZero(1);
// 	_fyy_global.setZero(1);	
	
	
	/// zmp-constraints	
//	_zmpx_ub.setConstant(_nsum,0.07);  _zmpx_lb.setConstant(_nsum,-0.03);
//	_zmpy_ub.setConstant(_nsum,0.05); _zmpy_lb.setConstant(_nsum,-0.05);
		
	
	// com-support range
// 	_comx_max.setConstant(1,0.06);
// 	_comx_min.setConstant(1,-0.04);  
// 	_comy_max.setConstant(1,0.6);  
// 	_comy_min.setConstant(1,0.02);
	
	// angle range
// 	_thetax_max.setConstant(1,10*M_PI/180);  
// 	_thetax_min.setConstant(1,-5*M_PI/180);
// 	_thetay_max.setConstant(1,10*M_PI/180);  
// 	_thetay_min.setConstant(1,-10*M_PI/180);
	
	// torque range
// 	_torquex_max.setConstant(1,80/_j_ini); 
// 	_torquex_min.setConstant(1,-60/_j_ini);
// 	_torquey_max.setConstant(1,80/_j_ini);  
// 	_torquey_min.setConstant(1,-80/_j_ini);	

	
	
	cout << "finish!!!!!!!!!!! initial for MPC parameters"<<endl;
///===========initiallize: preparation for MPC solution	
	// sulotion preparation		
//	_V_optimal.setZero(_Nt, _nsum);	
	
// 	 store n_vis
// 	_flag.setZero(_nsum);	
// 	_flag_global.setZero(_nsum);	
	
//	_Lx_ref.setZero(_nstep); _Ly_ref.setZero(_nstep); _Lz_ref.setZero(_nstep);
//	_V_ini.setZero(_Nt,1);
// 	_comx_center_ref.setZero(_nh,1);
// 	_comy_center_ref.setZero(_nh,1);
// 	_comz_center_ref.setZero(_nh,1);
	
/*	_thetax_center_ref.setZero(_nh,1); 
	_thetay_center_ref.setZero(_nh,1);*/	
	
        

	
// 	parameters for objective function======================	
// 	 _Rx = 1;           _Ry = 1;            _Rz =1;
// 	_alphax = 1;       _alphay = 1;        _alphaz = 10; 
// 	_beltax = 5000;   _beltay = 10;     _beltaz = 20000000;
// 	_gamax =  10000000; _gamay = 10000000;  _gamaz = 200;
// 	_Rthetax = 1; _Rthetay = 1;
// 	_alphathetax = 100; _alphathetay = 100;
// 	_beltathetax = 500000; _beltathetay = 500000;
	
//	_pvu_2 = _pvu.transpose()*_pvu;
//	_ppu_2 = _ppu.transpose()*_ppu;

	
//	_loop = 2;

///////////////////////////////////////////////////////////////
//////////// next code block just run once	
/*	  A_unit.setIdentity(_nh,_nh);
	  C_unit.setIdentity(_nstep,_nstep);*/		
	

	  // optimization objective function 	
/*	  _WX.setZero(_nh,_nh);
	  _WY.setZero(_nh,_nh);
	  _WZ.setZero(_nh,_nh);
	  _WthetaX.setZero(_nh,_nh);
	  _WthetaY.setZero(_nh,_nh);
	  _PHIX.setZero(_nstep,_nstep);
	  _PHIY.setZero(_nstep,_nstep);
	  _PHIZ.setZero(_nstep,_nstep);
	  _Q_goal.setZero(_Nt,_Nt);
	  _q_goal.setZero(_Nt,1);
	  _Q_goal1.setZero(_Nt,_Nt);
	  _q_goal1.setZero(_Nt,1);	
	  
	  _WX = _Rx*0.5 * A_unit + _alphax*0.5 * _pvu_2 + _beltax*0.5 * _ppu_2;	  
	  _WY = _Ry/2 * A_unit + _alphay/2 * _pvu_2 + _beltay/2 * _ppu_2;
	  _WZ = _Rz/2 * A_unit + _alphaz/2 * _pvu_2 + _beltaz/2 * _ppu_2;  
	  _WthetaX = _Rthetax/2 * A_unit + _alphathetax/2 * _pvu_2 + _beltathetax/2 * _ppu_2;
	  _WthetaY = _Rthetay/2 * A_unit + _alphathetay/2 * _pvu_2 + _beltathetay/2 * _ppu_2;
	  _PHIX  = _gamax/2 * C_unit;
	  _PHIY  = _gamay/2 * C_unit;
	  _PHIZ  = _gamaz/2 * C_unit;
		  
	  _Q_goal.block<_nh, _nh>(0, 0) = _WX;
	  _Q_goal.block< Dynamic, Dynamic>(_nh, _nh,_nh, _nh) = _WY;
	  _Q_goal.block< Dynamic, Dynamic>(2*_nh, 2*_nh,_nh, _nh) = _WZ;
	  _Q_goal.block< Dynamic, Dynamic>(3*_nh, 3*_nh,_nh, _nh) = _WthetaX;
	  _Q_goal.block< Dynamic, Dynamic>(4*_nh, 4*_nh,_nh, _nh) = _WthetaY;
	  _Q_goal.block< Dynamic, Dynamic>(5*_nh, 5*_nh,_nstep,_nstep) = _PHIX;
	  _Q_goal.block< Dynamic, Dynamic>(5*_nh+_nstep, 5*_nh+_nstep,_nstep, _nstep) = _PHIY;
	  _Q_goal.block< Dynamic, Dynamic>(5*_nh+2*_nstep, 5*_nh+2*_nstep,_nstep, _nstep) = _PHIZ;	  
	
	  _Q_goal1 = 2 * _Q_goal;*/	
	
	
	  // constraints
// 	  _Sjx.setZero(_nh,_Nt);
// 	  _Sjy.setZero(_nh,_Nt);
// 	  _Sjz.setZero(_nh,_Nt);
// 	  _Sjthetax.setZero(_nh,_Nt);
// 	  _Sjthetay.setZero(_nh,_Nt);
// 	  _Sjx.block<_nh, _nh>(0, 0) = A_unit;
// 	  _Sjy.block<_nh, _nh>(0, _nh) = A_unit;
// 	  _Sjz.block<_nh, _nh>(0, 2*_nh) = A_unit;	
// 	  _Sjthetax.block<_nh, _nh>(0, 3*_nh) = A_unit;
// 	  _Sjthetay.block<_nh, _nh>(0, 4*_nh) = A_unit;
	  
	  // ZMP boundary preparation
// 	  _H_q_upx.setZero(_nh,_Nt);
// 	  _F_zmp_upx.setZero(_nh,1);
// 	  _H_q_lowx.setZero(_nh,_Nt);
// 	  _F_zmp_lowx.setZero(_nh,1);
// 	  _H_q_upy.setZero(_nh,_Nt);
// 	  _F_zmp_upy.setZero(_nh,1);
// 	  _H_q_lowy.setZero(_nh,_Nt);
// 	  _F_zmp_lowy.setZero(_nh,1);

/*	  _phi_i_x_up.setZero(_Nt,_Nt);
	  _p_i_x_t_up.setZero(_Nt,_nh);
	  _del_i_x_up.setZero(1,_nh);
	  _phi_i_x_low.setZero(_Nt,_Nt);
	  _p_i_x_t_low.setZero(_Nt,_nh);
	  _del_i_x_low.setZero(1,_nh);
	  _phi_i_y_up.setZero(_Nt,_Nt);
	  _p_i_y_t_up.setZero(_Nt,_nh);
	  _del_i_y_up.setZero(1,_nh);
	  _phi_i_y_low.setZero(_Nt,_Nt);
	  _p_i_y_t_low.setZero(_Nt,_nh);
	  _del_i_y_low.setZero(1,_nh);	*/  
	  
	  // angle boundary preparation
/*	  _q_upx.setZero(_nh,_Nt);
	  _qq_upx.setZero(_nh,1);
	  _q_lowx.setZero(_nh,_Nt);
	  _qq_lowx.setZero(_nh,1);
	  _q_upy.setZero(_nh,_Nt);
	  _qq_upy.setZero(_nh,1);
	  _q_lowy.setZero(_nh,_Nt);
	  _qq_lowy.setZero(_nh,1);
	  
	  _qq1_upx.setZero(_nh,1);
	  _qq1_lowx.setZero(_nh,1);
	  _qq1_upy.setZero(_nh,1);
	  _qq1_lowy.setZero(_nh,1);*/	  

	  // torque bondary preparation
// 	  _t_upx.setZero(_nh,_Nt);
// 	  _tt_upx.setZero(_nh,1);
// 	  _t_lowx.setZero(_nh,_Nt);
// 	  _tt_lowx.setZero(_nh,1);
// 	  _t_upy.setZero(_nh,_Nt);
// 	  _tt_upy.setZero(_nh,1);
// 	  _t_lowy.setZero(_nh,_Nt);
// 	  _tt_lowy.setZero(_nh,1);
// 	  
// 	  _tt1_upx.setZero(_nh,1);
// 	  _tt1_lowx.setZero(_nh,1);
// 	  _tt1_upy.setZero(_nh,1);
// 	  _tt1_lowy.setZero(_nh,1);
	  
	  // CoM height boundary preparation
// 	  _H_h_upz.setZero(_nh,_Nt);
// 	  _F_h_upz.setZero(_nh,1);
// 	  _H_h_lowz.setZero(_nh,_Nt);
// 	  _F_h_lowz.setZero(_nh,1);
// 	  _delta_footz_up.setZero(_nh,1);
// 	  _delta_footz_low.setZero(_nh,1);

	  // CoM height acceleration boundary preparation
/*	  _H_hacc_lowz.setZero(_nh,_Nt);
	  _F_hacc_lowz.setZero(_nh,1);
	  _delta_footzacc_up.setZero(_nh,1);*/	  

	  
	  //swing foot velocity constraints
// 	  _Footvx_max.setZero(1,_Nt);
// 	  _Footvx_min.setZero(1,_Nt);
// 	  _Footvy_max.setZero(1,_Nt);
// 	  _Footvy_min.setZero(1,_Nt);
// 	  _footubxv.setZero(1,1);
// 	  _footlbxv.setZero(1,1);
// 	  _footubyv.setZero(1,1);
// 	  _footlbyv.setZero(1,1);
	  
	  // foot vertical location-equality constraints
// 	  _H_q_footz.setZero(1, _Nt);
// 	  _F_footz.setZero(1, 1);
	  
	  // CoMZ height-equality constraints
/*	  _h_h.setZero(_nh, _Nt);*/	  

	  // body inclination-equality constraints
// 	  _a_hx.setZero(_nh, _Nt);
// 	  _a_hxx.setZero(_nh, 1);
// 	  _a_hy.setZero(_nh, _Nt);
// 	  _a_hyy.setZero(_nh, 1);
	  
	  
	  // foot location constraints_ki, _k_yu
// 	  _Sfoot.setZero(1,2);
// 	  _Sfoot(0,0) = -1;
// 	  _Sfoot(0,1) = 1;
	  
/*	  _S1.setZero(1,_nh);
	  _S1(0,0) = 1;	*/  
	  
          
        cout << "Before initialization for ZMP constraints"<<endl;  
// 	// offline calulated the ZMP constraints coefficient: ERROR check check check!!!!!!!!!!
// 	vector <Eigen::MatrixXd> x_offline1(_nh)  ;
// 	  
// 	  
// 	ZMPx_constraints_offfline = x_offline1;
// // 	vector <Eigen::MatrixXd> ZMPy_constraints_offfline(_nh);
// 	ZMPy_constraints_offfline = x_offline1;
// 	
// 	
// // 	vector <Eigen::MatrixXd> ZMPx_constraints_half(_nh);
// // 	vector <Eigen::MatrixXd> ZMPy_constraints_half(_nh);
// 	ZMPx_constraints_half = x_offline1;
// 	
// 	ZMPy_constraints_half = x_offline1;
// 	
// 	
// 	for(int jxx=1; jxx<=_nh; jxx++)
// 	{
// 	  _Si.setZero(1,_nh);
// 	  _Si(0,jxx-1) = 1;
// 	  // ZMP constraints	      
// 		 
//          ZMPx_constraints_offfline[jxx-1] = (_Si * _ppu * _Sjx).transpose() * _Si * _pau * _Sjz - (_Si * _pau * _Sjx).transpose() * _Si * _ppu * _Sjz;
// 	 ZMPx_constraints_half[jxx-1] = - (_Si).transpose() * _Si * _pau * _Sjz;
// 	  
// 	  
//          ZMPy_constraints_offfline[jxx-1] = (_Si * _ppu * _Sjy).transpose() * _Si * _pau * _Sjz - (_Si * _pau * _Sjy).transpose() * _Si * _ppu * _Sjz;
// 	 ZMPy_constraints_half[jxx-1] = - (_Si).transpose() * _Si * _pau * _Sjz;
// 	      
// 	}

	cout << "After initialization for ZMP constraints"<<endl; 
   
	_nTdx =0;
	
	//////polynomial intepolation for lower level interpolation
	_AAA_inv.setZero();
	
	
	_j_count = 0;
	
	
	/////////////for ZMP distribution
	_F_R(2) = _F_L(2) = 0.5 * _robot_mass * _GGG;
	_Co_L.setZero(); _Co_R.setZero();
	
	_comxyzx.setZero(); _comvxyzx.setZero(); _comaxyzx.setZero(); 
	_thetaxyx.setZero(); _thetavxyx.setZero(); _thetaaxyx.setZero();
	_Lfootxyzx.setZero(); _Rfootxyzx.setZero();
	_ZMPxy_realx.setZero();
	
	_comxyzx(2) = _HCOM;		
	
        cout << "Finish initialization of WHOLE MPCCLASS"<<endl; 
  
}




///// main function for optimization
void NLPClass::step_timing_opti_loop(int i,Eigen::Matrix<double,18,1> estimated_state, Eigen::Vector3d _Rfoot_location_feedback, Eigen::Vector3d _Lfoot_location_feedback,double lamda, bool _stopwalking)
{
//   cout <<"i:"<<i<<endl;
  
//  clock_t _t_start,_t_finish;
//  _t_start = clock();
  
  //// step cycle number when (i+1)*dt fall into: attention that _tstep+1 falls ionto the next cycle      
  Indexfind((i+1)*_dt,xyz0);	   /// next one sampling time
  _periond_i = _j_period+1;      ///coincident with Matlab
  _j_period = 0;  
  
  
  
  ////================================================================
  /// judge if stop walking enable: if (from _bjx1 step, reference step length and step height will be set to be zero)
  if(_stopwalking)  
  {    
    for (int i_t = _bjx1; i_t < _footstepsnumber; i_t++) {
      if (i_t == _bjx1)
      {
	_steplength(i_t) = _steplength(i_t)/2;	
	_Lxx_ref(i_t) = _Lxx_ref(i_t)/2;	
      }
      else
      {
      _steplength(i_t) = 0;	
      _Lxx_ref(i_t) = 0;	
      }

      _footx_ref(i_t) = _footx_ref(i_t-1) + _steplength(i_t-1); 	      
    }	  
  }  
  
  
//   cout <<"k:"<<_periond_i<<endl;	


  //ZMP & ZMPv
  _px(0,0) = _footx_ref(_periond_i-1,0); //_zmpvx(0,i) = 0; 
  _py(0,0) = _footy_ref(_periond_i-1,0); //_zmpvy(0,i) = 0; 
  _zmpx_real(0,i) = _px(0,0);
  _zmpy_real(0,i) = _py(0,0);   
  
//   cout<<"_py:"<<_py(0,0)<<endl;  
  
  ///remaining step timing for the next stepping
  _ki = round(_tx(_periond_i-1,0)/_dt);
  _k_yu = i-_ki;                
  _Tk = _ts(_periond_i-1) - _k_yu*_dt; 
  
  //// reference remaining time and step length and step width
//    position tracking 
//   _Lxx_refx = _footx_offline(_periond_i)-_footx_ref(_periond_i-1);        
//   _Lyy_refy = _footy_offline(_periond_i)-_footy_ref(_periond_i-1); //%% tracking the step location
//   _Lxx_refx1 = _footx_offline(_periond_i+1)-_footx_offline(_periond_i);        
//   _Lyy_refy1 = _footy_offline(_periond_i+1)-_footy_offline(_periond_i); // tracking the step location

  //    velocity tracking 
   _Lxx_refx = _Lxx_ref(_periond_i-1);        
   _Lyy_refy = _Lyy_ref(_periond_i-1); //%% tracking relative location
   _Lxx_refx1 = _Lxx_ref(_periond_i);        
   _Lyy_refy1 = _Lyy_ref(_periond_i); // tracking relative location
   
  
  
   
  
  _tr1_ref = cosh(_Wn*_Tk);        _tr2_ref = sinh(_Wn*_Tk);
  _tr1_ref1 = cosh(_Wn*_ts(_periond_i));  _tr2_ref1 = sinh(_Wn*_ts(_periond_i));  
  

  
  // warm start
  if (i==1)
  {
    _vari_ini << _Lxx_refx,
                 _Lyy_refy,
		 _tr1_ref,
		 _tr2_ref,
		 _Lxx_refx1,
                 _Lyy_refy1,
		 _tr1_ref1,
		 _tr2_ref1;
  }
  else
  {
    _vari_ini = _Vari_ini.col(i-1);
  }
  
  
// step timing upper&lower boundaries  modification
  if ((_t_min -_k_yu*_dt)>=0.0001)
  {
    _tr1_min = cosh(_Wn*(_t_min-_k_yu*_dt));
    _tr2_min = sinh(_Wn*(_t_min-_k_yu*_dt));    
  }
  else
  {
    _tr1_min = cosh(_Wn*(0.0001));
    _tr2_min = sinh(_Wn*(0.0001));      
  }
 
  
 
  _tr1_max = cosh(_Wn*(_t_max-_k_yu*_dt));
  _tr2_max = sinh(_Wn*(_t_max-_k_yu*_dt));
  
  _tr11_min = cosh(_Wn*_t_min);
  _tr21_min = sinh(_Wn*_t_min);     
  _tr11_max = cosh(_Wn*_t_max);
  _tr21_max = sinh(_Wn*_t_max);    
    
//   cout <<_tr11_min<<endl;
//   cout <<_tr21_min<<endl;   
//   cout <<_tr11_max<<endl;
//   cout <<_tr21_max<<endl;  

  if (i==1)
  {
        _COMx_is(_periond_i-1) = _comx_feed(i-1)-_footx_ref(_periond_i-1);  
	_COMx_es.col(_periond_i-1) = _SS1*_vari_ini*0.5;  
	_COMvx_is(_periond_i-1)= (_COMx_es(_periond_i-1)-_COMx_is(_periond_i-1)*_SS3*_vari_ini)/(1/_Wn *_SS4*_vari_ini);
        _COMy_is(_periond_i-1) = _comy_feed(i-1)-_footy_ref(_periond_i-1);  
	_COMy_es.col(_periond_i-1) = _SS2*_vari_ini*0.5;  
	_COMvy_is(_periond_i-1)= (_COMy_es(_periond_i-1)-_COMy_is(_periond_i-1)*_SS3*_vari_ini)/(1/_Wn *_SS4*_vari_ini);  
        _comvx_endref= _Wn*_COMx_is(_periond_i-1)*_SS4*_vari_ini + _COMvx_is(_periond_i-1)*_SS3*_vari_ini;
        _comvy_endref= _Wn*_COMy_is(_periond_i-1)*_SS4*_vari_ini + _COMvy_is(_periond_i-1)*_SS3*_vari_ini;     
  }

  
// SEQUENCE QUADARTIC PROGRAMMING-step timing &step location optimization
  for (int xxxx=1; xxxx<=3; xxxx++)
  { 
    //// be careful the  divide / (one of the factor should be double: type)
// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// // %% optimal programme formulation/OBJECTIVE FUNCTION:
// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
// // %%%%% Lx0 = Lxk(:,i);Ly0 = Lyk(:,i);tr10 = Tr1k(:,i); tr20 = Tr2k(:,i);   
    step_timing_object_function(i);

// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// // %% constraints
// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    step_timing_constraints(i);
    
    ///// generated CoM final position
    if (_method_flag_nlp==2)
    {
      ////number of inequality constraints: 8+8+8+4+4+4
      solve_stepping_timing_twosteps();      
      
    }
    else
    {
      if (_method_flag_nlp==1)
      {
	solve_stepping_timing_onestep();
      }
      else
      {	
	solve_stepping_timing_non();
      }
      
    }

    

    
    if (_X.rows() == 8)
    {
      _vari_ini += _X;
    }   
    
    
  }
  
  
// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// // %% results postprocession
// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  _Vari_ini.col(i) = _vari_ini;  
  
  
  
// update the optimal parameter in the real-time: at the result, the effect of optimization reduced gradually
  /// exception catch: for 
  if(_robot_name == "coman")
  {
    if (_vari_ini(1)>=0.2)
    {
      _vari_ini(1) =0.2;
    }  
    else
    {
      if (_vari_ini(1)<=-0.2)
      {
	_vari_ini(1) =-0.2;
      }    
      
    }
    
    if (_vari_ini(0)>=0.35)
    {
      _vari_ini(0) =0.35;
    }  
    else
    {
      if (_vari_ini(0)<=-0.05)
      {
	_vari_ini(0) =-0.05;
      }    
      
    }    
    
  }
  else
  {
    if (_vari_ini(1)>=0.28)
    {
      _vari_ini(1) =0.28;
    }  
    else
    {
      if (_vari_ini(1)<=-0.28)
      {
	_vari_ini(1) =-0.28;
      }    
      
    }
    
    if (_vari_ini(0)>=0.35)
    {
      _vari_ini(0) =0.35;
    }  
    else
    {
      if (_vari_ini(0)<=-0.1)
      {
	_vari_ini(0) =-0.1;
      }    
      
    }     
  }
  
  
    
  
  _Lxx_ref(_periond_i-1) = _SS1*_vari_ini;
  _Lyy_ref(_periond_i-1) = _SS2*_vari_ini;
  _ts(_periond_i-1) = _k_yu*_dt+ log((_SS3+_SS4)*_vari_ini)/_Wn;   //check log function
  
  _Lxx_ref(_periond_i) = _SS5*_vari_ini;
  _Lyy_ref(_periond_i) = _SS6*_vari_ini;
  _ts(_periond_i) = log((_SS7+_SS8)*_vari_ini)/_Wn;    
  
  
  _Lxx_ref_real(i) = _Lxx_ref(_periond_i-1);
  _Lyy_ref_real(i) = _Lyy_ref(_periond_i-1);
  _Ts_ref_real(i) = _ts(_periond_i-1);  
  
  

  _COMx_is(_periond_i-1) = _comx_feed(i-1)-_footx_ref(_periond_i-1);  
  _COMx_es.col(_periond_i-1) = _SS1*_vari_ini*0.5;  
  _COMvx_is(_periond_i-1)= (_COMx_es(_periond_i-1)-_COMx_is(_periond_i-1)*_SS3*_vari_ini)/(1/_Wn *_SS4*_vari_ini);
  _COMy_is(_periond_i-1) = _comy_feed(i-1)-_footy_ref(_periond_i-1);  
  _COMy_es.col(_periond_i-1) = _SS2*_vari_ini*0.5;  
  _COMvy_is(_periond_i-1)= (_COMy_es(_periond_i-1)-_COMy_is(_periond_i-1)*_SS3*_vari_ini)/(1/_Wn *_SS4*_vari_ini);    


  
  
// update walking period and step location 
  for (int jxx = _periond_i+1; jxx<=_footstepsnumber; jxx++)
  {
    _tx(jxx-1) = _tx(jxx-2)+_ts(jxx-2);
  }
  
//   the real-time calcuated next one step location: similar with the the footx_real calulated in the NMPC ( _footx_real.row(_bjxx) = _V_ini.row(5*_nh);):
  _footx_ref(_periond_i)=_footx_ref(_periond_i-1)+_SS1*_vari_ini;
  _footy_ref(_periond_i)=_footy_ref(_periond_i-1)+_SS2*_vari_ini;    
  _footx_ref(_periond_i+1)=_footx_ref(_periond_i)+_SS5*_vari_ini;
  _footy_ref(_periond_i+1)=_footy_ref(_periond_i)+_SS6*_vari_ini;   
  
  
  _comvx_endref= _Wn*_COMx_is(_periond_i-1)*_SS4*_vari_ini + _COMvx_is(_periond_i-1)*_SS3*_vari_ini;
  _comvy_endref= _Wn*_COMy_is(_periond_i-1)*_SS4*_vari_ini + _COMvy_is(_periond_i-1)*_SS3*_vari_ini;     
      
// update CoM state  
   
  _comx(i)= _COMx_is(_periond_i-1)*cosh(_Wndt) + _COMvx_is(_periond_i-1)*1/_Wn*sinh(_Wndt)+_px(0);
  _comy(i)= _COMy_is(_periond_i-1)*cosh(_Wndt) + _COMvy_is(_periond_i-1)*1/_Wn*sinh(_Wndt)+_py(0);           
  _comvx(i)= _Wn*_COMx_is(_periond_i-1)*sinh(_Wndt) + _COMvx_is(_periond_i-1)*cosh(_Wndt);
  _comvy(i)= _Wn*_COMy_is(_periond_i-1)*sinh(_Wndt) + _COMvy_is(_periond_i-1)*cosh(_Wndt);     
  _comax(i)= pow(_Wn,2)*_COMx_is(_periond_i-1)*cosh(_Wndt) + _COMvx_is(_periond_i-1)*_Wn*sinh(_Wndt);
  _comay(i)= pow(_Wn,2)*_COMy_is(_periond_i-1)*cosh(_Wndt) + _COMvy_is(_periond_i-1)*_Wn*sinh(_Wndt);   
  
   _comz(i) =_hcom; 
   _comvz(i) =0; 
   _comaz(i) =0;
   


//  external disturbances!!! using feedback data:
 
  /// /// relative state to the actual foot lcoation: very good
  if (_periond_i % 2 == 0)  // odd : left support
  {
    estimated_state(0,0) =  estimated_state(0,0) - _Lfoot_location_feedback(0); //relative comx
    estimated_state(3,0) =  estimated_state(3,0) - _Lfoot_location_feedback(1);	// relative comy 	    
  }
  else
  {
    estimated_state(0,0) =  estimated_state(0,0) - _Rfoot_location_feedback(0);
    estimated_state(3,0) =  estimated_state(3,0) - _Rfoot_location_feedback(1);	  
  }  
  
 // feedback for t=0.6s, stepwidh = HALF_HIP_WIDTH()*2
//   double _lamda_comx = 0.75;
//   double _lamda_comvx = 0.75;
//   double _lamda_comy = 0.99;
//   double _lamda_comvy = 0.99;   
//   
//   if (_periond_i>2)
//   {
//     _lamda_comx = 0.01;
//     _lamda_comvx = 0.01;      
//     _lamda_comy = 0.79;
//     _lamda_comvy = 0.79;     
//   }


  double _lamda_comx = 1;
  double _lamda_comvx = 1;
  double _lamda_comy = 1;
  double _lamda_comvy = 1;   
  
//   if (_periond_i>2)
//   {
//     _lamda_comx = 0.01;
//     _lamda_comvx = 0.01;      
//     _lamda_comy = 0.85;
//     _lamda_comvy = 0.85;     
//   }

 
  _comx_feed(i) = (_lamda_comx*(_comx(i)-_px(0))+(1-_lamda_comx)*estimated_state(0,0))+_px(0);    
  _comvx_feed(i) = _lamda_comvx*_comvx(i) + (1-_lamda_comvx)*estimated_state(1,0);
  _comax_feed(i) = _lamda_comx*_comax(i) + (1-_lamda_comx)*estimated_state(2,0);
  _comy_feed(i) = (_lamda_comy*(_comy(i)-_py(0))+(1-_lamda_comy)*estimated_state(3,0))+_py(0);    
  _comvy_feed(i) = _lamda_comvy*_comvy(i) + (1-_lamda_comvy)*estimated_state(4,0);  
  _comay_feed(i) = _lamda_comy*_comay(i) + (1-_lamda_comy)*estimated_state(5,0);  
  


  
 // _t_finish = clock();  
  
  
  
  
  _t_f.setLinSpaced(_nh,(i+1)*_dt, (i+_nh)*_dt);
  
  Indexfind(i*_dt,xyz1);                   //// step cycle number when (i)*dt fall into : current sampling time
  _bjxx = _j_period+1;  //coincidence with matlab 
  _j_period = 0;  
  
  Indexfind(_t_f(0),xyz1);                /// step cycle number when (i+1)*dt fall into : current sampling time
  _bjx1 = _j_period+1;
  _j_period = 0;  
  
  _td = 0.2* _ts;
  
  _footxyz_real.row(0) = _footx_ref.transpose();
  _footxyz_real.row(1) = _footy_ref.transpose();  
  _footxyz_real.row(2) = _footz_ref.transpose(); 
  
  ///////// real_zmpx: zmpy generation during the next double support phase: 
  int t_d_period = 0;
  
  _nTdx = round(_td(_bjx1-1)/_dt)+2;  
  for (int jxx=2; jxx <=_nTdx; jxx++)
  {
    Indexfind(_t_f(jxx-1),xyz1);                
    t_d_period = _j_period+1;
    _j_period = 0;   
    _zmpx_real(0,i+jxx-1) = _footx_ref(t_d_period-1,0);
    _zmpy_real(0,i+jxx-1) = _footy_ref(t_d_period-1,0); 
   
  }
  

  
  
  
  
}



int NLPClass::Indexfind(double goalvari, int xyz)
{
  _j_period = 0;
  
  if (xyz<-0.5)
  {
      while (goalvari > (_tx(_j_period))+0.0001)
      {
	_j_period++;
      }    
      _j_period = _j_period-1;	    
  }
  else
  {
  
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

}


void NLPClass::step_timing_object_function(int i)
{
    _AxO(0) = _comx_feed(i-1)-_footx_ref(_periond_i-1); _BxO(0) = _comvx_feed(i-1)/_Wn; _Cx(0,0) =-0.5*_Lxx_refx; 
    _Axv = _Wn*_AxO;  _Bxv = _Wn*_BxO; _Cxv=_comvx_endref; 
    _AyO(0) = _comy_feed(i-1)-_footy_ref(_periond_i-1); _ByO(0) = _comvy_feed(i-1)/_Wn; _Cy(0,0) =-0.5*_Lyy_refy; 
    _Ayv = _Wn*_AyO;  _Byv = _Wn*_ByO; _Cyv=_comvy_endref;    
    
    _SQ_goal0(0,0)=0.5*_bbx;
    _SQ_goal0(1,1)=0.5*_bby;
    _SQ_goal0(2,2)=0.5*(_rr1+_aax*_AxO(0,0)*_AxO(0,0)+_aay*_AyO(0,0)*_AyO(0,0)+_aaxv*_Axv(0,0)*_Axv(0,0)+_aayv*_Ayv(0,0)*_Ayv(0,0));
    _SQ_goal0(2,3)=0.5*(     _aax*_AxO(0,0)*_BxO(0,0)+_aay*_AyO(0,0)*_ByO(0,0)+_aaxv*_Axv(0,0)*_Bxv(0,0)+_aayv*_Ayv(0,0)*_Byv(0,0));
    _SQ_goal0(3,2)=0.5*(     _aax*_BxO(0,0)*_AxO(0,0)+_aay*_ByO(0,0)*_AyO(0,0)+_aaxv*_Bxv(0,0)*_Axv(0,0)+_aayv*_Byv(0,0)*_Ayv(0,0));
    _SQ_goal0(3,3)=0.5*(_rr2+_aax*_BxO(0,0)*_BxO(0,0)+_aay*_ByO(0,0)*_ByO(0,0)+_aaxv*_Bxv(0,0)*_Bxv(0,0)+_aayv*_Byv(0,0)*_Byv(0,0));    
    _SQ_goal0(4,4)=0.5*_bbx1; 
    _SQ_goal0(5,5)=0.5*_bby1; 
    _SQ_goal0(6,6)=0.5*_rr11; 
    _SQ_goal0(7,7)=0.5*_rr21; 
    
    _SQ_goal = (_SQ_goal0+_SQ_goal0.transpose())/2.0;  
    _Sq_goal << -_bbx*_Lxx_refx,
	      -_bby*_Lyy_refy,
	      -_rr1*_tr1_ref+_aax*_AxO(0,0)*_Cx(0,0)+_aay*_AyO(0,0)*_Cy(0,0)+_aaxv*_Axv(0,0)*_Cxv(0,0)+_aayv*_Ayv(0,0)*_Cyv(0,0),
	      -_rr2*_tr2_ref+_aax*_BxO(0,0)*_Cx(0,0)+_aay*_ByO(0,0)*_Cy(0,0)+_aaxv*_Bxv(0,0)*_Cxv(0,0)+_aayv*_Byv(0,0)*_Cyv(0,0),
	      -_bbx1*_Lxx_refx1,
	      -_bby1*_Lyy_refy1,
	      -_rr11*_tr1_ref1,
	      -_rr21*_tr2_ref1;   
	      
    _SQ_goal1 = 2 * _SQ_goal;
    _Sq_goal1 = 2 * _SQ_goal * _vari_ini + _Sq_goal;
    
    
    
    _Ax= _SS3.transpose()*_AxO*_SS7+_SS4.transpose()*_BxO*_SS7-_SS1.transpose()*_SS7+_SS4.transpose()*_AxO*_SS8+_SS3.transpose()*_BxO*_SS8;   
    _Bx(0,0)= -0.5*_Lxx_refx1; 
    _Ay= _SS3.transpose()*_AyO*_SS7+_SS4.transpose()*_ByO*_SS7-_SS2.transpose()*_SS7+_SS4.transpose()*_AyO*_SS8+_SS3.transpose()*_ByO*_SS8;   
    _By(0,0)= -0.5*_Lyy_refy1;         ///check!!!!!!!!

    _ixi = _vari_ini.transpose()*_Ax*_vari_ini;
    _iyi = _vari_ini.transpose()*_Ay*_vari_ini;
    _SQ_goal20= _aax1/2.0*2*( 2*(_ixi(0,0)*_Ax.transpose() + 2*_Ax*_vari_ini*(_Ax*_vari_ini).transpose()) + 2*(2*_Bx(0,0)*_Ax.transpose()))+_aay1/2.0*2*( 2*(_iyi(0,0)*_Ay.transpose() + 2*_Ay*_vari_ini*(_Ay*_vari_ini).transpose()) + 2*(2*_By(0,0)*_Ay.transpose()));

    _SQ_goal2 = (_SQ_goal20.transpose()+_SQ_goal20)/2.0;    
    _Sq_goal2= _aax1/2.0*2*(2*_Ax*_vari_ini)*(_vari_ini.transpose()*_Ax*_vari_ini+_Bx)  +  _aay1/2.0*2*(2*_Ay*_vari_ini)*(_vari_ini.transpose()*_Ay*_vari_ini+_By);    
    
    _SQ_goal3 = _SQ_goal1+_SQ_goal2;
    _Sq_goal3 = _Sq_goal1+_Sq_goal2;   
  

}

void NLPClass::step_timing_constraints(int i)
{
  
// // %% constraints
// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// // %% remaining time constraints: inequality constraints    
    _trx1_up = _SS3;
    _det_trx1_up = -(_SS3*_vari_ini);
    _det_trx1_up(0,0) +=  _tr1_max; 

    _trx1_lp = -_SS3;
    _det_trx1_lp = -(-_SS3*_vari_ini); 
    _det_trx1_lp(0,0) -= _tr1_min;  	

    _trx2_up = _SS4;
    _det_trx2_up = -(_SS4*_vari_ini);
    _det_trx2_up(0,0) += _tr2_max;	

    _trx2_lp = -_SS4;
    _det_trx2_lp = -(-_SS4*_vari_ini);  
    _det_trx2_lp(0,0) -= _tr2_min;  	

    _trx3_up = _SS7;
    _det_trx3_up = -(_SS7*_vari_ini);
    _det_trx3_up(0,0) += _tr11_max;	

    _trx3_lp = -_SS7;
    _det_trx3_lp = -(-_SS7*_vari_ini);        
    _det_trx3_lp(0,0) -= _tr11_min;    
    
    _trx4_up = _SS8;
    _det_trx4_up = -(_SS8*_vari_ini);
    _det_trx4_up(0,0) += _tr21_max;
    
    _trx4_lp = -_SS8;
    _det_trx4_lp = -(-_SS8*_vari_ini);      
    _det_trx4_lp(0,0) -= _tr21_min;  
    
    _trx.row(0) = _trx1_up; _trx.row(1) = _trx1_lp; _trx.row(2) = _trx2_up; _trx.row(3) = _trx2_lp;
    _trx.row(4) = _trx3_up; _trx.row(5) = _trx3_lp; _trx.row(6) = _trx4_up; _trx.row(7) = _trx4_lp;	
    
    _det_trx.row(0) = _det_trx1_up; _det_trx.row(1) = _det_trx1_lp; _det_trx.row(2) = _det_trx2_up; _det_trx.row(3) = _det_trx2_lp;
    _det_trx.row(4) = _det_trx3_up; _det_trx.row(5) = _det_trx3_lp; _det_trx.row(6) = _det_trx4_up; _det_trx.row(7) = _det_trx4_lp;		
    
     
  
  
// tr1 & tr2: equation constraints:   
    _trx12 = (2*(_SS3.transpose()*_SS3-_SS4.transpose()*_SS4)*_vari_ini).transpose();
    _det_trx12 = -(_vari_ini.transpose()*(_SS3.transpose()*_SS3-_SS4.transpose()*_SS4)*_vari_ini);
    _det_trx12(0,0) +=1; 

    _trx121 = (2*(_SS7.transpose()*_SS7-_SS8.transpose()*_SS8)*_vari_ini).transpose();
    _det_trx121 = -(_vari_ini.transpose()*(_SS7.transpose()*_SS7-_SS8.transpose()*_SS8)*_vari_ini);
    _det_trx121(0,0) +=1; 

    _trxx.row(0) = _trx12;   _trxx.row(1) = _trx121;
    _det_trxx.row(0) = _det_trx12; _det_trxx.row(1) = _det_trx121;    

     
//  step location and step duration equality constraints 
    _Lxx_ref1_e = _SS5;     _Lyy_ref1_e = _SS6;  _tr1_ref1_e = _SS7;  _tr2_ref1_e = _SS8;  
     _Lxx_ref_e = _SS1;      _Lyy_ref_e = _SS2;   _tr1_ref_e = _SS3;   _tr2_ref_e = _SS4;	
    
    _Lxx_ref1_eq = -(_SS5*_vari_ini);   
    _Lyy_ref1_eq = -(_SS6*_vari_ini); 
    _tr1_ref1_eq = -(_SS7*_vari_ini);
    _tr2_ref1_eq = -(_SS8*_vari_ini);
    
    _Lxx_ref_eq = -(_SS1*_vari_ini);   
    _Lyy_ref_eq = -(_SS2*_vari_ini); 
    _tr1_ref_eq = -(_SS3*_vari_ini);
    _tr2_ref_eq = -(_SS4*_vari_ini);    
    
    
    _Lxx_ref1_eq(0,0) += _Lxx_refx1;  
    _Lyy_ref1_eq(0,0) +=_Lyy_refy1; 
    _tr1_ref1_eq(0,0) +=_tr1_ref1; 
    _tr2_ref1_eq(0,0) += _tr2_ref1;
    
    
     _Lxx_ref_eq(0,0) += _Lxx_refx;   
     _Lyy_ref_eq(0,0) +=_Lyy_refy;   
     _tr1_ref_eq(0,0) +=_tr1_ref;   
     _tr2_ref_eq(0,0) += _tr2_ref;    


    
    

//  foot location constraints      
    if (_periond_i % 2 == 0)
    {
      
	  if (i>=(round(2*_ts(1)/_dt))+1) ///update the footy_limit
	  {
	      
            _footy_min=-(0.3); 
	    _footy_max = -(0.11);
	  }
	  else
	  {
	    _footy_min=-(0.28); 
	    _footy_max=-(0.07); 	    
	  }  
    }
    else
    { 
       if (i>=(round(2*_ts(1)/_dt))+1)
       {
	_footy_max= 0.28; 
	_footy_min =0.11;
      }
       else
       {
	_footy_max=0.3; 
	_footy_min=  0.07; 	 
      }     
      
    }
    // only the next one step
    _h_lx_up = _SS1;
    _det_h_lx_up = -(_SS1*_vari_ini);
    _det_h_lx_up(0,0) += _footx_max;

    _h_lx_lp = -_SS1;
    _det_h_lx_lp = -(-_SS1*_vari_ini); 
    _det_h_lx_lp(0,0) -= _footx_min;

    _h_ly_up = _SS2;
    _det_h_ly_up = -(_SS2*_vari_ini);
    _det_h_ly_up(0,0) += _footy_max;

    _h_ly_lp = -_SS2;
    _det_h_ly_lp = -(-_SS2*_vari_ini);
    _det_h_ly_lp(0,0) -= _footy_min;    

    _h_lx_up1 = _SS5;
    _det_h_lx_up1 = -(_SS5*_vari_ini);
    _det_h_lx_up1(0,0) += _footx_max;    

    _h_lx_lp1 = -_SS5;
    _det_h_lx_lp1 = -(-_SS5*_vari_ini);
    _det_h_lx_lp1(0,0) -= _footx_min;     

    _h_ly_up1 = _SS6;
    _det_h_ly_up1 = -(_SS6*_vari_ini);
    _det_h_ly_up1(0,0) -= _footy_min;
    
    _h_ly_lp1 = -_SS6;
    _det_h_ly_lp1 = -(-_SS6*_vari_ini); 
    _det_h_ly_lp1(0,0) += _footy_max;     

    _h_lx_upx.row(0)= _h_lx_up;    _h_lx_upx.row(1)= _h_lx_lp;   _h_lx_upx.row(2)= _h_ly_up;     _h_lx_upx.row(3)= _h_ly_lp;
    _h_lx_upx.row(4)= _h_lx_up1;   _h_lx_upx.row(5)= _h_lx_lp1;  _h_lx_upx.row(6)= _h_ly_up1;    _h_lx_upx.row(7)= _h_ly_lp1;
    _det_h_lx_upx.row(0)=_det_h_lx_up;  _det_h_lx_upx.row(1)=_det_h_lx_lp;  _det_h_lx_upx.row(2)=_det_h_ly_up;  _det_h_lx_upx.row(3)=_det_h_ly_lp;
    _det_h_lx_upx.row(4)=_det_h_lx_up1; _det_h_lx_upx.row(5)=_det_h_lx_lp1; _det_h_lx_upx.row(6)=_det_h_ly_up1; _det_h_lx_upx.row(7)=_det_h_ly_lp1;    

 

    
    
// swing foot velocity boundary
    if (_k_yu ==0)
    {
	_h_lvx_up.setZero();  _h_lvx_lp.setZero(); _h_lvy_up.setZero(); _h_lvy_lp.setZero();
	_h_lvx_up1.setZero(); _h_lvx_lp1.setZero(); _h_lvy_up1.setZero(); _h_lvy_lp1.setZero();
	_det_h_lvx_up.setZero();_det_h_lvx_lp.setZero();_det_h_lvy_up.setZero();_det_h_lvy_lp.setZero();
	_det_h_lvx_up1.setZero();_det_h_lvx_lp1.setZero();_det_h_lvy_up1.setZero();_det_h_lvy_lp1.setZero();    
    }                
    else
    {   
	_h_lvx_up = _SS1;
	_det_h_lvx_up(0,0) = -(_SS1*_vari_ini-_Lxx_ref(_periond_i-1)- _footx_vmax*_dt);

	_h_lvx_lp = -_SS1;
	_det_h_lvx_lp(0,0) = _SS1*_vari_ini-_Lxx_ref(_periond_i-1)-_footx_vmin*_dt; 

	_h_lvy_up = _SS2;
	_det_h_lvy_up(0,0) = -(_SS2*_vari_ini-_Lyy_ref(_periond_i-1)- _footy_vmax*_dt);

	_h_lvy_lp = -_SS2;
	_det_h_lvy_lp(0,0) = _SS2*_vari_ini-_Lyy_ref(_periond_i-1)-_footy_vmin*_dt;  

	
	_h_lvx_up1.setZero(); 
	_h_lvx_lp1.setZero(); 
	_h_lvy_up1.setZero(); 
	_h_lvy_lp1.setZero();	
	
	_det_h_lvx_up1(0,0) = 0.001;

	_det_h_lvx_lp1(0,0) = 0.001; 

	_det_h_lvy_up1(0,0) = 0.001;

	_det_h_lvy_lp1(0,0) = 0.001;  

    }                                   
    _h_lvx_upx.row(0)= _h_lvx_up;    _h_lvx_upx.row(1)= _h_lvx_lp;  _h_lvx_upx.row(2)= _h_lvy_up;   _h_lvx_upx.row(3)= _h_lvy_lp;
    _h_lvx_upx.row(4)= _h_lvx_up1;   _h_lvx_upx.row(5)= _h_lvx_lp1; _h_lvx_upx.row(6)= _h_lvy_up1;  _h_lvx_upx.row(7)= _h_lvy_lp1;
    _det_h_lvx_upx.row(0)=_det_h_lvx_up; _det_h_lvx_upx.row(1)=_det_h_lvx_lp; _det_h_lvx_upx.row(2)=_det_h_lvy_up; _det_h_lvx_upx.row(3)=_det_h_lvy_lp;
    _det_h_lvx_upx.row(4)=_det_h_lvx_up1;_det_h_lvx_upx.row(5)=_det_h_lvx_lp1;_det_h_lvx_upx.row(6)=_det_h_lvy_up1;_det_h_lvx_upx.row(7)=_det_h_lvy_lp1;    
  
    


////////////////////////// CoM position relative to the current support center
// CoM accelearation boundary    
    
    _AA= _Wn*sinh(_Wn*_dt); _CCx = _comx_feed(0,i-1)-_footx_ref(_periond_i-1,0); _BBx = pow(_Wn,2)*_CCx*cosh(_Wn*_dt); 
		            _CCy = _comy_feed(0,i-1)-_footy_ref(_periond_i-1,0); _BBy = pow(_Wn,2)*_CCy*cosh(_Wn*_dt);

    _AA1x = _AA*_Wn; _AA2x = -2* _AA*_CCx*_Wn;  _AA3x = 2*_BBx; 
    _AA1y = _AA*_Wn; _AA2y = -2* _AA*_CCy*_Wn;  _AA3y = 2*_BBy;


    _CoM_lax_up = _AA1x*_SS1+_AA2x*_SS3+(_AA3x-2*_comax_max)*_SS4;
    _det_CoM_lax_up = -(_AA1x*_SS1+_AA2x*_SS3+(_AA3x-2*_comax_max)*_SS4)*_vari_ini;

    _CoM_lax_lp = -_AA1x*_SS1-_AA2x*_SS3-(_AA3x-2*_comax_min)*_SS4;
    _det_CoM_lax_lp = -(-_AA1x*_SS1-_AA2x*_SS3-(_AA3x-2*_comax_min)*_SS4)*_vari_ini; 

    _CoM_lay_up = _AA1y*_SS2+_AA2y*_SS3+(_AA3y-2*_comay_max)*_SS4;
    _det_CoM_lay_up = -(_AA1y*_SS2+_AA2y*_SS3+(_AA3y-2*_comay_max)*_SS4)*_vari_ini;

    _CoM_lay_lp = -_AA1y*_SS2-_AA2y*_SS3-(_AA3y-2*_comay_min)*_SS4;
    _det_CoM_lay_lp = -(-_AA1y*_SS2-_AA2y*_SS3-(_AA3y-2*_comay_min)*_SS4)*_vari_ini;      
    
    _CoM_lax_upx.row(0) = _CoM_lax_up; _CoM_lax_upx.row(1) = _CoM_lax_lp; 
    _CoM_lax_upx.row(2) = _CoM_lay_up; _CoM_lax_upx.row(3) = _CoM_lay_lp;
    _det_CoM_lax_upx.row(0) = _det_CoM_lax_up;  _det_CoM_lax_upx.row(1) = _det_CoM_lax_lp; 
    _det_CoM_lax_upx.row(2) = _det_CoM_lay_up;  _det_CoM_lax_upx.row(3) = _det_CoM_lay_lp;
    
//     cout <<"_CoM_lax_upx:"<<endl<<_CoM_lax_upx<<endl;
//     cout <<"_det_CoM_lax_upx:"<<endl<<_det_CoM_lax_upx<<endl;       
    
    
    
//   CoM velocity_inremental boundary  
    _VAA= cosh(_Wn*_dt); _VCCx = _comx_feed(0,i-1)-_footx_ref(_periond_i-1,0); _VBBx = _Wn*_VCCx*sinh(_Wn*_dt); 
		         _VCCy = _comy_feed(0,i-1)-_footy_ref(_periond_i-1,0); _VBBy = _Wn*_VCCy*sinh(_Wn*_dt);

    _VAA1x = _VAA*_Wn; _VAA2x = -2* _VAA*_VCCx*_Wn; _VAA3x = 2*_VBBx - 2*_comvx_feed(0,i-1); 
    _VAA1y = _VAA*_Wn; _VAA2y = -2* _VAA*_VCCy*_Wn; _VAA3y = 2*_VBBy - 2*_comvy_feed(0,i-1);

    _CoM_lvx_up = _VAA1x*_SS1+_VAA2x*_SS3+(_VAA3x-2*_comax_max*_dt)*_SS4;
    _det_CoM_lvx_up = -(_VAA1x*_SS1+_VAA2x*_SS3+(_VAA3x-2*_comax_max*_dt)*_SS4)*_vari_ini;

    _CoM_lvx_lp = -_VAA1x*_SS1-_VAA2x*_SS3-(_VAA3x-2*_comax_min*_dt)*_SS4;
    _det_CoM_lvx_lp = -(-_VAA1x*_SS1-_VAA2x*_SS3-(_VAA3x-2*_comax_min*_dt)*_SS4)*_vari_ini; 

    _CoM_lvy_up = _VAA1y*_SS2+_VAA2y*_SS3+(_VAA3y-2*_comay_max*_dt)*_SS4;
    _det_CoM_lvy_up = -(_VAA1y*_SS2+_VAA2y*_SS3+(_VAA3y-2*_comay_max*_dt)*_SS4)*_vari_ini;

    _CoM_lvy_lp = -_VAA1y*_SS2-_VAA2y*_SS3-(_VAA3y-2*_comay_min*_dt)*_SS4;
    _det_CoM_lvy_lp = -(-_VAA1y*_SS2-_VAA2y*_SS3-(_VAA3y-2*_comay_min*_dt)*_SS4)*_vari_ini;      
    
    _CoM_lvx_upx.row(0) = _CoM_lvx_up; _CoM_lvx_upx.row(1) = _CoM_lvx_lp; 
    _CoM_lvx_upx.row(2) = _CoM_lvy_up; _CoM_lvx_upx.row(3) = _CoM_lvy_lp;
    _det_CoM_lvx_upx.row(0) = _det_CoM_lvx_up;  _det_CoM_lvx_upx.row(1) = _det_CoM_lvx_lp; 
    _det_CoM_lvx_upx.row(2) = _det_CoM_lvy_up;  _det_CoM_lvx_upx.row(3) = _det_CoM_lvy_lp;  
    
    
    
//   CoM intial velocity boundary: check   
    _VAA1x1 = _Wn; _VAA2x1 = -2*_VCCx*_Wn; _VAA3x1 = - 2*_comvx_feed(0,i-1); 
    _VAA1y1 = _Wn; _VAA2y1 = -2*_VCCy*_Wn; _VAA3y1 = - 2*_comvy_feed(0,i-1);

    /// modified!!!
    _CoM_lvx_up1 = _VAA1x1*_SS1+_VAA2x1*_SS3+(_VAA3x1-2*_comax_max*_dt)*_SS4;
    _det_CoM_lvx_up1 = -(_VAA1x1*_SS1+_VAA2x1*_SS3+(_VAA3x1-2*_comax_max*_dt/2.0)*_SS4)*_vari_ini;

    _CoM_lvx_lp1 = -_VAA1x1*_SS1-_VAA2x1*_SS3-(_VAA3x1-2*_comax_min*_dt)*_SS4;
    _det_CoM_lvx_lp1 = -(-_VAA1x1*_SS1-_VAA2x1*_SS3-(_VAA3x1-2*_comax_min*_dt/2.0)*_SS4)*_vari_ini; 

    _CoM_lvy_up1 = _VAA1y1*_SS2+_VAA2y1*_SS3+(_VAA3y1-2*_comay_max*_dt)*_SS4;
    _det_CoM_lvy_up1 = -(_VAA1y1*_SS2+_VAA2y1*_SS3+(_VAA3y1-2*_comay_max*_dt/2.0)*_SS4)*_vari_ini;

    _CoM_lvy_lp1 = -_VAA1y1*_SS2-_VAA2y1*_SS3-(_VAA3y1-2*_comay_min*_dt)*_SS4;
    _det_CoM_lvy_lp1 = -(-_VAA1y1*_SS2-_VAA2y1*_SS3-(_VAA3y1-2*_comay_min*_dt/2.0)*_SS4)*_vari_ini;      
    
    _CoM_lvx_upx1.row(0) = _CoM_lvx_up1; _CoM_lvx_upx1.row(1) = _CoM_lvx_lp1; 
    _CoM_lvx_upx1.row(2) = _CoM_lvy_up1; _CoM_lvx_upx1.row(3) = _CoM_lvy_lp1;
    _det_CoM_lvx_upx1.row(0) = _det_CoM_lvx_up1;  _det_CoM_lvx_upx1.row(1) = _det_CoM_lvx_lp1; 
    _det_CoM_lvx_upx1.row(2) = _det_CoM_lvy_up1;  _det_CoM_lvx_upx1.row(3) = _det_CoM_lvy_lp1;      
    
  
  
}




Eigen::MatrixXd  NLPClass::Matrix_ps(Eigen::MatrixXd a, int nh,Eigen::MatrixXd cxps)
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


Eigen::MatrixXd NLPClass::Matrix_pu(Eigen::MatrixXd a, Eigen::MatrixXd b, int nh, Eigen::MatrixXd cxpu)
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






void NLPClass::solve_stepping_timing_twosteps()
{
/*  int nVars = 8;
  int nEqCon = 2;
  int nIneqCon = 24 + 12;
  resizeQP(nVars, nEqCon, nIneqCon);*/	    

  _G = _SQ_goal3;
  _g0 = _Sq_goal3;
  _X = _vari_ini;

// min 0.5 * x G x + g0 x
// _s.t.
// 		CE^T x + ce0 = 0   ///// equality constraints
// 		CI^T x + ci0 >= 0  //// inequality constraints
  _CI.block<8,8>(0,0) = _trx.transpose() * (-1);
  _CI.block<8,8>(0,8) = _h_lx_upx.transpose() * (-1);
  _CI.block<8,8>(0,16) = _h_lvx_upx.transpose() * (-1);
  _CI.block<8,4>(0,24) = _CoM_lax_upx.transpose() * (-1);
  _CI.block<8,4>(0,28) = _CoM_lvx_upx.transpose() * (-1);
  _CI.block<8,4>(0,32) = _CoM_lvx_upx1.transpose() * (-1);
  
  _ci0.block<8,1>(0, 0) = _det_trx;
  _ci0.block<8,1>(8, 0) = _det_h_lx_upx;
  _ci0.block<8,1>(16, 0) = _det_h_lvx_upx;
  _ci0.block<4,1>(24, 0) = _det_CoM_lax_upx;
  _ci0.block<4,1>(28, 0) = _det_CoM_lvx_upx;
  _ci0.block<4,1>(32, 0) = _det_CoM_lvx_upx1;

  
  _CE = _trxx.transpose()*(-1);
  _ce0 = _det_trxx;
  
  
  
  Solve();  

}


void NLPClass::solve_stepping_timing_onestep()
{
/*  int nVars = 8;
  int nEqCon = 2;
  int nIneqCon = 24 + 12;
  resizeQP(nVars, nEqCon, nIneqCon);*/	    

  _G = _SQ_goal3;
  _g0 = _Sq_goal3;
  _X = _vari_ini;

// min 0.5 * x G x + g0 x
// _s.t.
// 		CE^T x + ce0 = 0   ///// equality constraints
// 		CI^T x + ci0 >= 0  //// inequality constraints
  _CI.block<8,8>(0,0) = _trx.transpose() * (-1);
  _CI.block<8,8>(0,8) = _h_lx_upx.transpose() * (-1);
  _CI.block<8,8>(0,16) = _h_lvx_upx.transpose() * (-1);
  _CI.block<8,4>(0,24) = _CoM_lax_upx.transpose() * (-1);
  _CI.block<8,4>(0,28) = _CoM_lvx_upx.transpose() * (-1);
  _CI.block<8,4>(0,32) = _CoM_lvx_upx1.transpose() * (-1);
  
  _ci0.block<8,1>(0, 0) = _det_trx;
  _ci0.block<8,1>(8, 0) = _det_h_lx_upx;
  _ci0.block<8,1>(16, 0) = _det_h_lvx_upx;
  _ci0.block<4,1>(24, 0) = _det_CoM_lax_upx;
  _ci0.block<4,1>(28, 0) = _det_CoM_lvx_upx;
  _ci0.block<4,1>(32, 0) = _det_CoM_lvx_upx1;

  
  _CE.block<8,2>(0, 0) = _trxx.transpose()*(-1);
  _CE.block<8,1>(0, 2) = _Lxx_ref1_e.transpose()*(-1);  
  _CE.block<8,1>(0, 3) = _Lyy_ref1_e.transpose()*(-1);  
  _CE.block<8,1>(0, 4) = _tr1_ref1_e.transpose()*(-1);  
  _CE.block<8,1>(0, 5) = _tr2_ref1_e.transpose()*(-1);    
  
  
  _ce0.block<2,1>(0, 0) = _det_trxx;
  _ce0.block<1,1>(2, 0) = _Lxx_ref1_eq; 
  _ce0.block<1,1>(3, 0) = _Lyy_ref1_eq;
  _ce0.block<1,1>(4, 0) = _tr1_ref1_eq;
  _ce0.block<1,1>(5, 0) = _tr2_ref1_eq;  

 
  Solve();  

}


void NLPClass::solve_stepping_timing_non()
{    

  _G = _SQ_goal3;
  _g0 = _Sq_goal3;
  _X = _vari_ini;

// min 0.5 * x G x + g0 x
// _s.t.
// 		CE^T x + ce0 = 0   ///// equality constraints
// 		CI^T x + ci0 >= 0  //// inequality constraints
  _CI.block<8,8>(0,0) = _trx.transpose() * (-1);
  _CI.block<8,8>(0,8) = _h_lx_upx.transpose() * (-1);
  _CI.block<8,8>(0,16) = _h_lvx_upx.transpose() * (-1);
  _CI.block<8,4>(0,24) = _CoM_lax_upx.transpose() * (-1);
  _CI.block<8,4>(0,28) = _CoM_lvx_upx.transpose() * (-1);
  _CI.block<8,4>(0,32) = _CoM_lvx_upx1.transpose() * (-1);
  
  _ci0.block<8,1>(0, 0) = _det_trx;
  _ci0.block<8,1>(8, 0) = _det_h_lx_upx;
  _ci0.block<8,1>(16, 0) = _det_h_lvx_upx;
  _ci0.block<4,1>(24, 0) = _det_CoM_lax_upx;
  _ci0.block<4,1>(28, 0) = _det_CoM_lvx_upx;
  _ci0.block<4,1>(32, 0) = _det_CoM_lvx_upx1;

  
//   _CE = _trxx.transpose()*(-1);
//   _ce0 = _det_trxx;
 
  _CE.block<8,2>(0, 0) = _trxx.transpose()*(-1);
  _CE.block<8,1>(0, 2) = _Lxx_ref1_e.transpose()*(-1);  
  _CE.block<8,1>(0, 3) = _Lyy_ref1_e.transpose()*(-1);  
  _CE.block<8,1>(0, 4) = _tr1_ref1_e.transpose()*(-1);  
  _CE.block<8,1>(0, 5) = _tr2_ref1_e.transpose()*(-1);    
  _CE.block<8,1>(0, 6) = _Lxx_ref_e.transpose()*(-1);  
  _CE.block<8,1>(0, 7) = _Lyy_ref_e.transpose()*(-1);  
  _CE.block<8,1>(0, 8) = _tr1_ref_e.transpose()*(-1);  
  _CE.block<8,1>(0, 9) = _tr2_ref_e.transpose()*(-1);   
  
  _ce0.block<2,1>(0, 0) = _det_trxx;
  _ce0.block<1,1>(2, 0) = _Lxx_ref1_eq; 
  _ce0.block<1,1>(3, 0) = _Lyy_ref1_eq;
  _ce0.block<1,1>(4, 0) = _tr1_ref1_eq;
  _ce0.block<1,1>(5, 0) = _tr2_ref1_eq;  
  _ce0.block<1,1>(6, 0) = _Lxx_ref_eq; 
  _ce0.block<1,1>(7, 0) = _Lyy_ref_eq;
  _ce0.block<1,1>(8, 0) = _tr1_ref_eq;
  _ce0.block<1,1>(9, 0) = _tr2_ref_eq;  
  
  
  
  Solve();  

}









void NLPClass::Solve()
{
// min 0.5 * x G x + g0 x
// _s.t.
// 		CE^T x + ce0 = 0
// 		CI^T x + ci0 >= 0
		solveQP();

}


//// each time: planer for foot_trajectory:
void NLPClass::Foot_trajectory_solve(int j_index, bool _stopwalking)
{
  
  
    
  _footxyz_real(1,0) = -_stepwidth(0);
 
  
  Eigen::Vector3d t_plan(0,0,0);  
  
//// judge if stop  
  if(_stopwalking)  
  {
    
    for (int i_t = _bjx1+1; i_t < _footstepsnumber; i_t++) {	  
      _lift_height_ref(i_t) = 0;  
    }	  

  }    
  
  
//   foot trajectory generation:
  if (_bjx1 >= 2)
  {
//     cout << "_bjx1 >= 2"<<endl;
    if (_bjx1 % 2 == 0)           //odd:left support
    {
   
      _Lfootx(j_index) = _Lfootx(j_index-1);
      _Lfooty(j_index) = _Lfooty(j_index-1);
      _Lfootz(j_index) = _Lfootz(j_index-1);
      
      _Lfootx(j_index+1) = _Lfootx(j_index-1);
      _Lfooty(j_index+1) = _Lfooty(j_index-1);
      _Lfootz(j_index+1) = _Lfootz(j_index-1); 
      
      if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double support
      {
	_Rfootx(j_index) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfooty(j_index) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfootz(j_index) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1);	
	
	_Rfootx(j_index+1) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfooty(j_index+1) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfootz(j_index+1) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1);	
	
	
      }
      else
      {
	
// 	cout << "ssp"<<endl;
	//initial state and final state and the middle state
	double t_des = (j_index  - round(_tx(_bjx1-1)/_dt) +1)*_dt;
//	Eigen::Vector3d t_plan;
	t_plan(0) = t_des - _dt;
// 	t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
// 	t_plan(2) = _ts(_bjx1-1) + 0.0001;
	t_plan(1) = round((_td(_bjx1-1) + _ts(_bjx1-1))/2/_dt)*_dt + _dt/2 ;
	t_plan(2) = round( _ts(_bjx1-1)/_dt)*_dt+ _dt/2;	
	
	
	if (abs(t_des - _ts(_bjx1-1)) <= (_dt + 0.0005))
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
// 	    if (abs(t_plan(0)-t_plan(1))<_dt/2)
// 	    {
// 	      if (t_plan(0)>t_plan(1))
// 	      {
// 		t_plan(1) = t_plan(1)-_dt/4;
// 	      }
// 	      else
// 	      {
// 		t_plan(1) = t_plan(1)+_dt/4;
// 	      }
// 	      
// 	    }
	  
// 	  if ((abs(t_plan(0))<=0.01)||(abs(t_plan(2))<=0.4)||((abs(t_plan(0)-t_plan(1))<=_dt)))
// 	  {
// // 	    _Rfootx(j_index) = _Rfootx(j_index-1);
// // 	    _Rfooty(j_index) = _Rfooty(j_index-1);
// // 	    _Rfootz(j_index) = _Rfootz(j_index-1);
// // 	    
// // 	    _Rfootx(j_index+1) = _Rfootx(j_index-1);
// // 	    _Rfooty(j_index+1) = _Rfooty(j_index-1);
// // 	    _Rfootz(j_index+1) = _Rfootz(j_index-1); 	    
// // 	    
// 	  }
// 	  else	      
// 	  {
	    Eigen::Matrix<double,7,7> AAA_inv = solve_AAA_inv_x(t_plan,j_index);	  
	    
	    
	    
		    
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
	    if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1)+_dt)
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
	    Rfootz_plan(4) = _footxyz_real(2,_bjxx);  Rfootz_plan(5) = 0;                   Rfootz_plan(6) = 0;	
	    
	    Eigen::Matrix<double, 7, 1> Rfootz_co;
	    Rfootz_co.setZero();
	    Rfootz_co = AAA_inv * Rfootz_plan;
	    
	    _Rfootz(j_index) = t_a_plan * Rfootz_co;
	    _Rfootvz(j_index) = t_a_planv * Rfootz_co;
	    _Rfootaz(j_index) = t_a_plana * Rfootz_co;	
	  
	    
	    _Rfootx(j_index+1) = _Rfootx(j_index)+_dt * _Rfootvx(j_index);
	    _Rfooty(j_index+1) = _Rfooty(j_index)+_dt * _Rfootvy(j_index);
	    _Rfootz(j_index+1) = _Rfootz(j_index)+_dt * _Rfootvz(j_index);	
	    
/*	    cout << "_t_plan:"<<t_plan.transpose()<<endl;
            cout << "_footxyz_real:"<<_footxyz_real.col(_bjxx)<<endl;
	    
	    cout << "_Rfootx(j_index):"<<_Rfootx(j_index)<<endl;  
	    cout << "_Rfooty(j_index):"<<_Rfooty(j_index)<<endl;   
	    cout << "_Rfootz(j_index):"<<_Rfootz(j_index)<<endl; */ 	    
	    
	    
/*	  }*/  
	  
	}
      }   
    }
    
    else                       //right support
    {
//       cout << "right support"<<endl;
/*      _Rfootx(j_index) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfooty(j_index) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfootz(j_index) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1);
      
      _Rfootx(j_index+1) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfooty(j_index+1) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfootz(j_index+1) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1); */  
      _Rfootx(j_index) = _Rfootx(j_index-1);
      _Rfooty(j_index) = _Rfooty(j_index-1);
      _Rfootz(j_index) = _Rfootz(j_index-1);
      
      _Rfootx(j_index+1) = _Rfootx(j_index-1);
      _Rfooty(j_index+1) = _Rfooty(j_index-1);
      _Rfootz(j_index+1) = _Rfootz(j_index-1); 


      
      if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double suppot
      {
// 	cout << "dsp"<<endl;
	_Lfootx(j_index) = _Lfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfooty(j_index) = _Lfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfootz(j_index) = _Lfootz(round(_tx(_bjx1-1)/_dt) -1-1);

	_Lfootx(j_index+1) = _Lfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfooty(j_index+1) = _Lfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfootz(j_index+1) = _Lfootz(round(_tx(_bjx1-1)/_dt) -1-1);
	
      }
      else
      {
// 	cout << "ssp"<<endl;
	//initial state and final state and the middle state
	double t_des = (j_index  - round(_tx(_bjx1-1)/_dt) +1)*_dt;
// 	Eigen::Vector3d t_plan;
	t_plan(0) = t_des - _dt;
	t_plan(0) = t_des - _dt;
// 	t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
// 	t_plan(2) = _ts(_bjx1-1) + 0.0001;
	t_plan(1) = round((_td(_bjx1-1) + _ts(_bjx1-1))/2/_dt)*_dt + _dt/2 ;
	t_plan(2) = round( _ts(_bjx1-1)/_dt)*_dt+ _dt/2;
	
	if (abs(t_des - _ts(_bjx1-1)) <= (_dt + 0.0005))
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
	  
// 	  if ((abs(t_plan(0))<=0.01)||(abs(t_plan(2))<=0.4)||((abs(t_plan(0)-t_plan(1))<=_dt)))
// 	  {
// /*	    _Lfootx(j_index) = _Lfootx(j_index-1);
// 	    _Lfooty(j_index) = _Lfooty(j_index-1);
// 	    _Lfootz(j_index) = _Lfootz(j_index-1);
// 	    
// 	    _Lfootx(j_index+1) = _Lfootx(j_index-1);
// 	    _Lfooty(j_index+1) = _Lfooty(j_index-1);
// 	    _Lfootz(j_index+1) = _Lfootz(j_index-1); */	    
// 	    
// 	  }
// 	  else
// 	  {
// 	    if (abs(t_plan(0)-t_plan(1))<_dt/2)
// 	    {
// 	      if (t_plan(0)>t_plan(1))
// 	      {
// 		t_plan(1) = t_plan(1)-_dt/4;
// 	      }
// 	      else
// 	      {
// 		t_plan(1) = t_plan(1)+_dt/4;
// 	      }
// 	      
// 	    }
	     
	  
	    Eigen::Matrix<double,7,7> AAA_inv = solve_AAA_inv_x(t_plan,j_index);          
		    
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
	    
  // 	  cout <<"AAA="<<endl<<AAA<<endl;
  // 	  cout <<"AAA_inverse="<<endl<<AAA.inverse()<<endl;	  
  // 	  cout <<"t_des="<<endl<<t_des<<endl;
  // 	  cout <<"t_plan="<<endl<<t_plan<<endl;
	    
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
	    if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1)+_dt)
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
	    Lfootz_plan(4) = _footxyz_real(2,_bjxx);  Lfootz_plan(5) = 0;                   Lfootz_plan(6) = 0;		  
	    
	    
	    Eigen::VectorXd Lfootz_co;
	    Lfootz_co.setZero(7);
	    Lfootz_co = AAA_inv * Lfootz_plan;
	    
	    _Lfootz(j_index) = t_a_plan * Lfootz_co;
	    _Lfootvz(j_index) = t_a_planv * Lfootz_co;
	    _Lfootaz(j_index) = t_a_plana * Lfootz_co;
	    
	    
	    _Lfootx(j_index+1) = _Lfootx(j_index)+_dt * _Lfootvx(j_index);
	    _Lfooty(j_index+1) = _Lfooty(j_index)+_dt * _Lfootvy(j_index);
	    _Lfootz(j_index+1) = _Lfootz(j_index)+_dt * _Lfootvz(j_index);	    
// 	  }
	  
// 	    cout << "_t_plan:"<<t_plan.transpose()<<endl;
// 	    cout << "_footxyz_real:"<<_footxyz_real.col(_bjxx)<<endl;	    
// 	    cout << "_Lfootx(j_index):"<<_Lfootx(j_index)<<endl;  
// 	    cout << "_Lfooty(j_index):"<<_Lfooty(j_index)<<endl;   
// 	    cout << "_Lfootz(j_index):"<<_Lfootz(j_index)<<endl;  
	  

	  
	  
	}
      }

    }
      
  }
  else
  {
    _Rfooty(j_index) = -_stepwidth(0);
    _Lfooty(j_index) = _stepwidth(0);
  }
    

   
//    cout << "_t_plan:"<<t_plan.transpose()<<endl;
   
   if ((abs(_Lfootx(j_index))>1)||(abs(_Lfooty(j_index))>1)||(abs(_Lfootz(j_index))>1)||(abs(_Rfootx(j_index))>1)||(abs(_Rfooty(j_index))>1)||(abs(_Rfootz(j_index))>1))
   {
/*   cout << "j_index:"<<j_index<<endl;*/ 
//    cout << "_Lfootx(j_index):"<<_Lfootx(j_index)<<endl;  
//    cout << "_Lfooty(j_index):"<<_Lfooty(j_index)<<endl;   
//    cout << "_Lfootz(j_index):"<<_Lfootz(j_index)<<endl;   
//    cout << "_Rfootx(j_index):"<<_Rfootx(j_index)<<endl;  
//    cout << "_Rfooty(j_index):"<<_Rfooty(j_index)<<endl;   
//    cout << "_Rfootz(j_index):"<<_Rfootz(j_index)<<endl;    
//     cout << "_t_plan:"<<t_plan<<endl;
//     cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!1:"<<j_index<<endl;
   
   }  
  
}


//// each time: planer for CoM_height_trajectory: only the CoM_height_generation
void NLPClass::CoM_height_solve(int j_index, bool _stopwalking)
{     
  
//   CoM trajectory generation:
  if (_bjx1 >= 2)
  {      
    if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt >= _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double suppot
    {
      _comz(j_index) = _footxyz_real(2,_bjx1-1)+_hcom;
      _comz(j_index+1) = _footxyz_real(2,_bjx1-1)+_hcom;
           
    }
    else
    {
      //initial state and final state and the middle state
      double t_des = (j_index +1 - round(_tx(_bjx1-1)/_dt) )*_dt;
      Eigen::Vector3d t_plan;
//       t_plan(0) = t_des - _dt;
      t_plan(0) = 0.0001;
      t_plan(1) = _td(_bjx1-1)/2 + 0.0001;
      t_plan(2) = _td(_bjx1-1) + 0.0001;
      

      
      Eigen::Matrix<double,7,7> AAA_inv = solve_AAA_inv_x(t_plan, j_index);         
	      
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

      Eigen::Matrix<double, 7, 1> comz_plan;
      comz_plan.setZero();			
//       comz_plan(0) = _comvz(j_index-1);               comz_plan(1) = _comaz(j_index-1); comz_plan(2) = _comz(j_index-1); comz_plan(3) = (_comz(j_index)+_footxyz_real(2,_bjx1-1)+_hcom)/2;
      comz_plan(0) = 0;                 comz_plan(1) = 0;                 
      comz_plan(2) = _footxyz_real(2,_bjx1-2)+_hcom; 
      comz_plan(3) = (_footxyz_real(2,_bjx1-2)+_footxyz_real(2,_bjx1-1))/2+_hcom;      
      comz_plan(4) = _footxyz_real(2,_bjx1-1)+_hcom;  
      comz_plan(5) = 0;                 comz_plan(6) = 0;		  
      
      
      Eigen::VectorXd Lfootz_co;
      Lfootz_co.setZero(7);
      Lfootz_co = AAA_inv * comz_plan;
      
      _comz(j_index) = t_a_plan * Lfootz_co;
      _comvz(j_index) = t_a_planv * Lfootz_co;
      _comaz(j_index) = t_a_plana * Lfootz_co;
      
      
      _comx(j_index+1) = _comx(j_index)+_dt * _comvx(j_index);
      _comy(j_index+1) = _comy(j_index)+_dt * _comvy(j_index);
      _comz(j_index+1) = _comz(j_index)+_dt * _comvz(j_index);      
 
      
    }

      
  }
  else
  {
    _comz(j_index) =_hcom;
  }
    

  
}

////
Eigen::Vector3d NLPClass::X_CoM_position_squat(int walktime, double dt_sample)
{
  Eigen::Vector3d com_inte;
  com_inte.setZero();


  double t_des;
  t_des = walktime * dt_sample;

  Eigen::Vector3d t_plan;
  t_plan(0) = 0.00001;
  t_plan(1) = _height_squat_time/2+0.0001;
  t_plan(2) = _height_squat_time+0.0001;


  if (t_des<=_height_squat_time)
  {
    int j_x =1;
    
    Eigen::Matrix<double,7,7> AAA_inv = solve_AAA_inv_x(t_plan,j_x);
	    
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










int NLPClass::Get_maximal_number_reference()
{
  int nsum_max;
  nsum_max = (_nsum -_n_loop_omit-1);  
  return nsum_max;
}

int NLPClass::Get_maximal_number(double dtx)
{
  int nsum_max;
  nsum_max = (_nsum -_n_loop_omit-1)*floor(_dt/dtx);
  
  return nsum_max;
}



////====================================================================================================================
/////////////////////////// using the lower-level control-loop  sampling time as the reference: every 5ms;  at the same time: just using the next one position + next one velocity

Eigen::Vector3d NLPClass::XGetSolution_CoM_position(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
  //reference com position
	_comz(0) = _HCOM-_height_offset;
	_comz(1) = _HCOM-_height_offset;
	_comz(2) = _HCOM-_height_offset;
	_comz(3) = _HCOM-_height_offset;
	_comz(4) = _HCOM-_height_offset;	
	
	Eigen::Vector3d com_inte;	
	
	if (walktime>=2)
	{
	  int t_int = floor(walktime* dt_sample/ _dt);	  
	  ///// chage to be relative time
	  double t_cur = walktime * dt_sample ;
	  

	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = t_cur - 2*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(1) = t_cur - 1*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);
	  
//	  cout<< "t_plan"<<t_plan<<endl;
          
	  solve_AAA_inv(t_plan);
	   
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_cur-( t_cur - 2*dt_sample), 3);   
	  t_a_plan(1) = pow(t_cur-( t_cur - 2*dt_sample), 2);   
	  t_a_plan(2) = pow(t_cur-( t_cur - 2*dt_sample), 1);  
	  t_a_plan(3) = pow(t_cur-( t_cur - 2*dt_sample), 0); 


	  Eigen::Matrix<double, 1, 4> t_a_planv;
	  t_a_planv.setZero();
	  t_a_planv(0) = 3*pow(2*dt_sample, 2);   t_a_planv(1) = 2*pow(2*dt_sample, 1);   
	  t_a_planv(2) = 1;   t_a_planv(3) = 0; 	  
	  
	  // COM&&foot trajectory interpolation	  	  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = body_in1(0); temp(1) = body_in2(0); temp(2) =  _comx(t_int); temp(3) = _comvx(t_int);	  
	  com_inte(0) = t_a_plan * (_AAA_inv)*temp;
	  _comxyzx(0) = com_inte(0);
	  _comvxyzx(0) = t_a_planv * (_AAA_inv)*temp;	  
	  
	  
	  temp(0) = body_in1(1); temp(1) = body_in2(1); temp(2) =  _comy(t_int); temp(3) = _comvy(t_int);	  
	  com_inte(1) = t_a_plan * (_AAA_inv)*temp;
	  _comxyzx(1) = com_inte(1);
	  _comvxyzx(1) = t_a_planv * (_AAA_inv)*temp;	  
	  
	  
	  temp(0) = body_in1(2); temp(1) = body_in2(2); temp(2) =  _comz(t_int); temp(3) = _comvz(t_int);	  
	  com_inte(2) = t_a_plan *(_AAA_inv)*temp;
	  _comxyzx(2) = com_inte(2);
	  _comvxyzx(2) = t_a_planv * (_AAA_inv)*temp;
	  
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
	//cout<<"com_height"<< com_inte(2)<<endl;
	
}

Eigen::Vector3d NLPClass::XGetSolution_body_inclination(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
		
	Eigen::Vector3d com_inte(0,0,0);		
	if (walktime>=2)
	{
	  int t_int = floor(walktime / (_dt / dt_sample) );
	 
	  
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



Eigen::Vector3d NLPClass::XGetSolution_Foot_positionR(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{	
	int t_int= floor(walktime* dt_sample/ _dt  );
	
	///////////// 4th order interpolation
	Eigen::Vector3d com_inte;	
	
	if (t_int>=1)
	{
	 
	  
	  double t_cur = walktime * dt_sample ;
	  

	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = t_cur - 2*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(1) = t_cur - 1*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);
	  
//	  cout<< "t_plan"<<t_plan<<endl;
          
	  solve_AAA_inv(t_plan);	  
	  
	  	  
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_cur-( t_cur - 2*dt_sample), 3);  
	  t_a_plan(1) = pow(t_cur-( t_cur - 2*dt_sample), 2);  
	  t_a_plan(2) = pow(t_cur-( t_cur - 2*dt_sample), 1);  
	  t_a_plan(3) = pow(t_cur-( t_cur - 2*dt_sample), 0); 


	  
	  
	  // COM&&foot trajectory interpolation  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = body_in1(0); temp(1) = body_in2(0); temp(2) = _Rfootx(t_int); temp(3) = _Rfootvx(t_int);	  
	  com_inte(0) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(1); temp(1) = body_in2(1); temp(2) = _Rfooty(t_int); temp(3) = _Rfootvy(t_int);

	  com_inte(1) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(2); temp(1) = body_in2(2); temp(2) = _Rfootz(t_int); temp(3) = _Rfootvz(t_int);	  
	  com_inte(2) = t_a_plan *(_AAA_inv)*temp;
	  
	  
	  
	  /////// linear intepolation
	  com_inte(0) = (_Rfootx(t_int) - _Rfootx(t_int-1))*((t_cur -t_int*_dt) / _dt) +_Rfootx(t_int-1);
	  com_inte(1) = (_Rfooty(t_int) - _Rfooty(t_int-1))*((t_cur -t_int*_dt) / _dt) +_Rfooty(t_int-1);
	  com_inte(2) = (_Rfootz(t_int) - _Rfootz(t_int-1))*((t_cur -t_int*_dt) / _dt) +_Rfootz(t_int-1);
	  
	  
	  
	  
	  
	}
	else
	{
/*	  com_inte(0) = body_in3(0);	  
	  com_inte(1) = body_in3(1);	  	  
	  com_inte(2) = body_in3(2);*/ 
	  com_inte(0) = 0;	  
	  com_inte(1) = -_half_hip_width;	  	  
	  com_inte(2) = 0;  	  
	}

	 
	_Rfootxyzx(0) = com_inte(0);  
	_Rfootxyzx(1) = com_inte(1);   
	_Rfootxyzx(2) = com_inte(2);   	

 	return com_inte;
	
	
	
	
}

Eigen::Vector3d NLPClass::XGetSolution_Foot_positionL(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
	
	Eigen::Vector3d com_inte;

        int t_int = floor(walktime* dt_sample/ _dt  );	
	
	if (t_int>=1)
	{
	  

	  double t_cur = walktime * dt_sample ;
	  

	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = t_cur - 2*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(1) = t_cur - 1*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);
	  
//	  cout<< "t_plan"<<t_plan<<endl;
          
	  solve_AAA_inv(t_plan);
	  

	  
	  Eigen::RowVector4d t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_cur-( t_cur - 2*dt_sample), 3);  
	  t_a_plan(1) = pow(t_cur-( t_cur - 2*dt_sample), 2);   
	  t_a_plan(2) = pow(t_cur-( t_cur - 2*dt_sample), 1); 
	  t_a_plan(3) = pow(t_cur-( t_cur - 2*dt_sample), 0); 


	  
	  // COM&&foot trajectory interpolation

	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = body_in1(0); temp(1) = body_in2(0); temp(2) = _Lfootx(t_int); temp(3) = _Lfootvx(t_int);	  
	  com_inte(0) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(1); temp(1) = body_in2(1); temp(2) = _Lfooty(t_int); temp(3) = _Lfootvy(t_int);

	  com_inte(1) = t_a_plan * (_AAA_inv)*temp;
	  temp(0) = body_in1(2); temp(1) = body_in2(2); temp(2) = _Lfootz(t_int); temp(3) = _Lfootvz(t_int);	  
	  com_inte(2) = t_a_plan *(_AAA_inv)*temp;
	  
	  /////// linear intepolation
	  com_inte(0) = (_Lfootx(t_int) - _Lfootx(t_int-1))*((t_cur -t_int*_dt) / _dt) +_Lfootx(t_int-1);
	  com_inte(1) = (_Lfooty(t_int) - _Lfooty(t_int-1))*((t_cur -t_int*_dt) / _dt) +_Lfooty(t_int-1);
	  com_inte(2) = (_Lfootz(t_int) - _Lfootz(t_int-1))*((t_cur -t_int*_dt) / _dt) +_Lfootz(t_int-1);	  
	  
	  
	}
	else
	{
// 	  com_inte(0) = body_in3(0);	  
// 	  com_inte(1) = body_in3(1);	  	  
// 	  com_inte(2) = body_in3(2);
	  com_inte(0) = 0;	  
	  com_inte(1) = _half_hip_width;	  	  
	  com_inte(2) = 0;  	  
	}
	
	_Lfootxyzx(0) = com_inte(0);  
	_Lfootxyzx(1) = com_inte(1);   
	_Lfootxyzx(2) = com_inte(2);   
	

 	return com_inte;
	
}



//////////////////////////////////////////////////////
// Eigen::Matrix<double,1,6> NLPClass::XGetSolution_Foot_position_KMP(int walktime, double dt_sample,int j_index, bool _stopwalking)
// {
//   
//   ///////walktime=====>ij;   int j_index====>i;  dt_sample========>dtx;   
//    Eigen::Matrix<double,1,6> com_inte;	
//  
// //// judge if stop  
//   if(_stopwalking)  
//   {
//     
//     for (int i_t = _bjx1+1; i_t < _footstepsnumber; i_t++) {	
//       if (i_t == _bjx1)
//       {
//       _lift_height_ref(i_t) = _lift_height_ref(i_t)/2;	
//       }
//       else
//       {
//       _lift_height_ref(i_t) = 0.01;  	
//       }
// 
//     }	  
// 
//   }  
//  
//  
//  
// //   
//   //// three via_points: time, mean, sigma 
//   vec via_point1 =  zeros<vec>(43);
//   vec via_point2 =  zeros<vec>(43);
//   vec via_point3 =  zeros<vec>(43);
//   
//   via_point1(7) =0.00000000001; via_point1(14)=0.00000000001; via_point1(21)=0.00000000001;	
//   via_point1(28)=0.00000000001; via_point1(35)=0.00000000001; via_point1(42)=0.00000000001;  
//   
//   via_point2(0) =0.65/2;
//   via_point2(7) =0.00000000001; via_point2(14)=0.00000000001; via_point2(21)=0.00000000001;	
//   via_point2(28)=0.00000000001; via_point2(35)=0.00000000001; via_point2(42)=0.00000000001;
//   
//   via_point3(0) =0.65;
//   via_point3(7) =0.00000000001; via_point3(14)=0.00000000001; via_point3(21)=0.00000000001;	
//   via_point3(28)=0.00000000001; via_point3(35)=0.00000000001; via_point3(42)=0.00000000001;      
//   
//   
//   
//   double t_des;      /////////desired time during the current step
//   t_des = (walktime+4)*dt_sample - (_tx(_bjx1-1)+_td(_bjx1-1));
//   
// 
//   if (_bjx1 >= 2)
//   {      
//     
//     if (_bjx1 % 2 == 0)           //odd:left support
//     {
//       _Lfootx_kmp(0) = _footxyz_real(0,_bjx1-1);
//       _Lfooty_kmp(0) = _footxyz_real(1,_bjx1-1);
//       _Lfootz_kmp(0) = _footxyz_real(2,_bjx1-1);
//       
//       _Lfootvx_kmp(0) = 0;
//       _Lfootvy_kmp(0) = 0;
//       _Lfootvz_kmp(0) = 0;    
//       
//       if (t_des<=0)  // j_index and _bjx1 coincident with matlab: double support
//       {
// 
// 	_Rfootx_kmp(0) = _footxyz_real(0,_bjx1-2);
// 	_Rfooty_kmp(0) = _footxyz_real(1,_bjx1-2);
// 	_Rfootz_kmp(0) = _footxyz_real(2,_bjx1-2);
// 	
// 	_Rfootvx_kmp(0) = 0;
// 	_Rfootvy_kmp(0) = 0;
// 	_Rfootvz_kmp(0) = 0;  
// 	
//       }
//       else
//       {
// 	
// 	//initial state and final state and the middle state
// 	double t_des_k;
// 	t_des_k = (0.65)/((_ts(_bjx1-1)-_td(_bjx1-1)))*t_des;
// 	
// 	/////////////first sampling time of the current walking cycle: initialize the KMP_data
// 	if (t_des<=dt_sample)
// 	{
// 	  kmp_leg_R.kmp_initialize(_data_kmp, _inDim_kmp, _outDim_kmp, _pvFlag_kmp,_lamda_kmp, _kh_kmp);
// 	//// %%% initial state and final state and the middle state
// 	/////%%%%% be careful the x position start from zero, the y
// 	//// %%%%% position start from -RobotParaClass::HALF_HIP_WIDTH()
// 	  ////// add point************ current status****************////////
// 	  via_point1(0) = (0.65)/((_ts(_bjx1-1)-_td(_bjx1-1)))*(t_des-dt_sample);
// 	  via_point1(1) = _Rfootx_kmp_old(0)-_footxyz_real(0,_bjx1-2);
// 	  via_point1(2) = _Rfooty_kmp_old(0)-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	  via_point1(3) = _Rfootz_kmp_old(0)-_footxyz_real(2,_bjx1-2);
// 	  via_point1(4) = _Rfootvx_kmp_old(0);
// 	  via_point1(5) = _Rfootvy_kmp_old(0);
// 	  via_point1(6) = _Rfootvz_kmp_old(0);
// 	  kmp_leg_R.kmp_insertPoint(via_point1);  // insert point into kmp
// 	  
// 	  ////// add point************ middle point***********////////
// // 	    via_point2(1) = _footxyz_real(0,_bjx1-1)-_footxyz_real(0,_bjx1-2);
// 	  via_point2(1) = (_footxyz_real(0,_bjx1)-_footxyz_real(0,_bjx1-2))/2;
// 	  via_point2(2) = (_footxyz_real(1,_bjx1-2)+_footxyz_real(1,_bjx1))/2-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	  via_point2(3) = (_footxyz_real(2,_bjx1-1)+_lift_height_ref(_bjx1-1))-_footxyz_real(2,_bjx1-2);
// 	  via_point2(4) = (_footxyz_real(0,_bjx1)-_footxyz_real(0,_bjx1-2))/0.65*1.1;
// /*	  if (_bjx1==2)
// 	  {
// 	    via_point2(4) = (_footxyz_real(0,_bjx1)-_footxyz_real(0,_bjx1-2))/0.65*0.8;
// 	  }*/	    
// 	  via_point2(5) = (_footxyz_real(1,_bjx1)-_footxyz_real(1,_bjx1-2))/0.65;
//  	  via_point2(6) = 0;
// /*	  if (_footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)<=0)  ////downstairs
// 	  {
// 	    via_point2(6) = -0.1;
// 	  }
// 	  else
// 	  {
// 	    via_point2(6) = 0;
// 	  }*/	    	    
// 	  kmp_leg_R.kmp_insertPoint(via_point2);  // insert point into kmp
// 
// 	  ////// add point************ final status************////////	  
// 	  via_point3(1) = _footxyz_real(0,_bjx1)-_footxyz_real(0,_bjx1-2)+0.00;
// 	  via_point3(2) = _footxyz_real(1,_bjx1)-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// // 	  via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)-0.0030;
// 	  if (_bjx1<=4)
// 	  {
// 	    via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)-0.002;
// 	  }
// 	  else
// 	  {
// 	    if (_bjx1<=15)
// 	    {
// 	      via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)-0.003;
// 	    }
// 	    else
// 	    {
// // 	      ///obstacle avoidance
// // 	      via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)-0.004;
// // 	      ///obstacle avoidance_m
// 	      via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2);
// 	    }
// 	  }
// 	  via_point3(4) = 0;
// 	  via_point3(5) = 0;
//  	  via_point3(6) = 0.025;
// // 	  if (via_point3(3)<0)  ////downstairs
// // 	  {
// // 	    via_point3(6) = 0.15;
// // 	  }
// // 	  else
// // 	  {
// // 	    via_point3(6) = 0.025;
// // 	  }
// 	  
// 	  
// 	  kmp_leg_R.kmp_insertPoint(via_point3);  // insert point into kmp	
// 	  
// 	  kmp_leg_R.kmp_estimateMatrix();
// 	}
// 	else
// 	{
// 	  if ((abs(_Lxx_ref_real(j_index)-_Lxx_ref_real(j_index-1))>=0.005)||(abs(_Lyy_ref_real(j_index)-_Lyy_ref_real(j_index-1))>=0.005)||(abs(_Ts_ref_real(j_index)-_Ts_ref_real(j_index-1))>=0.005))
// 	  {
// 	    ////// add point************ current status****************////////
// 	    via_point1(0) = (0.65)/((_ts(_bjx1-1)-_td(_bjx1-1)))*(t_des-dt_sample);
// 	    via_point1(1) = _Rfootx_kmp_old(0)-_footxyz_real(0,_bjx1-2);
// 	    via_point1(2) = _Rfooty_kmp_old(0)-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	    via_point1(3) = _Rfootz_kmp_old(0)-_footxyz_real(2,_bjx1-2);
// 	    via_point1(4) = _Rfootvx_kmp_old(0);
// 	    via_point1(5) = _Rfootvy_kmp_old(0);
// 	    via_point1(6) = _Rfootvz_kmp_old(0);
// 	    kmp_leg_R.kmp_insertPoint(via_point1);  // insert point into kmp	
// 	    
// 	    kmp_leg_R.kmp_estimateMatrix();
// 	    
// // 	    cout<<"Lx_real_error"<<endl<<_Lxx_ref_real(j_index)-_Lxx_ref_real(j_index-1)<<endl;
// // 	    cout<<"Ly_real_error"<<endl<<_Lyy_ref_real(j_index)-_Lyy_ref_real(j_index-1)<<endl;
// // 	    cout<<"Ts_real_error"<<endl<<_Ts_ref_real(j_index)-_Ts_ref_real(j_index-1)<<endl;
// 	    
// 	  }
// 	}
// 	
// 	_query_kmp(0) = t_des_k;
// 	kmp_leg_R.kmp_prediction(_query_kmp,_mean_kmp);   ////////time& predictive values
// 	if (t_des_k<0.01){
// /*	cout<<"_query_kmp:"<<endl<<trans(_query_kmp)<<endl;
// 	cout<<"kmp:"<<endl<<trans(_mean_kmp)<<endl;
// 	cout<<"error:"<<trans(_mean_kmp)-trans(via_point1(span(1,6)))<<endl<<endl;*/	      	      
// 	}	  
//        
//       _Rfootx_kmp(0) = _mean_kmp(0)+_footxyz_real(0,_bjx1-2);
//       _Rfooty_kmp(0) = _mean_kmp(1)+(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
//       _Rfootz_kmp(0) = _mean_kmp(2)+_footxyz_real(2,_bjx1-2);
//       
//       _Rfootvx_kmp(0) = _mean_kmp(3);
//       _Rfootvy_kmp(0) = _mean_kmp(4);
//       _Rfootvz_kmp(0) = _mean_kmp(5);
//       
//       }   
//     }
//     
//     else                       //right support
//     {
//       _Rfootx_kmp(0) = _footxyz_real(0,_bjx1-1);
//       _Rfooty_kmp(0) = _footxyz_real(1,_bjx1-1);
//       _Rfootz_kmp(0) = _footxyz_real(2,_bjx1-1);
//       
//       _Rfootvx_kmp(0) = 0;
//       _Rfootvy_kmp(0) = 0;
//       _Rfootvz_kmp(0) = 0;    
//       
//       if (t_des<=0)  // j_index and _bjx1 coincident with matlab: double support
//       {
// 
// 	_Lfootx_kmp(0) = _footxyz_real(0,_bjx1-2);
// 	_Lfooty_kmp(0) = _footxyz_real(1,_bjx1-2);
// 	_Lfootz_kmp(0) = _footxyz_real(2,_bjx1-2);
// 	
// 	_Lfootvx_kmp(0) = 0;
// 	_Lfootvy_kmp(0) = 0;
// 	_Lfootvz_kmp(0) = 0;  
// 	
//       }
//       else
//       {
// 	
// 	//initial state and final state and the middle state
// 	double t_des_k;
// 	t_des_k = (0.65)/((_ts(_bjx1-1)-_td(_bjx1-1)))*t_des;
// 	
// 	if (t_des<=dt_sample)
// 	{	
// 	  kmp_leg_L.kmp_initialize(_data_kmp, _inDim_kmp, _outDim_kmp, _pvFlag_kmp,_lamda_kmp, _kh_kmp);
//          //// %%% initial state and final state and the middle state
//          /////%%%%% be careful the x position start from zero, the y
//          //// %%%%% position start from -RobotParaClass::HALF_HIP_WIDTH()
// 	  ////// add point************ current status****************////////
// 	  via_point1(0) = (0.65)/((_ts(_bjx1-1)-_td(_bjx1-1)))*(t_des-dt_sample);
// 	  via_point1(1) = _Lfootx_kmp_old(0)-_footxyz_real(0,_bjx1-2);
// 	  via_point1(2) = _Lfooty_kmp_old(0)-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	  via_point1(3) = _Lfootz_kmp_old(0)-_footxyz_real(2,_bjx1-2);
// 	  via_point1(4) = _Lfootvx_kmp_old(0);
// 	  via_point1(5) = _Lfootvy_kmp_old(0);
// 	  via_point1(6) = _Lfootvz_kmp_old(0);
//  	  kmp_leg_L.kmp_insertPoint(via_point1);  // insert point into kmp
// 	  
// 	  ////// add point************ middle point***********////////
// // 	  via_point2(1) = _footxyz_real(0,_bjx1-1)-_footxyz_real(0,_bjx1-2);
// 	  via_point2(1) = (_footxyz_real(0,_bjx1)-_footxyz_real(0,_bjx1-2))/2;
// 	  via_point2(2) = (_footxyz_real(1,_bjx1-2)+_footxyz_real(1,_bjx1))/2-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	  via_point2(3) = (_footxyz_real(2,_bjx1-1)+_lift_height_ref(_bjx1-1))-_footxyz_real(2,_bjx1-2);
// 	  via_point2(4) = (_footxyz_real(0,_bjx1)-_footxyz_real(0,_bjx1-2))/0.65*1.1;
// 	  via_point2(5) = (_footxyz_real(1,_bjx1)-_footxyz_real(1,_bjx1-2))/0.65;
// 	  via_point2(6) = 0;
// /*	  if (_footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)<=0)  ////downstairs
// 	  {
// 	    via_point2(6) = -0.1;	      
// 	  }
// 	  else
// 	  {
// 	    via_point2(6) = 0;
// 	  }*/	  	  	  
//  	  kmp_leg_L.kmp_insertPoint(via_point2);  // insert point into kmp
// 
//           ////// add point************ final status************////////	  
// 	  via_point3(1) = _footxyz_real(0,_bjx1)-_footxyz_real(0,_bjx1-2)+0.00;
// 	  via_point3(2) = _footxyz_real(1,_bjx1)-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	  if (_bjx1<=4)
// 	  {
// 	    via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)-0.002;
// 	  }
// 	  else
// 	  {
// 	    if (_bjx1<=15)
// 	    {
// 	      via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)-0.003;
// 	    }
// 	    else
// 	    {
// // 	      ///obstacle avoidance
// // 	      via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2)-0.004;
// // 	      ///obstacle avoidance_m
// 	      via_point3(3) = _footxyz_real(2,_bjx1)-_footxyz_real(2,_bjx1-2);
// 
// 	    }
// 	  }
// 	  via_point3(4) = 0;
// 	  via_point3(5) = 0;
//  	  via_point3(6) = 0.025;
// /*	  if (via_point3(3)<0)  ////downstairs
// 	  {
// 	     via_point3(6) = 0.15;
// 	  }
// 	  else
// 	  {
// 	    via_point3(6) = 0.025;
// 	  }*/	     	 
// 	  kmp_leg_L.kmp_insertPoint(via_point3);  // insert point into kmp
// 	  
// 	  
// 	  kmp_leg_L.kmp_estimateMatrix();
// 	   
// 	  
// 	}
// 	else
// 	{
// 	  if ((abs(_Lxx_ref_real(j_index)-_Lxx_ref_real(j_index-1))>=0.005)||(abs(_Lyy_ref_real(j_index)-_Lyy_ref_real(j_index-1))>=0.005)||(abs(_Ts_ref_real(j_index)-_Ts_ref_real(j_index-1))>=0.005))
// 	  {
// 	    ////// add point************ current status****************////////
// 	    via_point1(0) = (0.65)/((_ts(_bjx1-1)-_td(_bjx1-1)))*(t_des-dt_sample);
// 	    via_point1(1) = _Lfootx_kmp_old(0)-_footxyz_real(0,_bjx1-2);
// 	    via_point1(2) = _Lfooty_kmp_old(0)-(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	    via_point1(3) = _Lfootz_kmp_old(0)-_footxyz_real(2,_bjx1-2);
// 	    via_point1(4) = _Lfootvx_kmp_old(0);
// 	    via_point1(5) = _Lfootvy_kmp_old(0);
// 	    via_point1(6) = _Lfootvz_kmp_old(0);
// 	    kmp_leg_L.kmp_insertPoint(via_point1);  // insert point into kmp
// 	    
// 	    kmp_leg_L.kmp_estimateMatrix();
// 	    
// // 	    cout<<"Lx_real_error"<<endl<<_Lxx_ref_real(j_index)-_Lxx_ref_real(j_index-1)<<endl;
// // 	    cout<<"Ly_real_error"<<endl<<_Lyy_ref_real(j_index)-_Lyy_ref_real(j_index-1)<<endl;
// // 	    cout<<"Ts_real_error"<<endl<<_Ts_ref_real(j_index)-_Ts_ref_real(j_index-1)<<endl;
// 	    
// 	  }
// 	}	
// 	
// 	_query_kmp(0) = t_des_k;
// 	kmp_leg_L.kmp_prediction(_query_kmp,_mean_kmp);   ////////time& predictive values
// 	
// /*	if (t_des_k<0.01){
// 	cout<<"_query_kmp:"<<endl<<trans(_query_kmp)<<endl;
// 	cout<<"kmp:"<<endl<<trans(_mean_kmp)<<endl;
// 	cout<<"error:"<<endl<<trans(_mean_kmp)-trans(via_point1(span(1,6)))<<endl;	
// 	}*/		
// 	_Lfootx_kmp(0) = _mean_kmp(0)+_footxyz_real(0,_bjx1-2);
// 	_Lfooty_kmp(0) = _mean_kmp(1)+(_footxyz_real(1,_bjx1-2)-(-_half_hip_width));
// 	_Lfootz_kmp(0) = _mean_kmp(2)+_footxyz_real(2,_bjx1-2);
// 	
// 	_Lfootvx_kmp(0) = _mean_kmp(3);
// 	_Lfootvy_kmp(0) = _mean_kmp(4);
// 	_Lfootvz_kmp(0) = _mean_kmp(5); 	  	
// 	
//       }   
// 
//     }
//       
//     
//   }
//   
//   
//   /////Rfoot,xyz, Lfoot,XYZ
//   com_inte(0) = _Rfootx_kmp(0);
//   com_inte(1) = _Rfooty_kmp(0);
//   com_inte(2) = _Rfootz_kmp(0);
// 
//   com_inte(3) = _Lfootx_kmp(0);
//   com_inte(4) = _Lfooty_kmp(0);
//   com_inte(5) = _Lfootz_kmp(0);
//   
//   _Rfootx_kmp_old = _Rfootx_kmp;
//   _Rfooty_kmp_old = _Rfooty_kmp;
//   _Rfootz_kmp_old = _Rfootz_kmp;
//   _Rfootvx_kmp_old = _Rfootvx_kmp;
//   _Rfootvy_kmp_old = _Rfootvy_kmp;
//   _Rfootvz_kmp_old = _Rfootvz_kmp;
//   
//   _Rfootx_kmp_old = _Lfootx_kmp;
//   _Lfooty_kmp_old = _Lfooty_kmp;
//   _Lfootz_kmp_old = _Lfootz_kmp;
//   _Lfootvx_kmp_old = _Lfootvx_kmp;
//   _Lfootvy_kmp_old = _Lfootvy_kmp;
//   _Lfootvz_kmp_old = _Lfootvz_kmp;
//   
//   
//   
//   _Lfootxyzx(0) = com_inte(3);  
//   _Lfootxyzx(1) = com_inte(4);   
//   _Lfootxyzx(2) = com_inte(5);   
//   _Rfootxyzx(0) = com_inte(0);  
//   _Rfootxyzx(1) = com_inte(1);   
//   _Rfootxyzx(2) = com_inte(2);   
//   
//   return com_inte;
//   
// 
// 
// }





////solve the inverse matrix of 4*4 coefficient matrices
void NLPClass::solve_AAA_inv(Eigen::Matrix<double, 4, 1> t_plan)
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
Eigen::Matrix<double, 7, 7> NLPClass::solve_AAA_inv_x(Eigen::Vector3d t_plan, int j_index)
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
  
//   	  if ((j_index>135)&&(j_index<140))
// 	  {
// 	  cout <<"AAA="<<endl<<AAA<<endl;	  	    
// 	  }
	  
  
  
  
  Eigen::Matrix<double,7,7> AAA_inv = AAA.inverse(); 
  
//   cout <<"t_plan:" <<t_plan<<endl;
//   cout <<"AAA:" <<AAA<<endl;
  
  return AAA_inv;
  
}










//////////////////////////////////////////=============================ZMP optimal distribution for lower level adimittance control=================================
////reference_force_torque_distribution========================================

void NLPClass::Zmp_distributor(int walktime, double dt_sample)
{
// //// judge if stop  
  
  int j_index = floor(walktime / (_dt / dt_sample));
  
  zmp_interpolation(j_index,walktime,dt_sample);  

// reference_force_torque_distribution 
  if (_bjx1 >= 2)
  {
      if (_bjx1 % 2 == 0)           //odd:left support
      {  
	/// right swing
	if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))
	{
	  int nTx_n = round(_tx(_bjx1-1)/_dt);
	  int nTx_n_dsp = round((_tx(_bjx1-1)+_td(_bjx1-1))/_dt);
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
	if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double suppot
	{
	  int nTx_n = round(_tx(_bjx1-1)/_dt);
	  int nTx_n_dsp = round((_tx(_bjx1-1)+_td(_bjx1-1))/_dt);
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
    {/// left swing
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


void NLPClass::zmp_interpolation(int t_int,int walktime, double dt_sample)
{
  //// calculate by the nonlinear model:
  
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


void NLPClass::Force_torque_calculate(Eigen::Vector3d comxyzx1,Eigen::Vector3d comaxyzx1,Eigen::Vector3d thetaaxyx1,Eigen::Vector3d Lfootxyz1,Eigen::Vector3d Rfootxyz1)
{
  Eigen::Vector3d gra;
  gra << 0,0, -_ggg;
  
  Eigen::Vector3d F_total = _robot_mass * (comaxyzx1 - gra);
  
  Eigen::Vector3d the3a;
  the3a << thetaaxyx1(0),thetaaxyx1(1),0;
  Eigen::Vector3d L_total = _j_ini * the3a;
  
  _F_R = _Co_R * F_total;
  _F_L = _Co_L * F_total;
  
  Eigen::Vector3d R_det_foot_com = Rfootxyz1 -  comxyzx1;
  
  Eigen::Vector3d L_det_foot_com = Lfootxyz1 -  comxyzx1;
  
  Eigen::Vector3d M_total = L_total - _F_R.cross(R_det_foot_com) - _F_L.cross(L_det_foot_com);
  
  _M_R = _Co_R*M_total;  
  _M_L = _Co_L*M_total;
  
}

