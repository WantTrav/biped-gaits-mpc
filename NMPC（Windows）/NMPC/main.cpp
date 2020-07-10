#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<Eigen/Dense>
#include "test_class.h"
//#include "RobotState\RobotStateClass.h"
//#include "MPC\MPCClass.h"
#include "RTControl\MpcRTControlClass.h"
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>

using namespace std;


void main()
{
	
	clock_t t_start,t_finish;
	t_start = clock();

	MpcRTControlClass rt_control;


	int totnumber = rt_control._walkdtime_max;
	
	Eigen::MatrixXd nlp_result;
	nlp_result.setZero(27,totnumber);
	
	
	for (int j =0; j<totnumber;j++)
	{
	  
	  //nlp_result.col(j) = rt_control.WalkingReactStepping(j);
	}
	
	
	cout<< "save file!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	

    t_finish = clock();
	double _tcpu_prepara = (double)(t_finish - t_start)/CLOCKS_PER_SEC ;
	cout<<"time cost:"<<"\n"<<_tcpu_prepara<<endl;
	system("pause");
}