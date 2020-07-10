#include "MPC/MPCClass.h"
#include <iostream>
#include <Eigen/Dense>
#include <math.h>


using namespace Eigen;
using namespace std;


int main()
{
  
        MPCClass mpc;	
        cout << "1"<<endl;
  	int stepnumber =15;
	double stepwidthinput = 0.0726*2;
	double steplengthinput = 0.2;
	double stepheightinput = 0.0;
// 	MPCClass mpc(stepnumber, steplengthinput, stepwidthinput,stepheightinput);
	


	// pass inside the default foot steps
	mpc.FootStepInputs(stepnumber,stepwidthinput,steplengthinput,stepheightinput);
	

// 	// optimize the footsteps,COM aod swing foot trajectory
        mpc.CoM_foot_trajection_generation();
	
	int walktime = 1;
	// get the optimized footsteps
// 	std::vector<Eigen::Vector3d> optmized_footstep = mpc.GetSolution_CoM_position();
	Eigen::Vector3d CoM_position = mpc.GetSolution_CoM_position(walktime );
	Eigen::Vector2d CoM_inclination = mpc.GetSolution_CoM_inclination(walktime );
	Eigen::Vector3d foot_location = mpc.GetSolution_foot_location(walktime );
	Eigen::Vector3d Left_foot_position = mpc.GetSolution_Foot_positionL(walktime );
        Eigen::Vector3d Right_foot_position = mpc.GetSolution_Foot_positionR(walktime );

	// display the optimized footsteps
	std::cout << "Optimized foot location: " <<endl<< foot_location<< std::endl;	
	std::cout << "Optimized CoM_position: " <<endl<< CoM_position<< std::endl;	
	std::cout << "Optimized CoM_inclination: " <<endl<< CoM_inclination<< std::endl;	
	std::cout << "Optimized Left_foot_position: " <<endl<< Left_foot_position<< std::endl;	
	std::cout << "Optimized Right_foot_position: " <<endl<< Right_foot_position<< std::endl;	

	
	return 0;
}




