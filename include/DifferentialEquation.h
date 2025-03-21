#pragma once
#include "stdafx.h"

class ParameterSpace;

class ParameterSpace {
public:
	// Fixed paramters.
	double a[6] = { 0.1894, -0.02866, 0.002669, -0.00005479, 0.0000005124, -0.000000001797 };
	
	// Parameters needed optimization.
	double AHKs_ini = 60;
	double AHKs_P_iniratio = 0.0;

	double TCSn_ini = 30;

	double BARRs_tot = 160;
	double BARRs_P_iniratio = 0;

	double ARRs_ini = 100;
	double ARRs_P_iniratio = 0;

	double p_AHKs = 6;
	double p_AHKs_base = 0.36;
	double t_AHKs = 15;
	double d_AHKs = 0.025;

	double t_ARRs = 0.5;

	double p_ARRs = 6;
	double d_ARRs = 0.06;
	double d_ARRs_P = 1;

	double KCyt = 0.2;
	double Kbp = 0.01;
	double Km = 5;
	double Ka = 240;

	double p_TCSn = 3.0;
	double d_TCSn = 0.1;

	double dt = 0.001;

	ParameterSpace() {}
	ParameterSpace(std::vector<double> Parameters);



	bool ode(std::vector<double> &dstate, std::vector<double> &state, double &t);
	bool RK4_singlestep(std::vector<double>& state, double& t, bool (ParameterSpace::*ode)(std::vector<double>&,std::vector<double>&,double&));
	bool RK4_singlestep(std::vector<double>& state, double& t, double& maxdifratio, std::vector<double>& dstate, bool (ParameterSpace::* ode)(std::vector<double>&, std::vector<double>&, double&));

	bool MakeOutput(std::vector<double> timepoints, Eigen::MatrixXd& output);
	bool MakeOutput(std::vector<double> timepoints, Eigen::ArrayXd& output);
	bool MakeOutput(std::vector<double> timepoints, Eigen::ArrayXd& output, double& maxdifratio);
};

double objective_function_all(unsigned n, const double* Parameters, double* grad, void* data);
