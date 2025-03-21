#include "stdafx.h"
#include "DifferentialEquation.h"
#include "readCSV.h"

int main(void) {


	
	std::string all_file = "fitting_genes_revise_simple.csv";

	Eigen::MatrixXd all = readGenes_Minimum(all_file);

	std::cout << "KOKO" << std::endl;


	std::cout << all << std::endl;

	std::vector<std::vector<double>> all_vec(/*expand_*/all.rows(), std::vector<double>(all.cols()));

	double* grad = nullptr;
	std::vector<double> grad2;
	// BARRs_tot, ARRs_ini, AHKs_ini, TCSn_ini, AHKs_P_iniratio, ARRs_P_iniratio,
	// p_AHKs, p_AHKs_base, t_AHKs, d_AHKs,
	// t_ARRs, d_ARRs_P/t_ARRs, p_ARRs, d_ARRs/p_ARRs,
	// KCyt, Kbp, Km, Ka,
	// p_TCSn, d_TCSn
	
	double AllParameters[20] =
	{
		76.5287,50.0592,60,30,0.000259232,0.00023976,
		253.271,0.000312319,19.1776,0.837665,
		0.0923777,1.71492,8.45846,6.67826e-08,
		2.57727e-05,85.7202,0.00310143,37.1887,
		2.82286,0.10176
	};
	std::vector<double> AllModelParameters(std::begin(AllParameters), std::end(AllParameters));
	

	for (int i = 0; i < all.rows(); ++i) {
		for (int j = 0; j < all.cols(); ++j) {
			all_vec[i][j] = all(i, j);
		}
	}
	
	ParameterSpace ps = ParameterSpace::ParameterSpace(AllModelParameters);
	std::vector<double> timepoint;
	for (int i = 0; i <= 120; i += 24/*2*/) {
		timepoint.push_back(double(i));
	}
	Eigen::MatrixXd sim_data = Eigen::MatrixXd(3, 6);
	ps.MakeOutput(timepoint, sim_data);

	std::cout << "sim_data" << std::endl;
	std::cout << sim_data << std::endl;

	std::cout << "Residue:" << objective_function_all(20, AllParameters, grad, &all_vec) << std::endl;

	nlopt::opt opt(nlopt::GN_ISRES, 20);  // optimization methods.
	opt.set_min_objective(objective_function_all, &all_vec);

	// Find optimized parameters with non-negative values.
	//
	// BARRs_tot, ARRs_ini, AHKs_P_iniratio, ARRs_P_iniratio,
	// p_AHKs, p_AHKs_base, t_AHKs, d_AHKs,
	// t_ARRs, d_ARRs_P/t_ARRs, p_ARRs, d_ARRs/p_ARRs,
	// KCyt, Kbp, Km, Ka,
	// p_TCSn, d_TCSn
	std::vector<double> lower_bounds =
	{ 40, 50, 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0, 
		0, 0, 0, 0,
		0, 0
	};
	std::vector<double> upper_bounds =
	{ 170, 150, 100, 100, 1, 1,
		10000, 100, 10000, 100,
		100, 100, 10000, 100, 
		1000, 100, 1000, 100,
		10000, 100
	};
	opt.set_lower_bounds(lower_bounds);
	opt.set_upper_bounds(upper_bounds);

	opt.set_xtol_rel(1e-4);

	double minf;
	std::vector<double> optimized_parameters(std::begin(AllParameters), std::end(AllParameters));
		
	try {
		nlopt::result result = opt.optimize(optimized_parameters, minf);
		std::cout << "Optimal coefficients: " << std::endl;
		for (size_t j = 0; j < optimized_parameters.size(); j++) {
			std::cout << optimized_parameters[j] << std::endl;
		}

	}
	catch (std::exception& e) {
		std::cerr << "Optimization failed: " << e.what() << "\n";
	}
	
	std::ofstream output_file("optimized_parameters.csv");
	if (output_file.is_open()) {
		for (size_t i = 0; i < optimized_parameters.size(); ++i) {
			output_file << optimized_parameters[i];
			if (i < optimized_parameters.size() - 1) {
				output_file << ",";
			}
		}
		output_file << "\n";
		output_file.close();
		std::cout << "optimized_parameters has been written to optimized_parameters.csv" << std::endl;
	}
	else {
		std::cerr << "Failed to open the file for writing." << std::endl;
	}
	
	return 0;
}