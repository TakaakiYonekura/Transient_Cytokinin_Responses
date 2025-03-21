#include "stdafx.h"
#include "DifferentialEquation.h"
#include "difAve.h"



ParameterSpace::ParameterSpace(std::vector<double> Parameters) {
	if (Parameters.size() != 20)std::cout << "ParameterSpace.Parameters.size() should be 18." << std::endl;
	else {
		// BARRs_tot, AHKs_P_iniratio, TARRs_P_iniratio,
		// p_AHKs, p_AHKs_base, t_AHKs, d_AHKs,
		// t_ARRs, d_ARRs_P/t_ARRs, p_ARRs, d_ARRs/p_ARRs,
		// KCyt, Kbp, Km, Ka,
		// p_LDBs, d_LBDs/p_LDBs, KL

		BARRs_tot = Parameters[0];
		ARRs_ini = Parameters[1];
		AHKs_ini = Parameters[2];
		TCSn_ini = Parameters[3];


		AHKs_P_iniratio = Parameters[4];
		ARRs_P_iniratio = Parameters[5];

		p_AHKs = Parameters[6];
		p_AHKs_base = Parameters[7];
		t_AHKs = Parameters[8];
		d_AHKs = Parameters[9];

		t_ARRs = Parameters[10];
		d_ARRs_P = t_ARRs * Parameters[11];

		p_ARRs = Parameters[12];
		d_ARRs = p_ARRs * Parameters[13];

		KCyt = Parameters[14];
		Kbp = Parameters[15];
		Km = Parameters[16];
		Ka = Parameters[17];

		p_TCSn = Parameters[18];
		d_TCSn = Parameters[19];

	}


}

bool ParameterSpace::ode(std::vector<double>& dstate, std::vector<double>& state, double& t) {
	if (dstate.size() != 6 || state.size() != 6) return false;
	//state = BARRs_P, AHKs, AHKs_P, ARRs, ARRs_P, TCSn.
	double& BARRs_P = state[0];
	double BARRs = BARRs_tot - state[0];
	double& AHKs = state[1];
	double& AHKs_P = state[2];
	double& ARRs = state[3];
	double& ARRs_P = state[4];
	double& TCSn = state[5];

	double& dBARRs_P_dt = dstate[0];
	double& dAHKs_dt = dstate[1];
	double& dAHKs_P_dt = dstate[2];
	double& dARRs_dt = dstate[3];
	double& dARRs_P_dt = dstate[4];
	double& dTCSn_dt = dstate[5];

	double Cyt = a[0] + a[1] * t + a[2] * (t * t) + a[3] * (t * t * t) + a[4] * (t * t * t * t) + a[5] * (t * t * t * t * t);
	double BARRs_P_act = BARRs_P;// / (1 + LBDs / Ki);

	dBARRs_P_dt = -d_ARRs_P * BARRs_P + t_ARRs * BARRs * AHKs_P / (Km + BARRs + ARRs) ;

	dAHKs_P_dt = t_AHKs * Cyt / (KCyt + Cyt) * AHKs - (t_ARRs * BARRs + t_ARRs * ARRs) * AHKs_P / (Km + BARRs + ARRs) ;
	dAHKs_dt = p_AHKs * BARRs_P_act / (Kbp + BARRs_P_act) + p_AHKs_base - dAHKs_P_dt - d_AHKs * (AHKs + AHKs_P) * AHKs;

	dARRs_P_dt = t_ARRs * ARRs * AHKs_P / (Km + BARRs + ARRs)- d_ARRs_P * ARRs_P;
	dARRs_dt = p_ARRs * BARRs_P_act / (Ka + BARRs_P_act) - dARRs_P_dt - d_ARRs * ARRs;

	dTCSn_dt = p_TCSn * BARRs_P_act - d_TCSn * TCSn;


	return true;
}


bool ParameterSpace::RK4_singlestep(std::vector<double>& state, double& t, bool (ParameterSpace::* ode)(std::vector<double>&, std::vector<double>&, double&)) {

	size_t n = state.size();
	std::vector<double> k1(n), k2(n), k3(n), k4(n), state_temp(n);

	double dht = 0.5 * dt;

	for (size_t i = 0; i < n; i++) state_temp[i] = state[i];
	if (!(this->*ode)(k1, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) state_temp[i] = state[i] + dht * k1[i];
	t += dht;
	if (!(this->*ode)(k2, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) state_temp[i] = state[i] + dht * k2[i];
	if (!(this->*ode)(k3, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) state_temp[i] = state[i] + dt * k3[i];
	t += dht;
	if (!(this->*ode)(k4, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) state[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * dt / 6.0;

	for (size_t i = 0; i < n; i++) {
		if (state[i] < 0) state[i] = 0;
	}

	return true;

}

bool ParameterSpace::RK4_singlestep(std::vector<double>& state, double& t, double& maxdifratio, std::vector<double>& dstate, bool (ParameterSpace::* ode)(std::vector<double>&, std::vector<double>&, double&)) {

	size_t n = state.size();
	std::vector<double> k1(n), k2(n), k3(n), k4(n), state_temp(n);

	double dht = 0.5 * dt;

	for (size_t i = 0; i < n; i++) state_temp[i] = state[i];
	if (!(this->*ode)(k1, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) state_temp[i] = state[i] + dht * k1[i];
	t += dht;
	if (!(this->*ode)(k2, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) state_temp[i] = state[i] + dht * k2[i];
	if (!(this->*ode)(k3, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) state_temp[i] = state[i] + dt * k3[i];
	t += dht;
	if (!(this->*ode)(k4, state_temp, t)) return false;
	for (size_t i = 0; i < n; i++) {
		double dif = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
		if (t != dt) {
			if (fabs(dif - dstate[i]) / dt > maxdifratio) {
				maxdifratio = fabs(dif - dstate[i]) / dt;
			}
		}
		dstate[i] = dif;
		state[i] += (dif * dt);
	}

	for (size_t i = 0; i < n; i++) {
		if (state[i] < 0) state[i] = 0;
	}

	return true;

}


bool ParameterSpace::MakeOutput(std::vector<double> timepoints, Eigen::MatrixXd& output) {
	if (output.cols() != timepoints.size() || output.rows() != 3) return false;
	double t = 0;
	std::vector<double> state = { BARRs_tot * BARRs_P_iniratio, AHKs_ini * (1 - AHKs_P_iniratio), AHKs_ini * AHKs_P_iniratio, ARRs_ini * (1 - ARRs_P_iniratio), ARRs_ini * ARRs_P_iniratio, TCSn_ini };
	int i = 0;
	for (double time : timepoints) {
		while (t < time) {
			RK4_singlestep(state, t, &ParameterSpace::ode);
		}
		output(0, i) = state[1] + state[2];
		output(1, i) = state[5];
		output(2, i) = state[3] + state[4];
		i++;
	}

	return true;

}

bool ParameterSpace::MakeOutput(std::vector<double> timepoints, Eigen::ArrayXd& output) {
	if (output.size() != 3 * timepoints.size()) return false;
	double t = 0;
	std::vector<double> state = { BARRs_tot * BARRs_P_iniratio, AHKs_ini * (1 - AHKs_P_iniratio), AHKs_ini * AHKs_P_iniratio, ARRs_ini * (1 - ARRs_P_iniratio), ARRs_ini * ARRs_P_iniratio, TCSn_ini };
	int i = 0;
	int c = timepoints.size();
	for (double time : timepoints) {
		while (t < time) {
			RK4_singlestep(state, t, &ParameterSpace::ode);
		}
		output(i) = state[1] + state[2];
		output(c + i) = state[5];
		output(2 * c + i) = state[3] + state[4];
		i++;
	}

	return true;

}

bool ParameterSpace::MakeOutput(std::vector<double> timepoints, Eigen::ArrayXd& output, double& maxdifratio) {
	if (output.size() != 3 * timepoints.size()) return false;
	double t = 0;
	std::vector<double> state = { BARRs_tot * BARRs_P_iniratio, AHKs_ini * (1 - AHKs_P_iniratio), AHKs_ini * AHKs_P_iniratio, ARRs_ini * (1 - ARRs_P_iniratio), ARRs_ini * ARRs_P_iniratio, TCSn_ini };
	int i = 0;
	int c = timepoints.size();
	std::vector<double> dstate = { 0,0,0,0,0,0 };
	for (double time : timepoints) {
		while (t < time) {
			RK4_singlestep(state, t, maxdifratio, dstate, &ParameterSpace::ode);
		}
		output(i) = state[1] + state[2];
		output(c + i) = state[5];
		output(2 * c + i) = state[3] + state[4];
		i++;
	}

	return true;

}


double objective_function_all(unsigned n, const double* Parameters, double* grad, void* data) {

	std::vector<std::vector<double>>* observed_data = reinterpret_cast<std::vector<std::vector<double>>*>(data);
	std::vector<double> theta(Parameters, Parameters + n);

	size_t rows = observed_data->size();
	size_t cols = ((*observed_data)[0]).size();

	Eigen::ArrayXd o_matrix(rows * cols);

	// Building flatten data.
	for (size_t i = 0; i < rows; ++i) {
		if ((*observed_data)[i].size() != cols) {
			throw std::runtime_error("Inconsistent row sizes in input data.");
		}
		for (size_t j = 0; j < cols; ++j) {
			o_matrix(i * cols + j) = (*observed_data)[i][j];
		}
	}

	std::vector<double> timepoints;
	for (int i = 0; i <= 120; i += 24) {
		timepoints.push_back(double(i));
	}
	ParameterSpace ps = ParameterSpace(theta);


	Eigen::ArrayXd all_simdata = Eigen::ArrayXd(3 * 6);
	if (!ps.MakeOutput(timepoints, all_simdata)) std::cout << "ERROR" << std::endl;

	Eigen::ArrayXd weight_matrix = Eigen::ArrayXd(3 * 6);
	weight_matrix <<
		1, 1, 0.1, 
		1, 1, 0.5, 
		1, 1, 1, 
		1, 1, 1, 
		1, 1, 0.5, 
		1, 1, 0.1;

	double dif = ((o_matrix - all_simdata) * (o_matrix - all_simdata) * weight_matrix).sum();

	dif += 0.8 * (ps.BARRs_tot - 160) * (ps.BARRs_tot - 160);

	dif += std::accumulate(theta.begin() + 4, theta.end(), 0.0,
		[](double sum, double value) {
			return sum + 0.2 * value * value;
		});

	for (size_t i = 0; i < theta.size(); i++) {
		std::cout << theta[i] << ",";
	}

	std::cout << dif << std::endl;


	if (std::isnan(dif) || std::isinf(dif)) {
		return std::numeric_limits<double>::infinity();
	}

	return (dif);

}
