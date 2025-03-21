#include "stdafx.h"
#include "difAve.h"


double difAve(Eigen::MatrixXd& control, Eigen::MatrixXd& simulation) {
	if (control.rows() != simulation.rows() || control.cols() != simulation.cols()) {
		throw std::invalid_argument("Matrices A and B must have the same dimentions.");
	}
	double dif = 0;
	for (int i = 0; i < control.rows(); i++) {
		double ave = 0;
		double dif_r = 0;
		for (int j = 0; j < control.cols(); j++) {
			double cr = control(i, j) - simulation(i, j);
			dif_r += cr * cr;
			ave += control(i, j);
		}
		ave /= double(control.cols());
		dif += dif_r / (ave * ave);
	}

	return dif;
}


double difAve(Eigen::MatrixXd& control, Eigen::MatrixXd& simulation, Eigen::MatrixXd& weight) {
	if (control.rows() != simulation.rows() || control.cols() != simulation.cols() || control.rows() != weight.rows() || control.cols() != weight.cols()) {
		throw std::invalid_argument("Matrices A and B must have the same dimentions.");
	}
	double dif = 0;
	for (int i = 0; i < control.rows(); i++) {
		double ave = 0;
		double dif_r = 0;
		for (int j = 0; j < control.cols(); j++) {
			double cr = control(i, j) - simulation(i, j);
			dif_r += weight(i, j) * cr * cr;
			ave += control(i, j);
		}
		ave /= double(control.cols());
		dif += dif_r / (ave * ave);
	}

	return dif;
}


double dif_w(Eigen::ArrayXd& control, Eigen::ArrayXd& simulation, Eigen::ArrayXd& weight) {
	if (control.rows() != simulation.rows() || control.cols() != simulation.cols() || control.rows() != weight.rows() || control.cols() != weight.cols()) {
		throw std::invalid_argument("Matrices A and B must have the same dimentions.");
	}
	Eigen::ArrayXd dif = control - simulation;
	dif = weight * dif * dif;

	return dif.sum();
}
