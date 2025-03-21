#pragma once
#include "stdafx.h"

double difAve(Eigen::MatrixXd& control, Eigen::MatrixXd& simulation);

double difAve(Eigen::MatrixXd& control, Eigen::MatrixXd& simulation, Eigen::MatrixXd& weight);


double dif_w(Eigen::ArrayXd& control, Eigen::ArrayXd& simulation, Eigen::ArrayXd& weight);