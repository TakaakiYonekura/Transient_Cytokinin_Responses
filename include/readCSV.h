#pragma once
#include "stdafx.h"

Eigen::MatrixXd readCytokinin(std::string& filename); // (Cytp Cyt)^T
Eigen::MatrixXd readGenes(std::string& filename); // (Cytp Cyt AHKs TCSn Type-A_ARRs) ^T
Eigen::MatrixXd readGenes_Mini(std::string& filename); // (Cytp Cyt AHKs TCSn Type-A_ARRs) ^T
Eigen::MatrixXd readGenes_Minimum(std::string& filename); // (Cytp Cyt AHKs TCSn Type-A_ARRs) ^T

Eigen::MatrixXd expansionMatrix(Eigen::MatrixXd m_input, int expand);