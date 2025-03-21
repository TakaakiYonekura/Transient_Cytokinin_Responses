#include "stdafx.h"
#include "readCSV.h"

Eigen::MatrixXd readCytokinin(std::string& filename) {
	std::ifstream file(filename);
	std::string line;

	Eigen::MatrixXd data = Eigen::MatrixXd(2, 6);
	int i = 0;
	while (std::getline(file, line)) {
		std::vector<std::string> row;
		std::stringstream linestream(line);
		std::string cell;
		while (std::getline(linestream, cell, ',')) {
			row.push_back(cell);
		}
		if (i == 0) {
			for (int j = 0; j < 6; j++) {
				data(1, j) = stod(row[j + 1]);
			}
		}
		else if (i == 2) {
			for (int j = 0; j < 6; j++) {
				data(0, j) = stod(row[j + 1]);
			}
		}
		i++;

	}

	return data;
}


Eigen::MatrixXd readGenes(std::string& filename) {
	std::ifstream file(filename);
	std::string line;

	Eigen::MatrixXd data = Eigen::MatrixXd(5, 6);
	int i = 0;
	while (std::getline(file, line)) {
		std::vector<std::string> row;
		std::stringstream linestream(line);
		std::string cell;
		while (std::getline(linestream, cell, ',')) {
			row.push_back(cell);
		}
		for (int j = 0; j < 6; j++) {
			data(i, j) = stod(row[j + 1]);
		}
		i++;

	}

	return data;
}

Eigen::MatrixXd readGenes_Mini(std::string& filename) {
	std::ifstream file(filename);
	std::string line;

	Eigen::MatrixXd data = Eigen::MatrixXd(4, 6);
	int i = 0;
	while (std::getline(file, line)) {
		std::vector<std::string> row;
		std::stringstream linestream(line);
		std::string cell;
		while (std::getline(linestream, cell, ',')) {
			row.push_back(cell);
		}
		for (int j = 0; j < 6; j++) {
			data(i, j) = stod(row[j + 1]);
		}
		i++;

	}

	return data;
}


Eigen::MatrixXd readGenes_Minimum(std::string& filename) {
	std::ifstream file(filename);
	std::string line;

	Eigen::MatrixXd data = Eigen::MatrixXd(3, 6);
	int i = 0;
	while (std::getline(file, line)) {
		std::vector<std::string> row;
		std::stringstream linestream(line);
		std::string cell;
		while (std::getline(linestream, cell, ',')) {
			row.push_back(cell);
		}
		for (int j = 0; j < 6; j++) {
			data(i, j) = stod(row[j + 1]);
		}
		i++;

	}

	return data;
}



Eigen::MatrixXd expansionMatrix(Eigen::MatrixXd m_input, int expand) {
	size_t new_cols = 1 + (m_input.cols() - 1) * expand;
	Eigen::MatrixXd m_output = Eigen::MatrixXd(m_input.rows(), new_cols);
	for (size_t i = 0; i < m_input.cols() - 1; i++) {
		for (size_t j = 0; j < m_input.rows(); j++) {
			for (int k = 0; k < expand; k++) {
				m_output(j, i * expand + k) = m_input(j, i) + (m_input(j, i + 1) - m_input(j, i)) * k / double(expand);
			}
		}
	}
	for (size_t j = 0; j < m_input.rows(); j++) {
		m_output(j, m_input.cols() * expand) = m_input(j, m_input.cols());
	}
	return m_output;
}
