#include <iostream>
#include "shared.h"

using std::vector;

vector<vector<double>> matrixMultiplication(const vector<vector<double>>& firstM, const vector<vector<double>>& secondM) {
	vector<vector<double>> resMatrix(firstM.size(),
		vector<double>(firstM.size(), 0));

	for (int i = 0; i < firstM.size(); ++i) {
		for (int j = 0; j < secondM.size(); ++j) {
			for (int k = 0; k < firstM.size(); ++k) {
				resMatrix[i][j] += firstM[i][k] * secondM[k][j];
			}
		}
	}

	return resMatrix;
}


vector<double> CalcQRmethod() {
	vector<vector<double>> matrix;
	vector<double> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);

	vector<double> resultVect(rightVect.size(), 0.0);

	outputMatrix(matrixMultiplication(matrix, matrix));

	return {};
}

vector<double> getQR() {
	vector<double> res = CalcQRmethod();
	outputVector(res);

	return res;
}