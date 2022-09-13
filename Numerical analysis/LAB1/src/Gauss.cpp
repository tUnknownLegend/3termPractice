#include <iostream>
#include "shared.h"

using std::vector;

vector<double> CalcGaussMethod() {
	vector<vector<double>> matrix;
	vector<double> rightVect;
	inputMatrix(matrix);
	inputVector(rightVect);

	vector<double> resultVect(rightVect.size(), 0.0);

	for (int k = 0; k < rightVect.size(); ++k) {
		//  partial selection
		double maxValInd = k;
		for (int i = k + 1; i < rightVect.size(); ++i) {
			maxValInd = matrix[i][k] > matrix[maxValInd][k] ? i : maxValInd;
		}
		
		if (abs(matrix[maxValInd][k]) > COMPARE_RATE) {
			std::swap(matrix[maxValInd], matrix[k]);
			std::swap(rightVect[maxValInd], rightVect[k]);

			//  straight
			for (int i = k + 1; i < rightVect.size(); ++i) {
				double c = (matrix[i][k] / matrix[k][k]);
				for (int j = k; j < rightVect.size(); ++j) {
					matrix[i][j] -= c * matrix[k][j];
				}

				rightVect[i] -= c * rightVect[k];
			}
		}
		else {
			std::cerr << "Matrix is singular";
			return {};
		}
	}

	//  reverse
	for (int i = rightVect.size() - 1; i >= 0; --i) {
		double sumJ = 0.0;
		for (int j = i + 1; j < rightVect.size(); ++j) {
			sumJ += matrix[i][j] * resultVect[j];
		}

		resultVect[i] = (rightVect[i] - sumJ) / matrix[i][i];
	}

	outputMatrix(matrix);
	//outputMatrix(10);
	//outputVector(resultVect);
	return resultVect;
	//outputVector(10);
}

vector<double> getGauss() {
	vector<double> res = CalcGaussMethod();
	outputVector(res);

	return res;
}