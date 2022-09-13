#include <iostream>
#include "shared.h"

using std::vector;

int main()
{
	vector<vector<double>> matrix;
	vector<double> vect;
	inputMatrix(matrix);
	inputVector(vect);
	outputMatrix(matrix);
	outputVector(vect);
	return 0;
}
