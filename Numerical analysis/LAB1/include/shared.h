#ifndef INC_LAB1_SHARED_H
#define INC_LAB1_SHARED_H

// input file for matrix
#define IN_FILE_MATRIX "../../../inputMatrix.txt"
// input file for vector
#define IN_FILE_VECTOR "../../../inputVector.txt"
// output file for matrix
#define OUT_FILE_MATRIX "../../../outputMatrix.txt"
// output file for vector
#define OUT_FILE_VECTOR "../../../outputVector.txt"
// compare for double
#define COMPARE_RATE 10e-8
// zero division error
#define DIVISTION_ERROR 5

#endif //INC_1_SHARED_H

#include <vector>
#include <algorithm>

//  This function generates a random double in [i, j]
double GetRandomDouble(double i, double j);

void inputMatrix(std::vector<std::vector<double>>& matrix);

void outputMatrix(const std::vector<std::vector<double>>& matrix);

void outputMatrix(int amtOfVertices);

void outputVector(int amtOfElements);

void inputVector(std::vector<double>& vect);

void outputVector(const std::vector<double>& vect);
