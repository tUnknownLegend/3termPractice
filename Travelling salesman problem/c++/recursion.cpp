#include <iostream>
#include <fstream>
#include <limits>
#include "shared.h"

using std::vector;

double CalcRec() {
    // matrix
    vector<vector<double>> matrix;
    inputMatrix(matrix);

    // creating sub matrix
    vector<unsigned int> vertexes = {};
    vertexes.reserve(matrix.size());
    for (int i = 1; i < matrix.size(); ++i) {
        vertexes.push_back(i);
    }

    double minWeight = std::numeric_limits<double>::max();

    do {
        double currWeight = 0.0;
        unsigned int currVertex = 0;
        for (auto &i: vertexes) {
            currWeight += matrix[currVertex][i];
            currVertex = i;
        }
        currWeight += matrix[currVertex][0];

        minWeight = std::min(currWeight, minWeight);

    } while (std::next_permutation(vertexes.begin(), vertexes.end()));

    return minWeight;
}

int main() {
    double time = MeasureFuncExecTime([]() { std::cout << CalcRec();  });
    std::cout << "\nSearch time: " << time;

    return 0;
}