#include <fstream>
#include <iostream>
#include <algorithm>
#include "shared.h"

using std::ifstream;
using std::vector;
using std::cerr;
using std::ofstream;
using std::cout;

void inputMatrix(vector<vector<double>> &matrix) {
    ifstream in_file(inFileNonEff);
    if (!in_file.is_open()) {
        cerr << "error // input.txt open\n";
        return;
    }

    int amtOfVertices = 0;
    // reading from file
    in_file >> amtOfVertices;

    //cout << "Start reading\n";
    {
        vector<double> str;
        double node = 0.0;
        for (int i = 0; i < amtOfVertices; ++i) {

            for (int j = 0; j < amtOfVertices; ++j) {
                in_file >> node;
                str.push_back(node);
            }
            matrix.push_back(str);
            str.clear();
        }
    }
    //cout << "End reading\n";
}

double MeasureFuncExecTime(const std::function<void()>& FuncToMeasure){
    unsigned int startingTime = clock();
    FuncToMeasure();
    unsigned int stopTime = clock();
    unsigned int searchTime = stopTime - startingTime;   //  exec time
    //cout << "\nSearch time: " << ((double) searchTime) / CLOCKS_PER_SEC << "\n";

    return (((double) searchTime) / CLOCKS_PER_SEC);
}