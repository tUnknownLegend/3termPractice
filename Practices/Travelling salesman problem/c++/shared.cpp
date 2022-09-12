#include <fstream>
#include <iostream>
#include <algorithm>
#include "bits/stdc++.h"
#include "shared.h"

using std::ifstream;
using std::vector;
using std::cerr;
using std::ofstream;
using std::cout;

//  This function generates a random double in [i, j]
double GetRandomDouble(double i, double j) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(i, j);
    return dis(gen);
}

//  This function generates a random int in [i, j]
unsigned int GetRandomInt(unsigned int i, unsigned int j) {
    std::random_device rd;  //  Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //  Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(i, j);
    return distrib(gen);
}

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

double MeasureFuncExecTime(const std::function<void()> &FuncToMeasure) {
    unsigned int startingTime = clock();
    FuncToMeasure();
    unsigned int stopTime = clock();
    unsigned int searchTime = stopTime - startingTime;   //  exec time
    //cout << "\nSearch time: " << ((double) searchTime) / CLOCKS_PER_SEC << "\n";

    return (((double) searchTime) / CLOCKS_PER_SEC);
}

void outputMatrix(int amtOfVertices) {
    ofstream out_file(outFileNonEff);
    if (!out_file.is_open()) {
        cerr << "error // output.txt open\n";
        return;
    }

    // reading from file
    out_file << amtOfVertices << std::endl;

    //cout << "Start reading\n";
    {
        double node = 0.0;
        for (int i = 0; i < amtOfVertices; ++i) {

            for (int j = 0; j < amtOfVertices; ++j) {
                if (i == j)
                    out_file << 0.0 << " ";
                else
                    out_file << GetRandomDouble(i, j) << " ";
            }
            out_file << std::endl;
        }
    }
    //cout << "End reading\n";
}