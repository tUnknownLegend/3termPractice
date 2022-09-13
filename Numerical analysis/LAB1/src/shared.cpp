#include <fstream>
#include <iostream>
#include <random>
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

void inputMatrix(vector<vector<double>>& matrix) {
	ifstream inFile(IN_FILE_MATRIX);
	if (!inFile.is_open()) {
		cerr << "error // input.txt open\n";
		return;
	}

	int amtOfVertices = 0;
	inFile >> amtOfVertices;
	matrix.reserve(amtOfVertices);

	{
		vector<double> str;
		double node = 0.0;
		for (int i = 0; i < amtOfVertices; ++i) {

			for (int j = 0; j < amtOfVertices; ++j) {
				inFile >> node;
				str.push_back(node);
			}
			matrix.push_back(str);
			str.clear();
		}
	}
	inFile.close();
}

void inputVector(vector<double>& vect) {
	ifstream inFile(IN_FILE_VECTOR);
	if (!inFile.is_open()) {
		cerr << "error // input.txt open\n";
		return;
	}

	int amtOfVertices = 0;
	inFile >> amtOfVertices;
	vect.reserve(amtOfVertices);

	{
		vector<double> str;
		double node = 0.0;
		for (int j = 0; j < amtOfVertices; ++j) {
			inFile >> node;
			str.push_back(node);
		}
		vect = std::move(str);
	}
	inFile.close();
}

void outputVector(int amtOfElements) {
	ofstream outFile(OUT_FILE_VECTOR);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << amtOfElements << std::endl;

	{
		const int leftBound = 1;
		const int rightBound = 10;
		double node = 0.0;
		for (int j = 0; j < amtOfElements; ++j) {
			outFile << GetRandomDouble(leftBound, rightBound) << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
}

void outputVector(const vector<double>& vect) {
	ofstream outFile(OUT_FILE_VECTOR);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << vect.size() << std::endl;

	{
		double node = 0.0;
		for (auto& el : vect) {
			outFile << el << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
}

void outputMatrix(const vector<vector<double>>& matrix) {
	ofstream outFile(OUT_FILE_MATRIX);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << matrix.size() << std::endl;

	{
		const int leftBound = 0;
		const int rightBound = 10;
		double node = 0.0;
		for (auto& raw : matrix) {
			for (auto& el : raw) {
				outFile << el << " ";
			}
			outFile << std::endl;
		}
	}
	outFile.close();
}

void outputMatrix(int amtOfVertices) {
	ofstream outFile(OUT_FILE_MATRIX);
	if (!outFile.is_open()) {
		cerr << "error // output.txt open\n";
		return;
	}

	outFile << amtOfVertices << std::endl;

	{
		const int leftBound = 0;
		const int rightBound = 10;
		double node = 0.0;
		for (int i = 0; i < amtOfVertices; ++i) {

			for (int j = 0; j < amtOfVertices; ++j) {
				outFile << GetRandomDouble(leftBound, rightBound) << " ";
			}
			outFile << std::endl;
		}
	}
	outFile.close();
}