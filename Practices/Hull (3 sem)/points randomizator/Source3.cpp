#ifndef DataStructures
#define DataStructures

#define FILE_COORDINATES_PATH "cor.txt"

#define FILE_SAVE_PATH_BRUTE_FORCE "ss.txt"

#define FILE_SAVE_PATH_KIRKPATRICK_SEIDEL "SAVE.txt"

#include <iostream>
#include <vector>
#include "Source.cpp"
#include "Source2.cpp"
using namespace std;

enum Result { isEmpty, isFilled };
enum PointPosition { oneSide, otherSide, within };

#endif