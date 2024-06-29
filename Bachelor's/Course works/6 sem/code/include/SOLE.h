#ifndef INC_LAB_SOLE_H
#define INC_LAB_SOLE_H

#include <vector>
#include "shared.h"

//const std::vector<TT> initPoints = {1.0, 0.0};
const std::vector<TT> initPoints = {1.0, 0.0, 0.0};
//const std::pair<TT, TT> range = {0.0, 1.0};
const std::pair<TT, TT> range = {0.0, 40.0};

const unsigned int multiplayer = 1;

const int numOfPoints = 10000 * (int) range.second * multiplayer;
const TT step = TT(range.second - range.first) / TT(numOfPoints);
const TT tau = step;

#endif //INC_LAB_SOLE_H
