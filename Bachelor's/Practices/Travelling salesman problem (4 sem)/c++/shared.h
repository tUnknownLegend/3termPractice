#ifndef INC_1_SHARED_H
#define INC_1_SHARED_H

// input for Quick
#define inFileQuick "../test11.txt"
// input for для NonEff
#define inFileNonEff "../test15.txt"
// output for Quick
#define outFileQuick "../out_Quick.txt"
// output for Noneff
#define outFileNonEff "../outTest.txt"
// compare for double
#define COMPARE_RATE 10e-5
// depth of for loop in sim alg. Probably optimal
#define CYCLE_RATE 100
// zero division error
#define DIVISTION_ERROR 5
// colling rate of temperature. Change to optimize results. 0.1 - fastest, 0.(9) - slowest
#define COOLING_RATE 0.9
// temperature max. Probably optimal
#define TEMP_MAX DBL_MAX
// temperature mim. Probably optimal
#define TEMP_MIN 0.001
// depth of for loop to repeat simAnnealing() exec. Change to optimize results, mainly time of execution
#define REPEAT_RATE 10
// depth of for loop to repeat GetSimAnneling(). For testing;
#define LOOP_RATE 100

#endif //INC_1_SHARED_H

#include <vector>
#include <algorithm>

void inputMatrix(std::vector<std::vector<double>> &matrix);

double MeasureFuncExecTime(const std::function<void()> &FuncToMeasure);

//  This function generates a random double in [i, j]
double GetRandomDouble(double i, double j);

//  This function generates a random int in [i, j]
unsigned int GetRandomInt(unsigned int i, unsigned int j);

void outputMatrix(int amtOfVertices);

