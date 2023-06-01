#include <cmath>
#include <iostream>
#include <cassert>
#include "shared.h"
#include "nonLinearSolve.h"
#include "SOLE.h"

using namespace std;


enum calcMethod {
    MexplicitEuler,
    MimplicitEuler,
    Msymmetric,
    MrungeKutta2,
    MrungeKutta4,
    MrungeKuttaRungeStep,
    MexplicitAdams,
    Mbdf2,
    Mbdf4
};

vector<vector<TT>> calcDiff(const vector<vector<TT>> &answer) {
    vector<vector<TT>> diff(answer);

    for (int i = 0; i < answer.size(); ++i) {
        diff[i][0] = fabs(cos(i * step) - answer[i][0]);
        diff[i][1] = fabs(-sin(i * step) - answer[i][1]);
    }
    return diff;
}

TT getMaxDiff(const vector<vector<TT>> &diff) {
    TT maxValue = 0.0;
    for (const auto &raw: diff) {
        for (const auto &el: raw) {
            maxValue = max(maxValue, el);
        }
    }

    return maxValue;
}

vector<vector<TT>> explicitEuler(const vector<TT> &cond, int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; ++i) {
        y[i + 1] = vectorOperation(y[i],
                                   vectorRDigit(step, f(y[i]), '*'), '+');
    }
    return y;
}

vector<vector<TT>> implicitEuler(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; ++i) {
        y[i + 1] = Newton("Euler", cond.size(), {y[i]});
    }
    return y;
}

vector<vector<TT>> symmetric(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    for (int i = 0; i < n - 1; ++i) {
        y[i + 1] = Newton("Symmetric", cond.size(), {y[i]});
    }
    return y;
}

void k_iSum(vector<TT> &temp, const vector<TT> &k1, const vector<TT> &k2, const vector<TT> &k3,
            const vector<TT> &k4, const vector<TT> &y_i, const TT Tchange) {
    temp = vectorOperation(vectorRDigit(2.0, k2, '*'), vectorRDigit(2.0, k3, '*'), '+');
    temp = vectorOperation(k1, temp, '+');
    temp = vectorOperation(k4, temp, '+');
    temp = vectorRDigit(Tchange / 6.0, temp, '*');
    temp = vectorOperation(y_i, temp, '+');
}

vector<TT> rungeKutta4Calc(const vector<TT> &yi, vector<TT> &tempVector, const TT localStep, const TT Tchange) {
    vector<TT> resVector(yi.size());
    // calc k_i
    vector<TT> k1(f(yi));
    resVector = vectorOperation(tempVector, vectorRDigit(0.5 * Tchange * localStep, k1, '*'), '+');
    vector<TT> k2(f(resVector));
    resVector = vectorOperation(tempVector, vectorRDigit(0.5 * Tchange * localStep, k2, '*'), '+');
    vector<TT> k3(f(resVector));
    resVector = vectorOperation(tempVector, vectorRDigit(Tchange * localStep, k3, '*'), '+');
    vector<TT> k4(f(resVector));

    // sum k_i
    k_iSum(resVector, k1, k2, k3, k4, yi, Tchange);

    return resVector;
}

vector<vector<TT>> rungeKutta2(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = cond;

    vector<TT> k1(numOfPoints);
    vector<TT> k2(numOfPoints);
    for (int i = 0; i < n - 1; ++i) {
        k1 = f(y[i]);
        vector<TT> temp = vectorOperation(y[i],
                                          vectorRDigit(step, k1, '*'), '+');
        k2 = f(temp);
        temp = vectorOperation(k1, k2, '+');
        y[i + 1] = vectorOperation(y[i], vectorRDigit((step / 2.0), temp, '*'), '+');
    }
    return y;
}

vector<vector<TT>> rungeKutta4(const vector<TT> &cond) {
    vector<TT> previous(cond.size());
    vector<TT> help(cond.size());
    vector<vector<TT>> K(4);

    vector<vector<TT>> x;
    x.emplace_back();
    for (const TT i: cond) {
        x[0].push_back(i);
    }

    for (int i = 0; i < 4; ++i) { K[i] = help; };
    for (int t = 0; t < numOfPoints; ++t) {
        previous = x[t];
        x.emplace_back();

        K[0] = f(previous);

        help = vectorOperation(previous, vectorRDigit(step / 2, K[0], '*'), '+');

        K[1] = f(help);

        help = vectorOperation(previous, vectorRDigit(step / 2, K[1], '*'), '+');

        K[2] = f(help);

        help = vectorOperation(previous, vectorRDigit(step, K[2], '*'), '+');

        K[3] = f(help);

        for (int j = 0; j < cond.size(); ++j) {
            x[t + 1].push_back(previous[j] + step / 6 * (K[0][j] + 2 * K[1][j] + 2 * K[2][j] + K[3][j]));
        }
    }

    return x;
}

vector<vector<TT>> bdf2(const vector<TT> &cond, int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = rungeKutta2(cond, numOfPoints)[0];
    y[1] = rungeKutta2(cond, numOfPoints)[1];

    for (int i = 1; i < n - 1; ++i) {
        y[i + 1] = Newton("BDF2", cond.size(), {y[i - 1], y[i]});
    }
    return y;
}

vector<vector<TT>> bdf4(const vector<TT> &cond, int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    {
        const auto firstYs = rungeKutta4(cond);
        for (int i = 0; i < 4; ++i) {
            y[i] = firstYs[i];
        }
    }

    for (int i = 3; i < n - 1; ++i) {
        y[i + 1] = Newton("BDF4", cond.size(), {y[i - 3], y[i - 2], y[i - 1], y[i]});
    }
    return y;
}

vector<vector<TT>> Adams(const vector<TT> &cond, const int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    vector<vector<TT>> y4(n, vector<TT>(cond.size()));
    vector<TT> temp(cond.size());

    y4 = rungeKutta4(cond);
    for (int i = 0; i < 4; ++i) {
        y[i] = y4[i];
    }

    for (int i = 3; i < n - 1; ++i) {
        for (int j = 0; j < cond.size(); ++j) {
            temp = vectorOperation(vectorRDigit(55.0, f(y[i]), '*'),
                                   vectorRDigit(59.0, f(y[i - 1]), '*'), '-');
        }
        temp = vectorOperation(temp,
                               vectorRDigit(37.0, f(y[i - 2]), '*'), '+');

        temp = vectorOperation(temp,
                               vectorRDigit(9.0, f(y[i - 3]), '*'), '-');

        y[i + 1] = vectorOperation(y[i],
                                   vectorRDigit(step,
                                                vectorRDigit(24.0, temp, '*'),
                                                '/'),
                                   '+');
    }
    return y;
}

void templateOutput(const calcMethod method) {
    vector<vector<TT>> result;

    switch (method) {
        case MexplicitEuler:
            result = explicitEuler(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMexplicitEuler.txt");
            break;
        case MimplicitEuler:
            result = implicitEuler(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMimplicitEuler.txt");
            break;
        case Msymmetric:
            result = symmetric(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMsymmetric.txt");
            break;
        case MrungeKutta2:
            result = rungeKutta2(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMrungeKutta2.txt");
            break;
        case MrungeKutta4:
            result = rungeKutta4(initPoints);
            outputMatrix(result, ADD_DOTS"data/outMrungeKutta4.txt");
            break;
        case MexplicitAdams:
            result = Adams(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMexplicitAdams.txt");
            break;
        case Mbdf2:
            result = bdf2(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMbdf2.txt");
            break;
        case Mbdf4:
            result = bdf4(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMbdf4.txt");
            break;
    }
    std::cout << "diff: " << getMaxDiff(calcDiff(result)) << "\n";
    assert(getMaxDiff(calcDiff(result)) < 0.01);
//    outputOnTheScreenMatrix(calcDiff(result));
}

void ImplicitEuler() {
    templateOutput(MexplicitEuler);
}

void ExplicitEuler() {
    templateOutput(MimplicitEuler);
}

void Symmetric() {
    templateOutput(Msymmetric);
}

void RungeKutta2() {
    templateOutput(MrungeKutta2);
}

void RungeKutta4() {
    templateOutput(MrungeKutta4);
}

void ExplicitAdams() {
    templateOutput(MexplicitAdams);
}

void BDF2() {
    templateOutput(Mbdf2);
}

void BDF4() {
    templateOutput(Mbdf4);
}
