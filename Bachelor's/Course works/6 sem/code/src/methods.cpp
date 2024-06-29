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
    MimplicitAdams,
    Mimplicit2Adams,
    Mbdf2,
    Mbdf4
};

void saveAxesData(const vector<vector<TT>> &answer, const unsigned int axe, const string &fileName) {
    vector<vector<TT>> values;
    values.reserve(answer.size());
    for (int i = 0; i < answer.size(); ++i) {
        values.push_back({i * step, answer[i][axe]});
    }
    outputMatrix(values, fileName);
}

void saveAllAxes(const vector<vector<TT>> &answer, const string &file1, const string &file2, const string &file3) {
    saveAxesData(answer, 0, file1);
    saveAxesData(answer, 1, file2);
    saveAxesData(answer, 2, file3);
}

vector<vector<TT>> calcDiff(const vector<vector<TT>> &answer) {
    vector<vector<TT>> diff(answer);

    for (int i = 0; i < answer.size(); ++i) {
        diff[i][0] = fabs(cos(i * step) - answer[i][0]);
        diff[i][1] = fabs(-sin(i * step) - answer[i][1]);
//        diff[i][0] = fabs(-exp(2 * i) + i * exp(i) + 2 * exp(i) - answer[i][0]);
//        diff[i][1] = fabs(i * exp(i) - answer[i][1]);
//        diff[i][2] = fabs(-exp(2 * i) + i * exp(i) + exp(i) - answer[i][1]);
    }
    return diff;
}

TT getP(const vector<vector<TT>> &answer) {
    const unsigned int i = 10 * multiplayer;
    return max(fabs(cos(i * step) - answer[i][0]), fabs(-sin(i * step) - answer[i][1]));
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

vector<vector<TT>> rungeKutta4(const vector<TT> &cond, const int n) {
    vector<TT> previous(cond.size());
    vector<TT> help(cond.size());
    vector<vector<TT>> K(4);

    vector<vector<TT>> x;
    x.emplace_back();
    for (const TT i: cond) {
        x[0].push_back(i);
    }

    for (int i = 0; i < 4; ++i) { K[i] = help; };
    for (int t = 0; t < n; ++t) {
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
        const auto firstYs = rungeKutta4(cond, n);
        for (int i = 0; i < 4; ++i) {
            y[i] = firstYs[i];
        }
    }

    for (int i = 3; i < n - 1; ++i) {
        y[i + 1] = Newton("BDF4", cond.size(), {y[i - 3], y[i - 2], y[i - 1], y[i]});
    }
    return y;
}

vector<vector<TT>> implicit2Adams(const vector<TT> &cond, int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    y[0] = rungeKutta2(cond, numOfPoints)[0];
    y[1] = rungeKutta2(cond, numOfPoints)[1];

    for (int i = 1; i < n - 1; ++i) {
        y[i + 1] = Newton("implicit2Adams", cond.size(), {y[i - 1], y[i]});
    }
    return y;
}

vector<vector<TT>> implicitAdams(const vector<TT> &cond, int n) {
    vector<vector<TT>> y(n, vector<TT>(cond.size()));
    {
        const auto firstYs = rungeKutta4(cond, n);
        for (int i = 0; i < 4; ++i) {
            y[i] = firstYs[i];
        }
    }

    for (int i = 3; i < n - 1; ++i) {
        y[i + 1] = Newton("implicitAdams", cond.size(), {y[i - 3], y[i - 2], y[i - 1], y[i]});
    }
    return y;
}

vector<vector<TT>> Adams(const vector<TT> &cond, const int n) {
    vector<vector<TT>> x(n, vector<TT>(cond.size()));
    vector<vector<TT>> y4 = rungeKutta4(cond, n);
    for (int i = 0; i < 4; ++i) {
        x[i] = y4[i];
    }

    for (int t = 3; t < n - 1; ++t) {
        const vector<TT> currF = f(x[t]);
        const vector<TT> F1 = f(x[t - 1]);
        const vector<TT> F2 = f(x[t - 2]);
        const vector<TT> F3 = f(x[t - 3]);

        for (int j = 0; j < cond.size(); ++j) {
            x[t + 1][j] = x[t][j] + step / 24 * (55 * currF[j] - 59 * F1[j]
                                                 + 37 * F2[j] - 9 * F3[j]);
        }
    }
    return x;
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
            result = rungeKutta4(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMrungeKutta4.txt");
            break;
        case MexplicitAdams:
            result = Adams(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMexplicitAdams.txt");
            break;
        case MimplicitAdams:
            result = implicitAdams(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMimplicitAdams.txt");
            break;
        case Mimplicit2Adams:
            result = implicit2Adams(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMimplicit2Adams.txt");
            saveAllAxes(result, ADD_DOTS"data/y1/outMimplicit2Adams.txt",
                        ADD_DOTS"data/y2/outMimplicit2Adams.txt",
                        ADD_DOTS"data/y3/outMimplicit2Adams.txt");
            break;
        case Mbdf2:
            result = bdf2(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMbdf2.txt");
            saveAllAxes(result, ADD_DOTS"data/y1/outMbdf2.txt",
                        ADD_DOTS"data/y2/outMbdf2.txt",
                        ADD_DOTS"data/y3/outMbdf2.txt");
            break;
        case Mbdf4:
            result = bdf4(initPoints, numOfPoints);
            outputMatrix(result, ADD_DOTS"data/outMbdf4.txt");
            saveAllAxes(result, ADD_DOTS"data/y1/outMbdf4.txt",
                        ADD_DOTS"data/y2/outMbdf4.txt",
                        ADD_DOTS"data/y3/outMbdf4.txt");
            break;
    }
//    std::cout << "diff: " << getMaxDiff(calcDiff(result)) << "\n";
//    assert(getMaxDiff(calcDiff(result)) < 0.01);
//    std::cout << "p: " << getP(result) << "\n";
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

void ImplicitAdams() {
    templateOutput(MimplicitAdams);
}

void Implicit2Adams() {
    templateOutput(Mimplicit2Adams);
}

void BDF2() {
    templateOutput(Mbdf2);
}

void BDF4() {
    templateOutput(Mbdf4);
}
