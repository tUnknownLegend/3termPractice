#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <fstream>
#include <limits>

#define initKappa 1e0
#define MIN_VAL std::numeric_limits<double>::max();

using namespace std;

void outputVector(const vector<double> &vect, const string &fileName) {
    ofstream outFile(fileName);
    if (!outFile.is_open()) {
        cerr << "error // output.txt open\n";
        return;
    }
    outFile << vect.size() << std::endl;
    {
        double node = 0.0;
        for (auto &el: vect) {
            outFile << el << endl;
        }
        outFile << std::endl;
    }
    outFile.close();
}

void outputPair(const vector<pair<double, double>> &vect, const string &fileName) {
    ofstream outFile(fileName);
    if (!outFile.is_open()) {
        cerr << "error // output.txt open\n";
        return;
    }
    outFile << vect.size() << std::endl;
    {
        double node = 0.0;
        for (auto &el: vect) {
            outFile << el.first << " " << el.second << endl;
        }
        outFile << std::endl;
    }
    outFile.close();
}

double scalMultiply(const pair<double, double> a, const pair<double, double> b) {
    return (a.first * b.first + a.second * b.second);
}

pair<double, double> Mminus(const pair<double, double> a, const pair<double, double> b) {
    return {a.first - b.first, a.second - b.second};
}

pair<double, double> Mplus(const pair<double, double> a, const pair<double, double> b) {
    return {a.first + b.first, a.second + b.second};
}

pair<double, double> Mmultiply(double a, const pair<double, double> b) {
    return {a * b.first, a * b.second};
}

bool all = true;

double func1(const double x, const double y, const char option) {
    switch (option) {
        case 1:  //  quadratic
            return ((5 * pow(x, 2) + 4 * x * y + 2 * pow(y, 2) + 4 * sqrt(5) * (x + y) + 51));
        case 2:  //  8
            return (8 * pow((pow(x, 2) - y), 2) + pow(x - 1, 2));
        case 3:  //  128
            return (128 * pow((pow(x, 2) - y), 2) + pow(x - 1, 2));
        default:
            cout << "wrong option!\n";
            return 0;
    }
}

pair<double, double> Grad1(const double x, const double y, const char option) {
    switch (option) {
        case 1:  //  quadratic
            return {(-4 * sqrt(5) - 10 * x - 4 * y), (-4 * sqrt(5) - 4 * x - 4 * y)};
        case 2:  //  8
            return {-(2 * (16 * pow(x, 3) - 16 * x * y + x - 1)), -(-16 * (pow(x, 2) - y))};
        case 3:  //  128
            return {-(2 * (256 * x * (pow(x, 2) - y) + x - 1)), -(-256 * (pow(x, 2) - y))};
        default:
            cout << "wrong option!\n";
            return {0, 0};
    }
}

double Guess1(const double x, const double y, const char option) {
    switch (option) {
        case 1:  //  quadratic
            return 24;
        case 2:  //  8
            return (32 + 512 * (pow(x, 2) - y));
        case 3:  //  128
            return (512 + 131072 * (pow(x, 2) - y));
        default:
            cout << "wrong option!\n";
            return 0;
    }
}

double vectorNorm(const pair<double, double> p) {
    return sqrt(pow(p.first, 2) + pow(p.second, 2));
}

double GoldRate(double eps, const function<double(double)> &func, unsigned int &amtFunc) {
    double leftBound = 0;
    double rightBound = 5;
    double it = 1;

    double nextLeft = leftBound + (1 - 1 / ((1 + sqrt(5)) / 2)) * (leftBound + rightBound);
    double nextRight = leftBound + (leftBound + rightBound) * (1 / ((1 + sqrt(5)) / 2));
    double leftFunc = func(nextLeft);
    double rightFunc = func(nextRight);
    amtFunc += 2;

    while ((rightBound - leftBound) > eps) {

        if (leftFunc > rightFunc) {
            leftBound = nextLeft;
            nextLeft = nextRight;
            leftFunc = rightFunc;
            nextRight = leftBound + rightBound - nextLeft;
            rightFunc = func(nextRight);
        } else {
            rightBound = nextRight;
            nextRight = nextLeft;
            rightFunc = leftFunc;
            nextLeft = leftBound + rightBound - nextRight;
            leftFunc = func(nextLeft);
        }
        ++amtFunc;
        ++it;
    }
    return (rightBound + leftBound) / 2;
}

pair<double, double> Descent(const double initX, const double initY, const double eps, const string &folder,
                             const char option, const char method) {
    pair<double, double> currDir = Grad1(initX, initY, option);
    double kappa = initKappa;
    pair<double, double> XY = {initX, initY};
    long int it = 0;
    vector<double> kappaVec{kappa};
    vector<pair<double, double>> XYvec{XY};
    vector<double> Fvec{func1(XY.first, XY.second, option)};
    vector<double> normVec{vectorNorm(currDir)};
    unsigned int amtFunc = 0;
    unsigned int amtGrad = 1;


    int amtDots = 0;
    double gamma = 0;
    pair<double, double> p = currDir;
    while (vectorNorm(currDir) > eps) {
        ++amtDots;
        kappa = GoldRate(eps, [=](double
                                  kp) {
                             return func1(XY.first + kp * currDir.first, XY.second + kp * currDir.second, option);
                         }, amtFunc
        );

        double h = Guess1(XY.first, XY.second, option);
        XY.first += kappa * currDir.first;
        XY.second += kappa * currDir.second;
        pair<double, double> nextDir = Grad1(XY.first, XY.second, option);
        ++amtGrad;

        if (all || abs(amtDots % 3) <= 10e-10) {
            switch (method) {
                case 1:
                    gamma = pow((pow(vectorNorm(nextDir), 2) / (vectorNorm(currDir))), 2);
                    break;
                case 2:
                    gamma = scalMultiply(Mminus(nextDir, currDir), nextDir) / pow(vectorNorm(currDir), 2);
                    break;
                case 3:
                    gamma = -scalMultiply(Mmultiply(h, p), nextDir) / scalMultiply(Mmultiply(h, p), p);
                    break;
                default:
                    cout << "wrong option!";
                    return {0, 0};
            }
        } else {
            gamma = 0;
        }
        currDir = nextDir;
        p = Mplus(Mmultiply(gamma, p), currDir);
        currDir = p;
        ++it;
        if (all || it % 10 == 0 || it < 100) {
            kappaVec.push_back(kappa);
            XYvec.push_back(XY);
            Fvec.push_back(func1(XY.first, XY.second, option));
            normVec.push_back(vectorNorm(currDir));
        }
    }

    outputVector(kappaVec, folder + "/kappa.txt");
    outputPair(XYvec, folder + "/points.txt");
    outputVector(Fvec, folder + "/func.txt");
    outputVector(normVec, folder + "/norm.txt");

    cout << "amtFunc: " << amtFunc << "; amtGrad: " << amtGrad;
    return XY;
}

int main() {
    const double X = -1.0;
    const double Y = -2.0;
    const char option = 3;
    const char method = 1;

    std::cout << "Fletcher:" << std::endl;
    pair<double, double> res = Descent(X, Y, 10e-3, "../fast2", option, method);
    cout << "2 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second, option) << "\n";
    res = Descent(X, Y, 10e-7, "../fast6", option, method);
    cout << "6 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second, option) << "\n";

    std::cout << "Polak:" << std::endl;
    res = Descent(X, Y, 10e-3, "../grad2", option, 2);
    cout << "2 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second, option) << "\n";
    res = Descent(X, Y, 10e-7, "../grad6", option, 2);
    cout << "6 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second, option) << "\n";

    std::cout << "H:" << std::endl;
    res = Descent(X, Y, 10e-3, "../mod2", option, 3);
    cout << "2 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second, option) << "\n";
    res = Descent(X, Y, 10e-7, "../mod6", option, 3);
    cout << "6 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second, option) << "\n";
    return 0;
}
