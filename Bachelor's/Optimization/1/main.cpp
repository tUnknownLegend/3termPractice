#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double func(double x) {
    return exp((pow(x, 4) + 2 * pow(x, 3) - 5 * x + 6)) +
           cosh(1 / (-15 * pow(x, 3) + 10 * x + 5 * sqrt(10) - 3));
}

double Dexot(double eps) {
    double delta = eps / 5;
    double leftBound = 0;
    double rightBound = 1;
    double it = 0;

    while ((rightBound - leftBound) > eps) {
        double nextLeft = (leftBound + rightBound) / 2 - delta;
        double nextRight = (leftBound + rightBound) / 2 + delta;
        double leftFunc = func(nextLeft);
        double rightFunc = func(nextRight);

        if (leftFunc <= rightFunc) {
            rightBound = nextRight;
        } else {
            leftBound = nextLeft;
        }
        ++it;
    }
    return func((rightBound + leftBound) / 2);
}

double GoldRate(double eps) {
    double leftBound = 0;
    double rightBound = 1;
    double it = 1;

    double nextLeft = leftBound + (1 - 1 / ((1 + sqrt(5)) / 2)) * (leftBound + rightBound);
    double nextRight = leftBound + (leftBound + rightBound) * (1 / ((1 + sqrt(5)) / 2));
    double leftFunc = func(nextLeft);
    double rightFunc = func(nextRight);


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
        ++it;
    }
    return func((rightBound + leftBound) / 2);
}

int main() {
    std::cout << "Dexot:" << std::endl;
    cout << "2: " << setprecision(5) << Dexot(10e-2) << "\n";
    cout << "6: "<< setprecision(9) << Dexot(10e-6) << "\n";
    cout << "17: "<< setprecision(20) << Dexot(10e-17) << "\n";

    std::cout << "Golden:" << std::endl;
    cout << "2: "<< setprecision(5) << GoldRate(10e-2) << "\n";
    cout << "6: "<< setprecision(9) << GoldRate(10e-6) << "\n";
    cout << "17: "<< setprecision(20) << GoldRate(10e-17) << "\n";
    return 0;
}
