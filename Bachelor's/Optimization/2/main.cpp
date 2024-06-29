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

bool all = true;

double func1(const double x, const double y) {
    //  func 2
    //int alpha = 128; //  8; 128
    //return alpha * pow((pow(x, 2) - y), 2) + pow(x - 1, 2);

    // func1
    return (5 * pow(x, 2) + 4 * x * y + 2 * pow(y, 2) + 4 * sqrt(5) * (x + y) + 51);
}

pair<double, double> Grad1(const double x, const double y) {
    //  8
    //return {-(2*(16 * pow(x, 3) - 16 * x * y + x - 1)), -(-16 * (pow(x, 2)  - y))};
    //  128
    //return {-(2 * (256 * x * (pow(x, 2) - y) + x - 1)), -(-256 * (pow(x, 2) - y))};

    return {(-4 * sqrt(5) - 10 * x - 4 * y), (-4 * sqrt(5) - 4 * x - 4 * y)};
}

double vectorNorm(const pair<double, double> p) {
    return sqrt(pow(p.first, 2) + pow(p.second, 2));
}

double GoldRate(double eps, const function<double(double)> &func, unsigned int &amtF) {
    double leftBound = 0;
    double rightBound = 5;
    double it = 1;
    amtF += 2;
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
        ++amtF;
        ++it;
    }
    return (rightBound + leftBound) / 2;
}


pair<double, double> FastDescent(const double initX, const double initY, const double eps,
                                 const function<pair<double, double>(double, double)> &func, const string &folder) {
    pair<double, double> currDir = func(initX, initY);
    double kappa = initKappa;
    pair<double, double> XY = {initX, initY};
    long int it = 0;
    unsigned int amtF = 0;
    unsigned int amtG = 1;
    vector<double> kappaVec{kappa};
    vector<pair<double, double>> XYvec{XY};
    vector<double> Fvec{func1(XY.first, XY.second)};
    vector<double> normVec{vectorNorm(currDir)};

    while (vectorNorm(currDir) > eps) {
//        if (XY.first == MIN_VAL || XY.second == MIN_VAL){
//            break;
//        }

        kappa = GoldRate(eps, [=](double
                                  kp) {
                             return func1(XY.first + kp * currDir.first, XY.second + kp * currDir.second);
                         }, amtF
        );
        XY.first += kappa * currDir.first;
        XY.second += kappa * currDir.second;
        currDir = func(XY.first, XY.second);
        ++amtG;

        ++it;
        if (all || it % 10 == 0 || it < 100) {
            kappaVec.push_back(kappa);
            XYvec.push_back(XY);
            Fvec.push_back(func1(XY.first, XY.second));
            normVec.push_back(vectorNorm(currDir));
        }
    }

    outputVector(kappaVec, folder + "/kappa.txt");
    outputPair(XYvec, folder + "/points.txt");
    outputVector(Fvec, folder + "/func.txt");
    outputVector(normVec, folder + "/norm.txt");

    cout << "amtFunc: " << amtF << "; amtGrad: " << amtG << "; ";
    return XY;
}

pair<double, double> GradDescent(const double initX, const double initY, const double eps,
                                 const function<pair<double, double>(double, double)> &func, const string &folder) {
    const double omega = 0.5;
    const double nu = 0.95;
    pair<double, double> currDir = func(initX, initY);
    double kappa = initKappa;
    vector<double> kappaVec{kappa};
    pair<double, double> XY = {initX, initY};
    vector<pair<double, double>> XYvec{XY};
    unsigned int amtF = 1;
    unsigned int amtG = 1;
    long int it = 0;
    double Fcurr = func1(XY.first, XY.second);
    vector<double> Fvec{Fcurr};
    vector<double> normVec{vectorNorm(currDir)};
    while (vectorNorm(currDir) > eps) {
//        if (XY.first <= MIN_VAL || XY.second <= MIN_VAL){
//            break;
//        }
        pair<double, double> XYcurr;
        XYcurr.first = XY.first + kappa * currDir.first;
        XYcurr.second = XY.second + kappa * currDir.second;
        while (Fcurr - func1(XYcurr.first, XYcurr.second) <= omega * kappa * pow(vectorNorm(currDir), 2)) {
            kappa *= nu;
            XYcurr.first = XY.first + kappa * currDir.first;
            XYcurr.second = XY.second + kappa * currDir.second;
            ++amtF;
        }
        XY = XYcurr;
        ++amtF;
        Fcurr = func1(XY.first, XY.second);
        ++amtG;
        currDir = func(XY.first, XY.second);

        ++it;
        if (all || it % 10 == 0) {
            kappaVec.push_back(kappa);
            XYvec.push_back(XY);
            Fvec.push_back(Fcurr);
            normVec.push_back(vectorNorm(currDir));
        }
    }

    outputVector(kappaVec, folder + "/kappa.txt");
    outputPair(XYvec, folder + "/points.txt");
    outputVector(Fvec, folder + "/func.txt");
    outputVector(normVec, folder + "/norm.txt");

    cout << "amtFunc: " << amtF << "; amtGrad: " << amtG << "; ";
    return XY;
}

int main() {
    const double X = -1.0;
    const double Y = -2.0;

    std::cout << "FastDescent:" << std::endl;
    pair<double, double> res = FastDescent(X, Y, 10e-3, Grad1, "../fast2");
    cout << "2 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second) << "\n";
    res = FastDescent(X, Y, 10e-7, Grad1, "../fast6");
    cout << "6 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second) << "\n";

    std::cout << "GradDescent:" << std::endl;
    res = GradDescent(X, Y, 10e-4, Grad1, "../grad2");
    cout << "2 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second) << "\n";
    res = GradDescent(X, Y, 10e-7, Grad1, "../grad6");
    cout << "6 - " << setprecision(5) << "x: " << res.first << "; y: "
         << res.second << "; f(x, y): " << func1(res.first, res.second) << "\n";
    return 0;
}
