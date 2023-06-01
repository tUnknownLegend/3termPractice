#include <iostream>
#include "methods.h"

int main() {
//    std::cout << "ImplicitEuler:\n";
//    ImplicitEuler();
//    std::cout << "\nExplicitEuler:\n";
//    ExplicitEuler();
//    std::cout << "\nSymmetric:\n";
//    Symmetric();
    std::cout << "\nRungeKutta, 2:\n";
    RungeKutta2();
    std::cout << "\nRungeKutta, 4:\n";
    RungeKutta4();
//    std::cout << "\nexplicit Adams, 4:\n";
//    ExplicitAdams();

    return 0;
}
