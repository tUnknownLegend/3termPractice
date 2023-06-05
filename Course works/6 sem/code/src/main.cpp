#include <iostream>
#include "methods.h"

int main() {
//    std::cout << "ImplicitEuler:\n";
//    ImplicitEuler();
//    std::cout << "\nExplicitEuler:\n";
//    ExplicitEuler();
//    std::cout << "\nSymmetric:\n";
//    Symmetric();
//    std::cout << "\nRungeKutta, 2:\n";
//    RungeKutta2();
//    std::cout << "\nRungeKutta, 4:\n";
//    RungeKutta4();
//    std::cout << "\nexplicit Adams, 4:\n";
//    ExplicitAdams();
//    std::cout << "\nimplicit Adams, 4:\n";
//    ImplicitAdams();

    std::cout << "\nimplicit Adams, 2:\n";
    Implicit2Adams();
    std::cout << "\nbdf, 2:\n";
    BDF2();
    std::cout << "\nbdf, 4:\n";
    BDF4();

    return 0;
}
