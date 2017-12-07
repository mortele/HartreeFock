//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "examples.h"

using std::cout;
using std::endl;

int main(int, char**) {
    //Examples::firstExample();
    //Examples::secondExample();
    //Examples::methane();
    //Examples::diberyllium();
    //Examples::H2();

    //                  (Z,     basis,               nElectrons,  maxIterations,  tollerance, outBasisFileName);
    Examples::SingleAtom(10,     "3-21G",        10,           1e4,            1e-10);//,      "He-STO-6G-HF");
    return 0;
}

