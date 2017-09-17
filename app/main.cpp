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
    //                  (Z,     basis,               nElectrons,  maxIterations,  tollerance, outBasisFileName);
    Examples::SingleAtom(2,     "6-311G(2df,2pd)",    2,          1e4,            1e-10);
    return 0;
}

