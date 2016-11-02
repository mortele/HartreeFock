//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <cassert>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Solvers/unrestrictedhartreefock.h"
#include "Atoms/atom.h"
#include "Atoms/hydrogen.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {

    vec nucleus1 {0, -0.2, 0};
    vec nucleus2 {0, 0, 3.321};

    System* system = new System(1);
    system->addAtom(new Hydrogen("6-311++G**", nucleus1));
    system->addAtom(new Hydrogen("6-31G**", nucleus2));

    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);
    double result = solver.solve(1e-14, 1e4);

    //assert(std::fabs(-1.131284349300591 - result) < 1e-15);
    return 0;
}

// -13.25093 eV  UHF H-
// -14.348   eV      H- (http://nist.gov/data/PDFfiles/jpcrd68.pdf)
// -13.6     eV      H
