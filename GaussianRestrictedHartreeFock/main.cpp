//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <cassert>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Atoms/Hydrogen/hydrogen_321G.h"
#include "Atoms/Hydrogen/hydrogen_321gplus.h"
#include "Atoms/Hydrogen/hydrogen_631g.h"
#include "Atoms/Hydrogen/hydrogen_631gss.h"

using arma::vec;
using std::cout;
using std::endl;

int main(int, char**) {
    vec nucleus1 {0, 0, 0};
    vec nucleus2 {0, 0, 1.4};

    System system;
    system.addAtom(new Hydrogen_631Gss(nucleus1));
    system.addAtom(new Hydrogen_631Gss(nucleus2));

    RestrictedHartreeFock solver(&system);
    double result = solver.solve(1e-14, 1e3);

    assert(std::fabs(-1.131284349300591 - result) < 1e-15);
    return 0;
}

