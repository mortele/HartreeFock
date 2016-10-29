//#define ARMA_MAT_PREALLOC 4
//#define ARMA_EXTRA_DEBUG

#include <iostream>
#include <iomanip>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Atoms/Hydrogen/hydrogen_321G.h"
#include "Atoms/Hydrogen/hydrogen_321gplus.h"
#include "Atoms/Hydrogen/hydrogen_631g.h"
#include "Atoms/Hydrogen/hydrogen_31gss.h"

using arma::vec;

int main(int, char**) {
    vec nucleus1 {0, 0, 0};
    vec nucleus2 {0, 0, 1.4};

    System system;
    system.addAtom(new Hydrogen_31Gss(nucleus1));
    system.addAtom(new Hydrogen_31Gss(nucleus2));

    RestrictedHartreeFock solver(&system);
    solver.solve();
    return 0;
}

