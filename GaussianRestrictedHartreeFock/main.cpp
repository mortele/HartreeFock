//#define ARMA_MAT_PREALLOC 4
//#define ARMA_EXTRA_DEBUG

#include <iostream>
#include <iomanip>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Atoms/Hydrogen/hydrogen_321G.h"
#include "Atoms/Hydrogen/hydrogen_321gplus.h"
#include "Atoms/Hydrogen/hydrogen_631g.h"

using std::cout;
using std::endl;
using std::setprecision;
using arma::vec;
using arma::zeros;


int main(int, char**) {

    vec nucleus1 {0, 0, 0};
    vec nucleus2 {0, 0, 1};

    System system(2);
    system.addAtom(new Hydrogen_631G(nucleus1));
    system.addAtom(new Hydrogen_631G(nucleus2));

    RestrictedHartreeFock solver(&system);

    return 0;
}

