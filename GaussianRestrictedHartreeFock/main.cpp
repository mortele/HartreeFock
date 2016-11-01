//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <cassert>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Atoms/atom.h"
#include "Atoms/Hydrogen/hydrogen_321G.h"
#include "Atoms/Hydrogen/hydrogen_321gplus.h"
#include "Atoms/Hydrogen/hydrogen_631g.h"
#include "Atoms/Hydrogen/hydrogen_631gss.h"
#include "Parsers/basissetparser.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

int main(int, char**) {

    //BasisSetParser parser;
    //Atom* atom = parser.newAtomFromBasisSetFile("H", "3-21G", zeros<vec>(3));
    //cout << *atom << endl;
    //for (ContractedGaussian* contracted : atom->getContractedGaussians()) {
    //    for (GaussianPrimitive* primitive : *contracted) {
    //        cout << primitive << endl;
    //    }
    //}
    //return 0;

    vec nucleus1 {0, 0, 0};
    //vec nucleus2 {0, 0, 1.4};

    System system;
    Atom* Hminus = new Hydrogen_631Gss(nucleus1);
    Hminus->setNumberOfElectrons(2);
    system.addAtom(Hminus);

    //system.addAtom(new Hydrogen_631Gss(nucleus1));
    //system.addAtom(new Hydrogen_631Gss(nucleus2));

    RestrictedHartreeFock solver(&system);
    double result = solver.solve(1e-14, 1e3);

    //assert(std::fabs(-1.131284349300591 - result) < 1e-15);
    return 0;
}

