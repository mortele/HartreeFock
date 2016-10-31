//#define ARMA_MAT_PREALLOC 4
//#define ARMA_EXTRA_DEBUG

#include <iostream>
#include <iomanip>
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
    system.addAtom(new Hydrogen_321G(nucleus1));
    system.addAtom(new Hydrogen_321G(nucleus2));

    //for (Atom* atom : system.getAtoms()) {
    //    for (ContractedGaussian* contracted : atom->getContractedGaussians()) {
    //        for (GaussianPrimitive* primitive : contracted->getPrimitives()) {
    //            cout << *primitive << endl;
    //        }
    //    }
    //}

    for (uint i=0; i<system.getNumberOfBasisFunctions(); i++) {
        for (uint j=0; j<system.getNumberOfBasisFunctions(); j++) {
            for (uint k=0; k<system.getNumberOfBasisFunctions(); k++) {
                for (uint l=0; l<system.getNumberOfBasisFunctions(); l++) {
                    //cout << i << "," << j << "," << k << "," << l << ": ";
                    cout << 2*system.electronElectronIntegral(i,j,k,l)-system.electronElectronIntegral(i,j,l,k)  << endl;
                    //cout << system.electronElectronIntegral(i,j,k,l) << endl;
                }
            }
        }
    }

    RestrictedHartreeFock solver(&system);
    solver.solveSilently(1e-14, 2);
    return 0;
}

