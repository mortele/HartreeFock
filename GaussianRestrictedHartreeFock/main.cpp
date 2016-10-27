//#define ARMA_MAT_PREALLOC 4
//#define ARMA_EXTRA_DEBUG

#include <iostream>
#include <iomanip>
#include "system.h"
#include "Atoms/hydrogen.h"

using std::cout;
using std::endl;
using std::setprecision;
using arma::vec;
using arma::zeros;


int main(int, char**) {

    vec nucleus1 {0, 0, 0};
    vec nucleus2 {0, 0, 1};

    System twoHydrogenAtom(2);
    twoHydrogenAtom.addAtom(new Hydrogen(nucleus1));
    twoHydrogenAtom.addAtom(new Hydrogen(nucleus2));

    for (int i=0; i<twoHydrogenAtom.getNumberOfBasisFunctions(); i++) {
        for (int j=0; j<twoHydrogenAtom.getNumberOfBasisFunctions(); j++) {
            for (int k=0; k<twoHydrogenAtom.getNumberOfBasisFunctions(); k++) {
                for (int l=0; l<twoHydrogenAtom.getNumberOfBasisFunctions(); l++) {

                    cout << "i,j,k,l: " << i << "," << j << "," << k << ","  << l << "    "
                         << twoHydrogenAtom.electronElectronIntegral(i,j,k,l) << endl;
                }
            }
        }
    }

    return 0;
}

