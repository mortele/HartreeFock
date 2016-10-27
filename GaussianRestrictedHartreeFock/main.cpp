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

    vec nucleus1 {0,0,0};
    vec nucleus2 {0,0,1.4};

    System singleHydrogenAtom(2);
    singleHydrogenAtom.addAtom(new Hydrogen(nucleus1));
    singleHydrogenAtom.addAtom(new Hydrogen(nucleus2));

    for (int i=0; i<singleHydrogenAtom.getNumberOfBasisFunctions(); i++) {
        for (int j=0; j<singleHydrogenAtom.getNumberOfBasisFunctions(); j++) {

            cout << "i,j: " << i << "," << j << "    " << singleHydrogenAtom.overlapIntegral(i,j) << endl;
        }
    }

    return 0;
}

