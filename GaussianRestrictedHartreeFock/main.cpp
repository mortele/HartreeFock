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

    vec nucleus {1,2,3};

    System singleHydrogenAtom(1);
    singleHydrogenAtom.addAtom(new Hydrogen(nucleus));
    cout << singleHydrogenAtom.overlapIntegral(0,0) << endl;
    cout << singleHydrogenAtom.overlapIntegral(0,1) << endl;
    cout << singleHydrogenAtom.overlapIntegral(1,0) << endl;



    return 0;
}

