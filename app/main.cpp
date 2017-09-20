//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "examples.h"

#include "system.h"
#include "Solvers/restricteddft.h"
#include "Atoms/hydrogen.h"
#include "Atoms/helium.h"
#include "Atoms/beryllium.h"
#include "Atoms/carbon.h"

using std::cout;
using std::endl;
using std::setprecision;

int main(int, char**) {
    System*         system = new System();

    Helium*         helium      = new Helium    ("3-21G", arma::vec{0,0,0});
    Hydrogen*       hydrogen1   = new Hydrogen  ("3-21G", arma::vec{0,0,0});
    Hydrogen*       hydrogen2   = new Hydrogen  ("3-21G", arma::vec{0,0,1.4});
    Beryllium*      beryllium   = new Beryllium ("3-21G", arma::vec{0,0,0});
    Carbon*         carbon      = new Carbon    ("3-21G", arma::vec{0,0,0});

    system->addAtom(helium);
    //system->addAtom(hydrogen1);
    //system->addAtom(hydrogen2);
    //system->addAtom(beryllium);
    //system->addAtom(carbon);

    RestrictedDFT*  solver = new RestrictedDFT(system);
    solver->solve(1e-10,100);
    cout << setprecision(15) << (3.0/4.0)*pow((3.0/M_PI),1.0/3.0) << endl;
    return 0;
}

