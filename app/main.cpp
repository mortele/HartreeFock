//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "examples.h"

#include "system.h"
#include "Solvers/restricteddft.h"
#include "Solvers/restrictedhartreefock.h"
#include "Atoms/hydrogen.h"
#include "Atoms/helium.h"
#include "Atoms/beryllium.h"
#include "Atoms/carbon.h"

using std::cout;
using std::endl;
using std::setprecision;


int main(int, char**) {

    System*         system = new System();

    //Helium*         helium      = new Helium    ("3-21G", arma::vec{0,0,0});
    //Helium*         helium      = new Helium    ("6-311+G**", arma::vec{0,0,0});
    //Helium*         helium      = new Helium    ("6-311G(2df,2pd)", arma::vec{0,0,0});
    Helium*         helium      = new Helium    ("toy", arma::vec{0,0,0});
    Hydrogen*       hydrogen1   = new Hydrogen  ("3-21G", arma::vec{0,0,0});
    Hydrogen*       hydrogen2   = new Hydrogen  ("3-21G", arma::vec{0,0,1.4});
    Beryllium*      beryllium   = new Beryllium ("3-21G", arma::vec{0,0,0});
    //Beryllium*      beryllium   = new Beryllium ("6-311++G**", arma::vec{0,0,0});
    Carbon*         carbon      = new Carbon    ("3-21G", arma::vec{0,0,0});

    //system->addAtom(hydrogen1);
    //system->addAtom(hydrogen2);
    system->addAtom(helium);
    //system->addAtom(beryllium);
    //system->addAtom(carbon);

    RestrictedHartreeFock*  rhf     = new RestrictedHartreeFock(system);
    RestrictedDFT*          rdft    = new RestrictedDFT(system);
    rdft->setFunctional("LDA");
    rdft->solve(1e-8,20);

    cout << rdft->m_coefficientMatrix << endl;
    rdft->m_coefficientMatrix(0,0) = 0.295500;
    rdft->m_coefficientMatrix(1,0) = 0.815618;

    //rhf->solve();

    return 0;
}

