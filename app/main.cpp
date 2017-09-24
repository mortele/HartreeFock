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
#include "Orbitals/gaussianprimitive.h"
#include "Orbitals/contractedgaussian.h"
#include "Atoms/hydrogen.h"
#include "Atoms/helium.h"
#include "Atoms/beryllium.h"
#include "Atoms/carbon.h"
#include "Integrators/numericalintegrator.h"

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
    rdft->solve(1e-8,0);

    // HF
    rdft->m_coefficientMatrix(0,0) = 0.300859;
    rdft->m_coefficientMatrix(1,0) = 0.811650;

    // DFT
    rdft->m_coefficientMatrix(0,0) = 0.295500;
    rdft->m_coefficientMatrix(1,0) = 0.815618;

    rdft->m_smoothing = false;
    rdft->computeDensityMatrix();
    arma::mat& S = rdft->m_overlapMatrix;
    arma::mat& P = rdft->m_densityMatrix;
    arma::mat& C = rdft->m_coefficientMatrix;
    rdft->m_densityMatrix = 2*C*C.t();
    GaussianPrimitive* pr1 = system->getBasis().at(0)->getPrimitives().at(0);
    GaussianPrimitive* pr2 = system->getBasis().at(1)->getPrimitives().at(0);
    //cout << *pr1 << endl;
    //cout << *pr2 << endl;

    cout << setprecision(10) << rdft->m_numericalIntegrator->testIntegral() << endl;

    //rhf->solve();

    return 0;
}

