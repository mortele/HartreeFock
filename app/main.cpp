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
#include "Integrators/numericalintegrator.h"

using std::cout;
using std::endl;

int main(int, char**) {
    System*         system = new System(1);

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
    NumericalIntegrator* integrator = new NumericalIntegrator(system);


    const arma::mat& C = solver->m_coefficientMatrix;
    const arma::mat& S = solver->m_overlapMatrix;
    cout << C << endl;
    cout << S << endl;
    cout << C.t() * S * C << endl;


    arma::mat D = arma::zeros<arma::mat>(2,2);

    D(0,0) = C(0,0);
    D(0,1) = C(0,0);
    D(1,0) = C(1,0);
    D(1,1) = C(1,0);

    D = 2*D*D.t();
    cout << D<< endl;
    cout << arma::trace(D) << endl;

    cout << integrator->integrateDensity(D)   << endl;

    //cout << integrator->testIntegral(solver->m_densityMatrix)       << endl;
    //cout << integrator->integrateDensity(solver->m_densityMatrix)   << endl;

    //cout << integrator->testIntegral(arma::eye(2,2)) << endl;
    //cout << integrator->integrateDensity(arma::eye(2,2)) << endl;
    //cout << I << endl;

    return 0;
}

