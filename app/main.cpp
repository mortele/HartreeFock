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

    //system->addAtom(helium);
    //system->addAtom(hydrogen1);
    //system->addAtom(hydrogen2);
    system->addAtom(beryllium);
    //system->addAtom(carbon);

    RestrictedDFT*  solver = new RestrictedDFT(system);
    solver->solve(1e-10,100);
    NumericalIntegrator* integrator = new NumericalIntegrator(system);


    const arma::mat& C = solver->m_coefficientMatrix;
    const arma::mat& S = solver->m_overlapMatrix;
    const arma::mat& D = solver->m_densityMatrix;

    int n = system->getNumberOfSpinUpElectrons()*2;
    int b = system->getBasis().size();

    arma::mat C2 = arma::zeros<arma::mat>(b,n);
    arma::mat D2 = arma::zeros<arma::mat>(n,n);

    cout << "C:" << endl;
    cout << C << endl;
    for (int i = 0; i < n; i+=2) {
        int j = i/2;
        C2.col(i)   = C.col(j);
        C2.col(i+1) = C.col(j);
    }
    cout << "C2:" << endl;
    cout << C2 << endl;

    cout << "S:" << endl;
    cout << S << endl;

    cout << "D:" << endl;
    cout << D << endl;

    cout << "D2:" << endl;
    D2 = D*S;
    cout << D2 << endl;

    cout << "TRACE:" << endl;

    cout << arma::trace(D2) << endl;

    //cout << integrator->integrateDensity(D)   << endl;

    cout << integrator->testIntegral(D2)       << endl;
    cout << integrator->integrateDensity(D2)   << endl;

    //cout << integrator->testIntegral(arma::eye(2,2)) << endl;
    //cout << integrator->integrateDensity(arma::eye(2,2)) << endl;
    //cout << I << endl;

    return 0;
}

