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
#include "Integrators/numericalintegrator.h"

using std::cout;
using std::endl;

int main(int, char**) {
    System*         system = new System(0);
    Helium*         helium = new Helium("3-21G", arma::vec{2,0,0});
    Hydrogen*       hydrogen1 = new Hydrogen("3-21G", arma::vec{1,0,0});
    Hydrogen*       hydrogen2 = new Hydrogen("3-21G", arma::vec{0,0,1.4});
    Beryllium*      beryllium = new Beryllium("3-21G", arma::vec{0,0,0});

    //system->addAtom(helium);

    //system->addAtom(hydrogen1);
    //system->addAtom(hydrogen2);

    system->addAtom(beryllium);

    RestrictedDFT*  solver = new RestrictedDFT(system);
    solver->solve(1e-10,100);
    NumericalIntegrator* integrator = new NumericalIntegrator(system);

    cout << solver->m_densityMatrix << endl;
    const double I = integrator->testIntegral(solver->m_densityMatrix);
    //const double I = integrator->integrateDensity(solver->m_densityMatrix);
    cout << I << endl;

    return 0;
}

