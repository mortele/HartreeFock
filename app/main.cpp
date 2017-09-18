//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "examples.h"

#include "system.h"
#include "Solvers/restricteddft.h"
#include "Atoms/helium.h"
#include "Integrators/numericalintegrator.h"

using std::cout;
using std::endl;

int main(int, char**) {
    System*         system = new System(2);
    Helium*         helium = new Helium("3-21G", arma::vec{0,0,0});
    system->addAtom(helium);
    NumericalIntegrator* integrator = new NumericalIntegrator(system);
    RestrictedDFT*  solver = new RestrictedDFT(system);
    solver->solveSilently(1e-10,100);
    //const double I = integrator->testIntegral();
    const double I = integrator->integrateDensity(solver->m_densityMatrix);
    cout << I << endl;

    return 0;
}

