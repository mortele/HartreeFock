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
    const double I = integrator->testIntegral();
    return 1;
    RestrictedDFT*  solver = new RestrictedDFT(system);
    solver->solve(1e-10,100);

    return 0;
}

