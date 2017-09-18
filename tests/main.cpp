#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <armadillo>
#include <cmath>

#include "system.h"
#include "Atoms/helium.h"
#include "Integrators/grid.h"
#include "Integrators/numericalintegrator.h"

using namespace arma;
using namespace std;


TEST_CASE("Simple grid integral", "[Numerical integration]") {
    System*         system = new System(2);
    Helium*         helium = new Helium("3-21G", arma::vec{0,0,0});
    system->addAtom(helium);

    NumericalIntegrator* integrator  = new NumericalIntegrator(system);
    const double rMax   = 5;
    const double pi     = acos(-1.0);
    const double I      = integrator->testIntegral(arma::ones<mat>(2,2));
    cout << "I " <<  I << endl;
    REQUIRE( I == Approx(4*pi*(2-(rMax*(rMax+2)+2)*exp(-rMax))) );
}
