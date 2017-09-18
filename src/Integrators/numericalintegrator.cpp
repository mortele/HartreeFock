#include "numericalintegrator.h"
#include "system.h"
#include "Integrators/grid.h"
#include <cmath>

using arma::vec;
using arma::mat;
using std::sqrt;
using std::exp;
using std::cout;
using std::endl;


NumericalIntegrator::NumericalIntegrator(System* system) {
    m_system = system;
    m_grid   = new Grid(system);
}

double NumericalIntegrator::testIntegral() {
    m_grid->createSimpleOneAtomGrid(100,150,5);

    const vec& w = m_grid->getWeights();
    const mat& p = m_grid->getPoints();

    double integral = 0;
    for (int i = 0; i < w.n_elem; i++) {
        const double x = p(i,0);
        const double y = p(i,1);
        const double z = p(i,2);
        const double r = sqrt(x*x + y*y + z*z);
        integral += exp(-r)*w(i);
    }
    return integral;
}
