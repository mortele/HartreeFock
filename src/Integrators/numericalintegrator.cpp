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
    m_grid->createSimpleOneAtomGrid(30,20,-1);
    cout << m_grid->getPoints(0,0) << endl;

    double integral = 0;
    for (int i = 0; i < 30; i++) {
        const double x = m_grid->getPoints(i,0);
        const double y = m_grid->getPoints(i,1);
        const double z = m_grid->getPoints(i,2);
        const double r = sqrt(x*x + y*y + z*z);
        integral += exp(-r)*m_grid->getWeights(i);
    }
    return integral;
}
