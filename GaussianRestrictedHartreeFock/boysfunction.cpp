#include "boysfunction.h"
#include <cmath>

double BoysFunction::directIntegration(double x, double n) {
    const double dt = 1.0 / m_numberOfIntegrationPoints;
    double t        = 0.0;
    double integral   = 0.0;
    for (int i=0; i < m_numberOfIntegrationPoints; i++) {
        integral    += integrand(x,t,n);
        t           += dt;
    }
    return integral;
}

double BoysFunction::integrand(double x, double t, double n) {
    double tPower2n = 1.0;
    for (int i = 0; i < 2*n; i++) {
        tPower2n *= t;
    }
    return exp(-x*t*t)*tPower2n;
}

BoysFunction::BoysFunction() {

}
