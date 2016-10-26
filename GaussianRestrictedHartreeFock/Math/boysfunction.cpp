#include "Math/boysfunction.h"
#include <cmath>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>

using boost::math::tgamma;
using boost::math::tgamma_lower;
using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;


double BoysFunction::directIntegration(double x, double n) {
    const double dt = 1.0 / m_numberOfIntegrationPoints;
    double t        = 0.0;
    double integral   = 0.0;
    for (int i=0; i < m_numberOfIntegrationPoints; i++) {
        integral    += integrand(x,t,n);
        t           += dt;
    }
    return integral*dt;
}


double BoysFunction::integrand(double x, double t, double n) {
    double tPower2n = 1.0;
    for (int i = 0; i < 2*n; i++) {
        tPower2n *= t;
    }
    return exp(-x*t*t)*tPower2n;
}


double BoysFunction::analyticalIncompleteGammaFunction(double x, double n) {
    /* boost::math::tgamma_lower(x,n) is the lower incomplete gamma function
     *
     *     x
     *   ╭    n-1  -t
     *   ┃  t     e   dt  =  gamma   (x, n)
     *   ╯                         IL
     *     0
     */
    const double nPlusOneHalf = n+0.5;
    return 1.0/(2*pow(x,nPlusOneHalf)) * boost::math::tgamma_lower(nPlusOneHalf,x);
}


double BoysFunction::analyticalCompleteGammaFunction(double x, double n) {
    /* boost::math::tgamma(x,n) is the upper incomplete gamma function.
     * boost::math::tgamma(n) is the "true" gamma function (not normalized).
     *
     *     oo
     *   ╭    n-1  -t
     *   ┃  t     e   dt  =  gamma   (x, n)
     *   ╯                         UL
     *     x
     *
     *     oo
     *   ╭    n-1  -t
     *   ┃  t     e   dt  =  gamma (n)
     *   ╯
     *     0
     *
     *
     */
    double nPlusOneHalf = n+0.5;
    return 1.0/(2*pow(x,nPlusOneHalf)) * (boost::math::tgamma(nPlusOneHalf) -
                                          boost::math::tgamma(nPlusOneHalf,x));
}




double BoysFunction::compute(double x, double n) {
    return analyticalIncompleteGammaFunction(x,n);
}

double BoysFunction::computeAndApplyDownwardRecurrence(double x, double n) {
    m_recurrenceValues = zeros<vec>(n+1);
    m_recurrenceValues(n) = compute(x, n);

    const double expMinusX = std::exp(-x);
    for (int m=n; m>0; m--) {
        m_recurrenceValues(m-1) = (2*x*m_recurrenceValues(m) + expMinusX) / (2.0*m-1.0);
    }
    return m_recurrenceValues(0);
}





