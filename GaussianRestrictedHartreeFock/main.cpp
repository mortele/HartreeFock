#include <iostream>
#include <iomanip>
#include <time.h>
#include <boost/math/special_functions/gamma.hpp>
#include "Orbitals/gaussianprimitive.h"
#include "Integrators/overlapintegrator.h"
#include "Integrators/kineticintegrator.h"
#include "Factorizations/hermitegaussian.h"
#include "Math/boysfunction.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::pow;
using arma::vec;
using arma::zeros;


double boysAnalyticalIncomplete(double x, double m) {
    double mPlusOneHalf = m+0.5;
    return 1.0/(2*pow(x,mPlusOneHalf)) * boost::math::tgamma_lower(mPlusOneHalf,x);
}

double boysAnalyticalComplete(double x, double m) {
    double mPlusOneHalf = m+0.5;
    return 1.0/(2*pow(x,mPlusOneHalf)) * (boost::math::tgamma(mPlusOneHalf) -
                                          boost::math::tgamma(mPlusOneHalf,x));
}

int main(int, char**) {

    vec A = {1,0,0};
    vec B = {0,0,1};

    int x1 = 0;
    int y1 = 1;
    int z1 = 0;

    int x2 = 2;
    int y2 = 1;
    int z2 = 0;

    GaussianPrimitive primitive1 = GaussianPrimitive(x1,y1,z1, 2.0, A);
    GaussianPrimitive primitive2 = GaussianPrimitive(x2,y2,z2, 13.0, B);

    //OverlapIntegrator integrator;
    KineticIntegrator integrator;
    //cout << integrator.computeIntegral(primitive1, primitive2) << endl;









    return 0;
}

