#include "gaussianprimitive.h"
#include <cmath>

using std::pow;
using std::exp;
using arma::vec;

GaussianPrimitive::GaussianPrimitive(int    i,
                                     int    j,
                                     int    k,
                                     double a,
                                     vec    nucleusPosition) :
        m_xExponent(i),
        m_yExponent(j),
        m_zExponent(k),
        m_exponent(a),
        m_nucleusPosition(nucleusPosition) {
}

void GaussianPrimitive::setCoefficient(double coefficient) {
    m_coefficient = coefficient;
}

double GaussianPrimitive::evaluate(vec r) {
    const double x              = r(0) - m_nucleusPosition(0);
    const double y              = r(1) - m_nucleusPosition(1);
    const double z              = r(2) - m_nucleusPosition(2);
    const double rNormSquared   = x*x + y*y + z*z;
    return pow(x,m_xExponent) *
           pow(y,m_yExponent) *
           pow(z,m_zExponent) *
           exp(-m_exponent * rNormSquared);
}
