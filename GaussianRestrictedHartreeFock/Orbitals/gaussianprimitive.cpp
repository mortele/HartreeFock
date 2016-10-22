#include "Orbitals/gaussianprimitive.h"
#include <cmath>

using std::pow;
using std::exp;
using std::max;
using arma::vec;



GaussianPrimitive::GaussianPrimitive(int    i,
                                     int    j,
                                     int    k,
                                     double a,
                                     vec    nucleusPosition,
                                     double constantTerm) :
        m_xExponent             (i),
        m_yExponent             (j),
        m_zExponent             (k),
        m_maximumAngularMomentum(max(max(i,j),k)),
        m_exponent              (a),
        m_constantTerm          (constantTerm),
        m_nucleusPosition       (nucleusPosition),
        m_coefficient           (1.0) {
}

void GaussianPrimitive::setCoefficient(double coefficient) {
    m_coefficient = coefficient;
}

double GaussianPrimitive::evaluate(vec& r) {
    const double x              = r(0) - m_nucleusPosition(0);
    const double y              = r(1) - m_nucleusPosition(1);
    const double z              = r(2) - m_nucleusPosition(2);
    const double rNormSquared   = x*x + y*y + z*z;
    return m_constantTerm *
           pow(x,m_xExponent) *
           pow(y,m_yExponent) *
           pow(z,m_zExponent) *
           exp(-m_exponent * rNormSquared);
}

GaussianPrimitive GaussianPrimitive::product(GaussianPrimitive& primitive1,
                                             GaussianPrimitive& primitive2) {
    int     xExponent       = primitive1.xExponent() + primitive2.xExponent();
    int     yExponent       = primitive1.yExponent() + primitive2.yExponent();
    int     zExponent       = primitive1.zExponent() + primitive2.zExponent();
    double  alpha           = primitive1. exponent();
    double  beta            = primitive2. exponent();
    double  exponent        = alpha + beta;
    double  nucleusWeight1  = alpha / (alpha + beta);
    double  nucleusWeight2  = beta  / (alpha + beta);
    vec nucleusCenter       = nucleusWeight1 * primitive1.nucleusPosition() +
                              nucleusWeight2 * primitive2.nucleusPosition();
    vec nucleonDistance     = primitive1.nucleusPosition() -
                              primitive2.nucleusPosition();
    double constantTerm     = exp( - alpha * beta / (alpha + beta) *
                                   dot(nucleonDistance, nucleonDistance));
    return GaussianPrimitive(xExponent,
                             yExponent,
                             zExponent,
                             exponent,
                             nucleusCenter,
                             constantTerm);
}

void GaussianPrimitive::adjustExponentX(int increment) {
    m_xExponent += increment;
    m_maximumAngularMomentum = max(m_xExponent, m_maximumAngularMomentum);
}

void GaussianPrimitive::adjustExponentY(int increment) {
    m_yExponent += increment;
    m_maximumAngularMomentum = max(m_yExponent, m_maximumAngularMomentum);
}

void GaussianPrimitive::adjustExponentZ(int increment) {
    m_zExponent += increment;
    m_maximumAngularMomentum = max(m_zExponent, m_maximumAngularMomentum);
}

void GaussianPrimitive::adjustExponentDimension(int increment, int dimension) {
    if (dimension==0) {
        adjustExponentX(increment);
    } else if (dimension==1) {
        adjustExponentY(increment);
    } else if (dimension==2) {
        adjustExponentZ(increment);
    }
}

int GaussianPrimitive::getExponentDimension(int dimension) const {
    if (dimension==0) {
        return xExponent();
    } else if (dimension==1) {
        return yExponent();
    } else if (dimension==2) {
        return zExponent();
    }
}

