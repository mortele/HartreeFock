#include "overlapintegrator.h"
#include <armadillo>
#include <cmath>

using std::exp;
using std::sqrt;
using arma::vec;
using arma::dot;

OverlapIntegrator::OverlapIntegrator() :
        m_hermiteGaussian(HermiteGaussian()) {
}

double OverlapIntegrator::computeIntegral(GaussianPrimitive& primitive1,
                                          GaussianPrimitive& primitive2) {

    m_hermiteGaussian.set(primitive1, primitive2);
    const double exponentSum    = m_hermiteGaussian.getExponentSum();
    const double sqrtPiOverP    = sqrt(M_PI / exponentSum);
    const double integralTerm   = sqrtPiOverP * sqrtPiOverP * sqrtPiOverP;
    const int    xExponent1     = primitive1.xExponent();
    const int    yExponent1     = primitive1.yExponent();
    const int    zExponent1     = primitive1.zExponent();
    const int    xExponent2     = primitive2.xExponent();
    const int    yExponent2     = primitive2.yExponent();
    const int    zExponent2     = primitive2.zExponent();
    const double Ex = m_hermiteGaussian.getCoefficientX(xExponent1, xExponent2);
    const double Ey = m_hermiteGaussian.getCoefficientY(yExponent1, yExponent2);
    const double Ez = m_hermiteGaussian.getCoefficientZ(zExponent1, zExponent2);

    return integralTerm * Ex * Ey * Ez;
}
