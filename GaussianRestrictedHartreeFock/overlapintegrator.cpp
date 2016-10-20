#include "overlapintegrator.h"
#include <armadillo>
#include <cmath>

using std::exp;
using std::sqrt;
using arma::vec;
using arma::dot;

OverlapIntegrator::OverlapIntegrator() :
        m_Ex(0),
        m_Ey(0),
        m_Ez(0),
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
    m_Ex = integralTerm*m_hermiteGaussian.getCoefficientX(xExponent1, xExponent2);
    m_Ey = integralTerm*m_hermiteGaussian.getCoefficientY(yExponent1, yExponent2);
    m_Ez = integralTerm*m_hermiteGaussian.getCoefficientZ(zExponent1, zExponent2);
    return  m_Ex * m_Ey * m_Ez;
}

double OverlapIntegrator::getIntegralDimension(int dimension) {
    if (dimension == 0) {
        return getIntegralX();
    } else if (dimension == 1) {
        return getIntegralY();
    } else if (dimension == 2) {
        return getIntegralZ();
    }
}

