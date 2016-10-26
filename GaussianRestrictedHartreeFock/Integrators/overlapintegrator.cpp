#include "overlapintegrator.h"
#include <armadillo>
#include <cmath>

using std::exp;
using std::sqrt;
using std::cout;
using std::endl;
using arma::vec;
using arma::cube;
using arma::dot;

OverlapIntegrator::OverlapIntegrator() :
        m_Ex(0),
        m_Ey(0),
        m_Ez(0),
        m_sqrtPiOverP(0),
        m_hermiteGaussian(HermiteGaussian()) {
}

double OverlapIntegrator::computeIntegral(GaussianPrimitive& primitive1,
                                          GaussianPrimitive& primitive2) {

    m_hermiteGaussian.setupCoefficients(primitive1, primitive2);
    const double exponentSum    = m_hermiteGaussian.getExponentSum();
    m_sqrtPiOverP               = sqrt(M_PI / exponentSum);
    const int    xExponent1     = primitive1.xExponent();
    const int    yExponent1     = primitive1.yExponent();
    const int    zExponent1     = primitive1.zExponent();

    const int    xExponent2     = primitive2.xExponent();
    const int    yExponent2     = primitive2.yExponent();
    const int    zExponent2     = primitive2.zExponent();

    m_Ex = m_sqrtPiOverP*m_hermiteGaussian.getCoefficientX(xExponent1, xExponent2);
    m_Ey = m_sqrtPiOverP*m_hermiteGaussian.getCoefficientY(yExponent1, yExponent2);
    m_Ez = m_sqrtPiOverP*m_hermiteGaussian.getCoefficientZ(zExponent1, zExponent2);
    return  m_Ex * m_Ey * m_Ez;
}

double OverlapIntegrator::getIntegralIndicesDimension(int i, int j, int dimension) {
    return m_sqrtPiOverP * m_hermiteGaussian.getCoefficientDimension(i,j,dimension);
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

