#include "hermitegaussianintegral.h"

using arma::cube;
using arma::zeros;

HermiteGaussianIntegral::HermiteGaussianIntegral() :
        m_t(0),
        m_u(0),
        m_v(0),
        m_nucleusPosition(zeros<vec>(3)),
        m_coefficients(zeros<cube>(5,5,5)) {
}

HermiteGaussianIntegral::setupCoefficients(GaussianPrimitive& primitive1,
                                           GaussianPrimitive& primitive2,
                                           arma::vec nucleusPosition) {
    m_nucleusPosition = nucleusPosition;
    double  alpha   = primitive1.exponent();
    double  beta    = primitive2.exponent();
    double  p       = alpha + beta;
    vec     P       = (alpha * primitive1.nucleusPosition() + beta * primitive2.nucleusPosition()) / p;
    vec     m_PC    = P - nucleusPosition;
    int     m_t     = primitive1.xExponent() + primitive2.xExponent();
    int     m_u     = primitive1.yExponent() + primitive2.yExponent();
    int     m_v     = primitive1.zExponent() + primitive2.zExponent();



}
