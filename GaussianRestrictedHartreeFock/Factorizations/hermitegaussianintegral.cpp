#include "Factorizations/hermitegaussianintegral.h"

using arma::cube;
using arma::zeros;
using arma::vec;
using arma::zeros;


HermiteGaussianIntegral::HermiteGaussianIntegral() :
        m_t(0),
        m_u(0),
        m_v(0),
        m_nucleusPosition(zeros<vec>(3)),
        m_maxExponents(0),
        m_boysFunction(BoysFunction()){
}

int HermiteGaussianIntegral::computeMaximumExponents(GaussianPrimitive& primitive1,
                                                     GaussianPrimitive& primitive2){
    int max = 0;
    for (int dimension=0; dimension<3; dimension++) {
        max = std::max(primitive1.getExponentDimension(dimension), max);
        max = std::max(primitive2.getExponentDimension(dimension), max);
    }
    return max;
}


void HermiteGaussianIntegral::setupCoefficients(GaussianPrimitive& primitive1,
                                                GaussianPrimitive& primitive2,
                                                vec nucleusPosition) {

    m_maxExponents  = computeMaximumExponents(primitive1, primitive2);
    int maxIndex    = 4*m_maxExponents;
    m_coefficients.set_size(maxIndex);
    for (int i=0; i<maxIndex; i++) {
        m_coefficients(i) = zeros<cube>(maxIndex+1, maxIndex+1, maxIndex+1);
    }

    m_nucleusPosition   = nucleusPosition;
    double  alpha       = primitive1.exponent();
    double  beta        = primitive2.exponent();
    double  p           = alpha + beta;
    vec     P           = (alpha * primitive1.nucleusPosition() +
                           beta  * primitive2.nucleusPosition()) / p;
    vec     m_PC        = P - nucleusPosition;
    int     m_t         = primitive1.xExponent() + primitive2.xExponent();
    int     m_u         = primitive1.yExponent() + primitive2.yExponent();
    int     m_v         = primitive1.zExponent() + primitive2.zExponent();
    int     m_tuv       = m_t + m_u + m_v;

    double x = p * arma::dot(m_PC, m_PC);
    m_coefficients(0)(0,0,0) = m_boysFunction.computeAndApplyDownwardRecurrence(x, 40);


}
