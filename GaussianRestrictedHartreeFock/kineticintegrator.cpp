#include "kineticintegrator.h"

using arma::zeros;
using arma::mat;
using arma::vec;

KineticIntegrator::KineticIntegrator() :
        m_overlapIntegrals        (zeros<vec>(3)),
        m_T                       (zeros<vec>(3)),
        m_adjustedOverlapIntegrals(zeros<mat>(3,2)),
        m_overlapIntegrator       (OverlapIntegrator()),
        m_primitive1              (GaussianPrimitive(0,0,0,0,zeros<vec>(3))),
        m_primitive2              (GaussianPrimitive(0,0,0,0,zeros<vec>(3))) {
}

double KineticIntegrator::computeIntegral(GaussianPrimitive& primitive1,
                                          GaussianPrimitive& primitive2) {

    m_primitive1 = primitive1;
    m_primitive2 = primitive2;

    m_overlapIntegrator.computeIntegral(primitive1, primitive2);
    m_overlapIntegrals(0) = m_overlapIntegrator.getIntegralX();
    m_overlapIntegrals(1) = m_overlapIntegrator.getIntegralY();
    m_overlapIntegrals(2) = m_overlapIntegrator.getIntegralZ();

    for (int dimension = 0; dimension < 3; dimension++) {
        for (int adjustment = -2; adjustment <= 4; adjustment+=4) {
            computeAdjustedOverlapIntegral(dimension, adjustment);
        }
    }
    for (int dimension = 0; dimension < 3; dimension++) {
        computeT(dimension);
    }
    return m_T(0)                * m_overlapIntegrals(1) * m_overlapIntegrals(2) +
           m_overlapIntegrals(0) * m_T(1)                * m_overlapIntegrals(2) +
           m_overlapIntegrals(0) * m_overlapIntegrals(1) * m_T(2);
}

void KineticIntegrator::computeAdjustedOverlapIntegral(int dimension,
                                                       int adjustment) {

    const int j = (adjustment == -2) ? 0 : 1;
    if (m_primitive2.getExponentDimension(dimension)+adjustment >= 0) {
        m_primitive2.adjustExponentDimension(adjustment, dimension);
        m_overlapIntegrator.computeIntegral(m_primitive1, m_primitive2);

        m_adjustedOverlapIntegrals(dimension,j) = m_overlapIntegrator.getIntegralDimension(dimension);
        m_primitive2.adjustExponentDimension(-adjustment, dimension);
    } else {
        m_adjustedOverlapIntegrals(dimension,j) = 0;
    }
}

void KineticIntegrator::computeT(int dimension) {
    double beta     = m_primitive2.exponent();
    int    j        = m_primitive2.getExponentDimension(dimension);
    m_T(dimension)  = 4*beta*beta    * m_adjustedOverlapIntegrals(dimension,1)   +
                      2*beta*(2*j+1) * m_overlapIntegrals(dimension)             +
                      j*(j-1)        * m_adjustedOverlapIntegrals(dimension,0);
}




