#include "Integrators/kineticintegrator.h"
#include "Integrators/overlapintegrator.h"
#include "Orbitals/gaussianprimitive.h"

using arma::zeros;
using arma::mat;
using arma::vec;
using std::cout;
using std::endl;

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
    const int ix = primitive1.xExponent();
    const int iy = primitive1.yExponent();
    const int iz = primitive1.zExponent();

    const int jx = primitive2.xExponent();
    const int jy = primitive2.yExponent();
    const int jz = primitive2.zExponent();

    primitive2.adjustExponentX(2);
    primitive2.adjustExponentY(2);
    primitive2.adjustExponentZ(2);
    m_overlapIntegrator.computeIntegral(primitive1, primitive2);
    primitive2.adjustExponentX(-2);
    primitive2.adjustExponentY(-2);
    primitive2.adjustExponentZ(-2);

    m_overlapIntegrals(0) = m_overlapIntegrator.getIntegralIndicesDimension(ix,jx,0);
    m_overlapIntegrals(1) = m_overlapIntegrator.getIntegralIndicesDimension(iy,jy,1);
    m_overlapIntegrals(2) = m_overlapIntegrator.getIntegralIndicesDimension(iz,jz,2);

    for (int dimension = 0; dimension < 3; dimension++) {
        for (int adjustment = -2; adjustment <= 4; adjustment+=4) {
            computeAdjustedOverlapIntegral(dimension, adjustment);
        }
    }
    for (int dimension = 0; dimension < 3; dimension++) {
        computeT(dimension);
    }

    return - 0.5 * (
           m_T(0)                * m_overlapIntegrals(1) * m_overlapIntegrals(2) +
           m_overlapIntegrals(0) * m_T(1)                * m_overlapIntegrals(2) +
           m_overlapIntegrals(0) * m_overlapIntegrals(1) * m_T(2));
}

void KineticIntegrator::computeAdjustedOverlapIntegral(int dimension,
                                                       int adjustment) {

    const int j = (adjustment == -2) ? 0 : 1;
    const int i2 = m_primitive2.getExponentDimension(dimension);
    if (i2+adjustment >= 0) {
        const int i = m_primitive1.getExponentDimension(dimension);
        m_adjustedOverlapIntegrals(dimension,j) = m_overlapIntegrator.getIntegralIndicesDimension(i,i2+adjustment,dimension);
    } else {
        m_adjustedOverlapIntegrals(dimension,j) = 0;
    }
}

void KineticIntegrator::computeT(int dimension) {
    double beta     = m_primitive2.exponent();
    int    j        = m_primitive2.getExponentDimension(dimension);
    m_T(dimension)  = 4*beta*beta    * m_adjustedOverlapIntegrals(dimension,1)   -
                      2*beta*(2*j+1) * m_overlapIntegrals(dimension)             +
                      j*(j-1)        * m_adjustedOverlapIntegrals(dimension,0);
}




