#include "Integrators/coulombnucleusintegrator.h"

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

CoulombNucleusIntegrator::CoulombNucleusIntegrator() :
        m_nucleusPosition(zeros<vec>(3)),
        m_hermiteGaussian(HermiteGaussian()),
        m_hermiteGaussianIntegral(HermiteGaussianIntegral()) {

}

void CoulombNucleusIntegrator::setNucleusPosition(arma::vec nucleusPosition) {
    m_nucleusPosition = nucleusPosition;
}

double CoulombNucleusIntegrator::computeIntegral(GaussianPrimitive& primitive1,
                                                 GaussianPrimitive& primitive2) {
    m_hermiteGaussianIntegral.setupCoefficients(primitive1, primitive2, m_nucleusPosition);
}

double CoulombNucleusIntegrator::computeIntegral(GaussianPrimitive& primitive1,
                                                 GaussianPrimitive& primitive2,
                                                 arma::vec          nucleusPosition) {
    setNucleusPosition(nucleusPosition);
    return computeIntegral(primitive1, primitive2);
}
