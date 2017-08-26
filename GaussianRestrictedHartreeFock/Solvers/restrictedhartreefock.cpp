#include "restrictedhartreefock.h"
#include <cassert>

using std::cout;
using std::endl;
using std::sqrt;
using arma::eye;
using arma::zeros;
using arma::mat;
using arma::vec;
using arma::field;


RestrictedHartreeFock::RestrictedHartreeFock(System* system) :
        HartreeFock(system) {
    m_epsilon                   = zeros(m_numberOfBasisFunctions);
    m_epsilonOld                = zeros(m_numberOfBasisFunctions);
    m_fockMatrix                = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTilde           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrix         = zeros(m_numberOfBasisFunctions, m_numberOfElectrons/2);
    m_coefficientMatrixTilde    = zeros(m_numberOfBasisFunctions, m_numberOfElectrons/2);
    m_densityMatrix             = 2 * m_coefficientMatrix * m_coefficientMatrix.t(); // zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
}

void RestrictedHartreeFock::setup() {
    assert(m_numberOfElectrons > 0 && m_numberOfElectrons % 2 == 0);

    setupOverlapMatrix();
    diagonalizeOverlapMatrix();
    setupOneBodyElements();
    setupTwoBodyElements();
    //computeDensityMatrix();
    m_nucleusNucleusInteractionEnergy = m_system->nucleusNucleusInteractionEnergy();
}

void RestrictedHartreeFock::computeFockMatrix() {
    for(int p = 0; p < m_numberOfBasisFunctions; p++)
    for(int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_fockMatrix(p,q) = m_oneBodyMatrixElements(p,q);

        for(int r = 0; r < m_numberOfBasisFunctions; r++)
        for(int s = 0; s < m_numberOfBasisFunctions; s++) {
            m_fockMatrix(p,q) += 0.5 * m_densityMatrix(s,r) * twoBodyMatrixElementsAntiSymmetric(p,q,r,s);
        }
    }
}

void RestrictedHartreeFock::diagonalizeFockMatrix() {
    const mat& A = m_transformationMatrix;
    m_fockMatrixTilde = A.t() * m_fockMatrix * A;
    arma::eig_sym(m_epsilon, m_coefficientMatrixTilde, m_fockMatrixTilde);
    m_coefficientMatrix = A * m_coefficientMatrixTilde.submat(0,0,m_numberOfBasisFunctions-1, m_numberOfElectrons/2-1);
    normalizeCoefficientMatrix();
}

void RestrictedHartreeFock::normalizeCoefficientMatrix() {
    for (int k = 0; k < m_numberOfElectrons/2; k++) {
        double normalizationFactor = 0;
        for (int p = 0; p < m_numberOfBasisFunctions; p++) {
            for (int q = 0; q < m_numberOfBasisFunctions; q++) {
                // Since the coefficient matrix is real and symmetric, we use
                // C(p,k) = C(k,p) here, even though in general this should be
                // C'(k,p), where ' denotes the Hermitian conjugate.
                normalizationFactor += m_coefficientMatrix(p,k) *
                                       m_coefficientMatrix(q,k) *
                                       m_overlapMatrix(p,q);
            }
        }
        normalizationFactor = sqrt(normalizationFactor);
        m_coefficientMatrix.col(k) = m_coefficientMatrix.col(k) /
                                     normalizationFactor;
    }
}

void RestrictedHartreeFock::selfConsistentFieldIteration() {
    computeFockMatrix();
    diagonalizeFockMatrix();
    computeDensityMatrix();
}

void RestrictedHartreeFock::computeHartreeFockEnergy() {
    m_hartreeFockEnergy = 0;

    for (int p = 0; p < m_numberOfBasisFunctions; p++)
    for (int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_hartreeFockEnergy += m_densityMatrix(p,q) * m_oneBodyMatrixElements(p,q);

        for (int r = 0; r < m_numberOfBasisFunctions; r++)
        for (int s = 0; s < m_numberOfBasisFunctions; s++) {
            m_hartreeFockEnergy += twoBodyMatrixElementsAntiSymmetric(p,q,r,s) *
                                   m_densityMatrix(p,q) *
                                   m_densityMatrix(s,r) *
                                   0.25;
        }
    }
    m_hartreeFockEnergy += m_nucleusNucleusInteractionEnergy;
}

void RestrictedHartreeFock::storeEnergy() {
    m_epsilonOld = m_epsilon;
}

double RestrictedHartreeFock::twoBodyMatrixElementsAntiSymmetric(int p,
                                                                 int q,
                                                                 int r,
                                                                 int s) {
    return 2 * m_twoBodyMatrixElements(p,r)(q,s) - m_twoBodyMatrixElements(p,r)(s,q);
}

double RestrictedHartreeFock::convergenceTest() {
    //return arma::abs(m_epsilon - m_epsilonOld).max();
    return arma::sum(arma::abs(m_epsilon - m_epsilonOld)) / m_epsilon.n_elem;
}

void RestrictedHartreeFock::computeDensityMatrix() {
    if (m_smoothing && !m_firstSetup) {
        m_firstSetup = false;
        double a = m_smoothingFactor;
        mat densityMatrixTmp    = 2 * m_coefficientMatrix * m_coefficientMatrix.t();
        m_densityMatrix         = a*densityMatrixTmp * m_densityMatrix + (1-a)*densityMatrixTmp;

    } else {
        m_densityMatrix = 2 * m_coefficientMatrix * m_coefficientMatrix.t();
    }
}





















