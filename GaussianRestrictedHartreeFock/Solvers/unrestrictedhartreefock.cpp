#include "unrestrictedhartreefock.h"

using std::cout;
using std::endl;
using arma::eye;
using arma::zeros;
using arma::mat;
using arma::vec;
using arma::field;


UnrestrictedHartreeFock::UnrestrictedHartreeFock(System* system) :
        HartreeFock(system) {

    m_numberOfSpinUpElectrons   = m_system->getNumberOfSpinUpElectrons();
    m_numberOfSpinDownElectrons = m_system->getNumberOfSpinDownElectrons();

    m_epsilonUp                   = zeros(m_numberOfSpinUpElectrons);
    m_epsilonOldUp                = zeros(m_numberOfSpinUpElectrons);
    m_fockMatrixUp                = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTildeUp           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixUp         = zeros(m_numberOfBasisFunctions, m_numberOfSpinUpElectrons);
    m_coefficientMatrixTildeUp    = zeros(m_numberOfBasisFunctions, m_numberOfSpinUpElectrons);
    m_densityMatrixUp             = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);

    m_epsilonDown                 = zeros(m_numberOfSpinDownElectrons);
    m_epsilonOldDown              = zeros(m_numberOfSpinDownElectrons);
    m_fockMatrixDown              = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTildeDown         = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixDown       = zeros(m_numberOfBasisFunctions, m_numberOfSpinDownElectrons);
    m_coefficientMatrixTildeDown  = zeros(m_numberOfBasisFunctions, m_numberOfSpinDownElectrons);
    m_densityMatrixDown           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
}

void UnrestrictedHartreeFock::setup() {
    setupOverlapMatrix();
    diagonalizeOverlapMatrix();
    setupOneBodyElements();
    setupTwoBodyElements();
    computeDensityMatrices();
    m_nucleusNucleusInteractionEnergy = m_system->nucleusNucleusInteractionEnergy();
}

void UnrestrictedHartreeFock::computeFockMatrices() {
    for(int p = 0; p < m_numberOfBasisFunctions; p++)
    for(int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_fockMatrixUp(p,q)   = m_oneBodyMatrixElements(p,q);
        m_fockMatrixDown(p,q) = m_oneBodyMatrixElements(p,q);

        for(int r = 0; r < m_numberOfBasisFunctions; r++)
        for(int s = 0; s < m_numberOfBasisFunctions; s++) {
            const double prqs = m_twoBodyMatrixElements(p,r)(q,s);
            const double prsq = m_twoBodyMatrixElements(p,r)(s,q);
            m_fockMatrixUp(p,q)   += m_densityMatrixUp(s,r)   * (prqs-prsq)
                                   + m_densityMatrixDown(s,r) * prqs;
            m_fockMatrixDown(p,q) += m_densityMatrixDown(s,r) * (prqs-prsq)
                                   + m_densityMatrixUp(s,r)   * prqs;
        }
    }
}

void UnrestrictedHartreeFock::computeDensityMatrices() {
    m_densityMatrixUp   = m_coefficientMatrixUp   * m_coefficientMatrixUp.t();
    m_densityMatrixDown = m_coefficientMatrixDown * m_coefficientMatrixDown.t();
}

void UnrestrictedHartreeFock::diagonalizeFockMatrices() {
    const mat& A = m_transformationMatrix;
    m_fockMatrixTildeUp   = A.t() * m_fockMatrixUp   * A;
    m_fockMatrixTildeDown = A.t() * m_fockMatrixDown * A;
    arma::eig_sym(m_epsilonUp,   m_coefficientMatrixTildeUp,   m_fockMatrixTildeUp);
    arma::eig_sym(m_epsilonDown, m_coefficientMatrixTildeDown, m_fockMatrixTildeDown);
    m_coefficientMatrixUp   = A * m_coefficientMatrixTildeUp.submat(0,0,m_numberOfBasisFunctions-1, m_numberOfSpinUpElectrons-1);
    m_coefficientMatrixDown = A * m_coefficientMatrixTildeDown.submat(0,0,m_numberOfBasisFunctions-1, m_numberOfSpinDownElectrons-1);
    normalizeCoefficientMatrices();
}

void UnrestrictedHartreeFock::normalizeCoefficientMatrices() {
    for (int k = 0; k < m_numberOfSpinUpElectrons; k++) {
        double normalizationFactor = 0;
        for (int p = 0; p < m_numberOfBasisFunctions; p++) {
            for (int q = 0; q < m_numberOfBasisFunctions; q++) {
                normalizationFactor += m_coefficientMatrixUp(p,k) *
                                       m_coefficientMatrixUp(q,k) *
                                       m_overlapMatrix(p,q);
            }
        }
        m_coefficientMatrixUp.col(k) = m_coefficientMatrixUp.col(k) /
                                       normalizationFactor;
    }
    for (int k = 0; k < m_numberOfSpinDownElectrons; k++) {
        double normalizationFactor = 0;
        for (int p = 0; p < m_numberOfBasisFunctions; p++) {
            for (int q = 0; q < m_numberOfBasisFunctions; q++) {
                normalizationFactor += m_coefficientMatrixDown(p,k) *
                                       m_coefficientMatrixDown(q,k) *
                                       m_overlapMatrix(p,q);
            }
        }
        m_coefficientMatrixDown.col(k) = m_coefficientMatrixDown.col(k) /
                                       normalizationFactor;
    }
}

void UnrestrictedHartreeFock::selfConsistentFieldIteration() {
    computeFockMatrices();
    diagonalizeFockMatrices();
    computeDensityMatrices();
}

void UnrestrictedHartreeFock::computeHartreeFockEnergy() {
    m_hartreeFockEnergy = 0;

    for (int p = 0; p < m_numberOfBasisFunctions; p++)
    for (int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_hartreeFockEnergy += m_oneBodyMatrixElements(p,q)
                               * (m_densityMatrixUp(p,q) + m_densityMatrixDown(p,q));
        m_hartreeFockEnergy += m_fockMatrixUp(p,q)   * m_densityMatrixUp(p,q);
        m_hartreeFockEnergy += m_fockMatrixDown(p,q) * m_densityMatrixDown(p,q);
    }
    m_hartreeFockEnergy *= 0.5;
    m_hartreeFockEnergy += m_nucleusNucleusInteractionEnergy;
}

void UnrestrictedHartreeFock::storeEnergy() {
    m_epsilonOldUp   = m_epsilonUp;
    m_epsilonOldDown = m_epsilonDown;
}

double UnrestrictedHartreeFock::convergenceTest() {
    double sumUp    = arma::sum(arma::abs(m_epsilonUp   - m_epsilonOldUp))   / m_epsilonUp.n_elem;
    double sumDown  = arma::sum(arma::abs(m_epsilonDown - m_epsilonOldDown)) / m_epsilonDown.n_elem;
    return 0.5 * (sumUp + sumDown);
}
