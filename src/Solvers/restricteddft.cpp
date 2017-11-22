#include "restricteddft.h"
#include "Integrators/numericalintegrator.h"
#include "ExchangeCorrelationFunctionals/localdensityapproximation.h"
#include <cassert>

using std::cout;
using std::endl;
using std::sqrt;
using std::string;
using arma::eye;
using arma::zeros;
using arma::mat;
using arma::vec;
using arma::randu;
using arma::field;


RestrictedDFT::RestrictedDFT(System* system) :
        HartreeFock(system) {
    m_dft = true;

    m_epsilon                   = zeros(m_numberOfBasisFunctions);
    m_epsilonOld                = zeros(m_numberOfBasisFunctions);
    m_fockMatrix                = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTilde           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrix         = eye(m_numberOfBasisFunctions, m_numberOfElectrons/2);
    m_coefficientMatrixTilde    = eye(m_numberOfBasisFunctions, m_numberOfElectrons/2);
    m_densityMatrix             = 2 * m_coefficientMatrix * m_coefficientMatrix.t();
    m_xcMatrix                  = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_numericalIntegrator       = new NumericalIntegrator(system, &m_densityMatrix);
}

void RestrictedDFT::setFunctional(std::string name) {
    if (name == "LDA") {
        m_xcFunctional = new LocalDensityApproximation(m_system, &m_densityMatrix);
    } else {
        cout << "Exchange-correlation functional <" << name << "> not found." << endl;
        exit(1);
    }
    m_numericalIntegrator->setFunctional(m_xcFunctional);
}



double RestrictedDFT::computeEnergyX() {
    return m_numericalIntegrator->integrateExchangeCorrelationEnergy();
}

void RestrictedDFT::setup() {
    assert(m_numberOfElectrons > 0 && m_numberOfElectrons % 2 == 0);

    setupOverlapMatrix();
    diagonalizeOverlapMatrix();
    setupOneBodyElements();
    setupTwoBodyElements();
    computeDensityMatrix();
    computeXcMatrix();
    m_nucleusNucleusInteractionEnergy = m_system->nucleusNucleusInteractionEnergy();

}

void RestrictedDFT::computeXcMatrix() {
    //m_densityMatrix = m_densityMatrix*m_overlapMatrix;
    for(int p = 0; p < m_numberOfBasisFunctions; p++) {
        for(int q = p; q < m_numberOfBasisFunctions; q++) {
            m_xcMatrix(p,q) = Vxc(p,q);
            if (p != q) {
                m_xcMatrix(q,p) = m_xcMatrix(p,q);
            }
        }
    }
    //m_densityMatrix = m_densityMatrix*inv(m_overlapMatrix);
}

void RestrictedDFT::computeFockMatrix() {
    for(int p = 0; p < m_numberOfBasisFunctions; p++)
    for(int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_fockMatrix(p,q) = m_oneBodyMatrixElements(p,q) + m_xcMatrix(p,q);
        //const double vxc = Vxc(p,q);
        for(int r = 0; r < m_numberOfBasisFunctions; r++)
        for(int s = 0; s < m_numberOfBasisFunctions; s++) {
            //m_fockMatrix(p,q) += 0.5 * m_densityMatrix(s,r) * (2*twoBodyMatrixElements(p,r,q,s) - twoBodyMatrixElements(p,r,s,q));
            m_fockMatrix(p,q) += m_densityMatrix(s,r) * twoBodyMatrixElements(p,r,q,s);
        }
    }
}

void RestrictedDFT::computeDensityMatrix() {
    if (m_smoothing) {
        double a = m_smoothingFactor;
        mat densityMatrixTmp = 2 * m_coefficientMatrix * m_coefficientMatrix.t();
        m_densityMatrix      = a * m_densityMatrix + (1.0 - a) * densityMatrixTmp;
    } else {
        m_densityMatrix = 2 * m_coefficientMatrix * m_coefficientMatrix.t();
    }
}

void RestrictedDFT::diagonalizeFockMatrix() {
    const mat& A = m_transformationMatrix;
    m_fockMatrixTilde = A.t() * m_fockMatrix * A;
    arma::eig_sym(m_epsilon, m_coefficientMatrixTilde, m_fockMatrixTilde);
    m_coefficientMatrix = A * m_coefficientMatrixTilde.submat(0,0,m_numberOfBasisFunctions-1, m_numberOfElectrons/2-1);
    normalizeCoefficientMatrix();
}

void RestrictedDFT::normalizeCoefficientMatrix() {
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

void RestrictedDFT::selfConsistentFieldIteration() {
    computeXcMatrix();
    computeFockMatrix();
    diagonalizeFockMatrix();
    computeDensityMatrix();
    //computeDensityIntegral();
}

void RestrictedDFT::computeHartreeFockEnergy() {
    m_xcEnergy           = m_numericalIntegrator->integrateExchangeCorrelationEnergy();
    m_hartreeFockEnergy  = 0;
    m_hartreeFockEnergy += m_xcEnergy;
    m_hartreeFockEnergy += m_nucleusNucleusInteractionEnergy;

    m_oneElectronEnergy = 0;
    m_twoElectronEnergy = 0;

    for (int p = 0; p < m_numberOfBasisFunctions; p++)
    for (int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_oneElectronEnergy += m_densityMatrix(p,q) * m_oneBodyMatrixElements(p,q);
        m_hartreeFockEnergy += m_densityMatrix(p,q) * m_oneBodyMatrixElements(p,q);
        for (int r = 0; r < m_numberOfBasisFunctions; r++)
        for (int s = 0; s < m_numberOfBasisFunctions; s++) {
            double pqrs = 0.5*twoBodyMatrixElements(p,r,q,s);
            m_hartreeFockEnergy += pqrs * m_densityMatrix(s,r) * m_densityMatrix(p,q);
            m_twoElectronEnergy += pqrs * m_densityMatrix(s,r) * m_densityMatrix(p,q);
        }
    }
}

void RestrictedDFT::storeEnergy() {
    m_epsilonOld = m_epsilon;
}

void RestrictedDFT::computeDensityIntegral() {
    m_densityIntegral = m_numericalIntegrator->integrateDensity();
}

double RestrictedDFT::twoBodyMatrixElements(int p,
                                            int q,
                                            int r,
                                            int s) {
    return m_twoBodyMatrixElements(p,q)(r,s);
}

double RestrictedDFT::twoBodyMatrixElementsAntiSymmetric(int p,
                                                         int q,
                                                         int r,
                                                         int s) {
    return 2 * m_twoBodyMatrixElements(p,r)(q,s) - m_twoBodyMatrixElements(p,r)(s,q);
}

double RestrictedDFT::Vxc(int p, int q) {
    return m_numericalIntegrator->integrateExchangeCorrelationPotential(p,q);
}

double RestrictedDFT::Exc(int p, int q) {
    //return m_numericalIntegrator->integrateExchangeCorrelationEnergy(p,q);
    //return m_numericalIntegrator->integrateExchangeCorrelationPotential(p,q);
}

double RestrictedDFT::convergenceTest() {
    //return arma::abs(m_epsilon - m_epsilonOld).max();
    return arma::sum(arma::abs(m_epsilon - m_epsilonOld)) / m_epsilon.n_elem;
}
