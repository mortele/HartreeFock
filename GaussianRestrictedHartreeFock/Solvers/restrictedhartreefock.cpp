#include "restrictedhartreefock.h"
#include <cassert>

using std::cout;
using std::endl;
using arma::eye;
using arma::zeros;
using arma::mat;
using arma::vec;
using arma::field;


RestrictedHartreeFock::RestrictedHartreeFock(System* system) :
        m_system(system) {

    m_numberOfBasisFunctions    = m_system->getNumberOfBasisFunctions();
    m_numberOfElectrons = 0;
    for (Atom* atom : system->getAtoms()) {
        m_numberOfElectrons += atom->getNumberOfElectrons();
    }
    m_basis                     = m_system->getBasis();
    m_epsilon                   = zeros(m_numberOfBasisFunctions);
    m_epsilonOld                = zeros(m_numberOfBasisFunctions);
    m_fockMatrix                = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTilde           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrix         = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixTilde    = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_densityMatrix             = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_overlapMatrix             = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_transformationMatrix      = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_oneBodyMatrixElements     = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);

    m_twoBodyMatrixElements.set_size(m_numberOfBasisFunctions,
                                     m_numberOfBasisFunctions);
    for (int i = 0; i < m_numberOfBasisFunctions; i++) {
        for (int j = 0; j < m_numberOfBasisFunctions; j++) {
            m_twoBodyMatrixElements(i,j) = zeros(m_numberOfBasisFunctions,
                                                 m_numberOfBasisFunctions);
        }
    }
}

void RestrictedHartreeFock::setup() {
    assert(m_numberOfElectrons > 0 && m_numberOfElectrons % 2 == 0);

    setupOverlapMatrix();
    diagonalizeOverlapMatrix();
    setupOneBodyMatrixElements();
    setupTwoBodyMatrixElements();
    computeDensityMatrix();
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

void RestrictedHartreeFock::diagonalizeOverlapMatrix() {
    vec s;
    mat A;
    arma::eig_sym(s, A, m_overlapMatrix);
    m_transformationMatrix = A * arma::diagmat(1.0 / sqrt(s));
}

void RestrictedHartreeFock::normalizeCoefficientMatrix() {
    for (int k = 0; k < m_numberOfElectrons/2; k++) {
        double normalizationFactor = 0;
        for (int p = 0; p < m_numberOfBasisFunctions; p++) {
            for (int q = 0; q < m_numberOfBasisFunctions; q++) {
                normalizationFactor += m_coefficientMatrix(p,k) *
                                       m_coefficientMatrix(q,k) *
                                       m_overlapMatrix(p,q);
            }
        }
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

void RestrictedHartreeFock::printInitialInfo() {
    printf(" ================== Starting SCF iterations ================= \n");
    printf(" => Maximum iterations:    %-10d \n", (int) m_maximumIterations);
    printf(" => Convergence criterion: %-10g \n", m_convergenceCriterion);
    printf(" => Total basis size:      %-10d \n", (int) m_system->getBasis().size());
    printf(" => Number of atoms:       %-10d \n", (int) m_system->getAtoms().size());
    printf("      ------------------------------------------------------- \n");
    for (Atom* atom : m_system->getAtoms()) {
        printf("      | %-20s (%5.3f, %5.3f, %5.3f)          | \n", atom->getInfo().c_str(),
                                                          atom->getPosition()(0),
                                                          atom->getPosition()(1),
                                                          atom->getPosition()(2));
    }
    printf("      ------------------------------------------------------- \n\n");
    printf(" ============================================================ \n");
    printf(" %15s %20s %20s \n", "Iteration", "Energy", "Convergence");
    printf(" ------------------------------------------------------------ \n");
}

void RestrictedHartreeFock::printIterationInfo(int iteration) {
    if (iteration != 1 && iteration % 20 == 0) {
        printf(" ------------------------------------------------------------ \n");
        printf(" %15s %20s %20s \n", "Iteration", "Energy", "Convergence");
        printf(" ------------------------------------------------------------ \n");
    }
    printf(" %15d %20.9g %20.9g \n", iteration, m_hartreeFockEnergy, m_convergenceTest);
}

void RestrictedHartreeFock::printFinalInfo() {
    printf(" ============================================================ \n");
    if (m_reachedSelfConsistency) {
        printf("\n Self consistency SUCCESFULLY reached. \n\n");
        printf(" => Iterations used:        %30d   \n",  m_iterationsUsed);
        printf(" => Final convergence test: %30.16g \n", m_convergenceTest);
        printf(" => Final energy:           %30.16g \n", m_hartreeFockEnergy);
    } else {
        printf("\n Self consistency -> NOT <- reached. \n\n");
        printf(" => Iterations used:        %30d    \n",   m_iterationsUsed);
        printf(" => Final convergence test: %30.16g  \n",   m_convergenceTest);
        printf(" => Final energy:           %30.16g  \n", m_hartreeFockEnergy);
    }
    printf(" ============================================================ \n");
}

double RestrictedHartreeFock::twoBodyMatrixElementsAntiSymmetric(int p,
                                                                 int q,
                                                                 int r,
                                                                 int s) {
    return 2 * m_twoBodyMatrixElements(p,r)(q,s) - m_twoBodyMatrixElements(p,r)(s,q);
}

void RestrictedHartreeFock::computeDensityMatrix() {
    m_densityMatrix = 2 * m_coefficientMatrix * m_coefficientMatrix.t();
}

void RestrictedHartreeFock::setupOverlapMatrix() {
    for (int p = 0; p < m_numberOfBasisFunctions; p++) {
        for (int q = 0; q < m_numberOfBasisFunctions; q++) {
            m_overlapMatrix(p,q) = m_system->overlapIntegral(p,q);
        }
    }
}


void RestrictedHartreeFock::setupTwoBodyMatrixElements() {
    for (int p = 0; p < m_numberOfBasisFunctions; p++)
    for (int q = 0; q < m_numberOfBasisFunctions; q++)
    for (int r = p; r < m_numberOfBasisFunctions; r++)
    for (int s = q; s < m_numberOfBasisFunctions; s++) {
        m_twoBodyMatrixElements(p,q)(r,s) = m_system->twoBodyElements(p,r,q,s);
    }

    for(int p = 0; p < m_numberOfBasisFunctions; p++)
    for(int q = 0; q < m_numberOfBasisFunctions; q++)
    for(int r = p; r < m_numberOfBasisFunctions; r++)
    for(int s = q; s < m_numberOfBasisFunctions; s++) {
        double pqrs = m_twoBodyMatrixElements(p,q)(r,s);
        m_twoBodyMatrixElements(r,s)(p,q) = pqrs;
        m_twoBodyMatrixElements(r,q)(p,s) = pqrs;
        m_twoBodyMatrixElements(p,s)(r,q) = pqrs;
        m_twoBodyMatrixElements(q,p)(s,r) = pqrs;
        m_twoBodyMatrixElements(s,p)(q,r) = pqrs;
        m_twoBodyMatrixElements(q,r)(s,p) = pqrs;
        m_twoBodyMatrixElements(s,r)(q,p) = pqrs;
    }
}

void RestrictedHartreeFock::setupOneBodyMatrixElements() {
    for (int p = 0; p < m_numberOfBasisFunctions; p++) {
        for (int q = 0; q < m_numberOfBasisFunctions; q++) {
            m_oneBodyMatrixElements(p,q) = m_system->oneBodyElements(p,q);
        }
    }
}

double RestrictedHartreeFock::solve(double  convergenceCriteria,
                                    int     maximumIterations) {
    m_convergenceCriterion      = convergenceCriteria;
    m_maximumIterations         = maximumIterations;
    m_reachedSelfConsistency    = false;

    setup();
    if (! m_silent) printInitialInfo();
    for (int iteration = 1; iteration < m_maximumIterations+1; iteration++) {
        selfConsistentFieldIteration();

        if (iteration != 0) {
            m_convergenceTest = arma::abs(m_epsilon - m_epsilonOld).max();
            //m_convergenceTest = arma::sum(arma::abs(m_epsilon - m_epsilonOld)) / m_epsilon.n_elem;
            if(m_convergenceTest < m_convergenceCriterion) {
                m_reachedSelfConsistency    = true;
                m_iterationsUsed            = iteration;
                if (! m_silent) printIterationInfo(iteration);
                break;
            }
        }
        computeHartreeFockEnergy();
        m_epsilonOld = m_epsilon;
        if (! m_silent) printIterationInfo(iteration);
    }
    if (! m_silent) printFinalInfo();
    return m_hartreeFockEnergy;
}

double RestrictedHartreeFock::solveSilently(double  convergenceCriterion,
                                            int     maximumIterations) {
    m_silent = true;
    return solve(convergenceCriterion, maximumIterations);
}





















