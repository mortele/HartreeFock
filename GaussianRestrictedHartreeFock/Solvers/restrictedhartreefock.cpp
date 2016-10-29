#include "restrictedhartreefock.h"

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
    m_numberOfElectrons         = m_system->getAtoms().size();
    m_basis                     = m_system->getBasis();
    m_epsilon                   = zeros<vec>(m_numberOfBasisFunctions);
    m_epsilonOld                = zeros<vec>(m_numberOfBasisFunctions);
    m_fockMatrix                = zeros<mat>(m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);
    m_U                         = zeros<mat>(m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);
    m_densityMatrix             = eye<mat>  (m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);
    m_overlapMatrix             = zeros<mat>(m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);
    m_oneBodyMatrixElements     = zeros<mat>(m_numberOfBasisFunctions,
                                             m_numberOfBasisFunctions);

    m_twoBodyMatrixElements.set_size(m_numberOfBasisFunctions,
                                     m_numberOfBasisFunctions);
    for (int i = 0; i < m_numberOfBasisFunctions; i++) {
        for (int j = 0; j < m_numberOfBasisFunctions; j++) {
            m_twoBodyMatrixElements(i,j) = zeros<mat>(m_numberOfBasisFunctions,
                                                      m_numberOfBasisFunctions);
        }
    }
}


void RestrictedHartreeFock::setup() {
    setupOverlapMatrix();
    setupOneBodyMatrixElements();
    setupTwoBodyMatrixElements();
    m_nucleusNucleusInteractionEnergy = m_system->nucleusNucleusInteractionEnergy();
}

void RestrictedHartreeFock::computeFockMatrix() {
    for(int p = 0; p < m_numberOfBasisFunctions; p++)
    for(int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_fockMatrix(p,q) = m_oneBodyMatrixElements(p,q);

        double rsSum = 0;
        for(int r = 0; r < m_numberOfBasisFunctions; r++)
        for(int s = 0; s < m_numberOfBasisFunctions; s++) {
            rsSum += m_densityMatrix(s,r) * m_twoBodyMatrixElements(p,r)(q,s);
        }
        m_fockMatrix(p,q) += rsSum;
    }
}

void RestrictedHartreeFock::diagonalizeFockMatrix() {
    m_U = m_overlapMatrix * m_U;
    arma::eig_sym(m_epsilon, m_U, m_fockMatrix);
}

void RestrictedHartreeFock::selfConsistentFieldIteration() {
    computeDensityMatrix();
    computeFockMatrix();
    diagonalizeFockMatrix();
}

void RestrictedHartreeFock::computeHartreeFockEnergy() {
    m_hartreeFockEnergy = 0;

    for(int p = 0; p < m_numberOfBasisFunctions; p++) {
        m_hartreeFockEnergy += 2 * m_epsilon(p);
    }
    for(int i = 0; i < m_numberOfBasisFunctions; i++)
    for(int p = 0; p < m_numberOfBasisFunctions; p++)
    for(int q = 0; q < m_numberOfBasisFunctions; q++) {
        double rsSum = 0;
        for(int r = 0; r < m_numberOfBasisFunctions; r++)
        for(int s = 0; s < m_numberOfBasisFunctions; s++) {
                rsSum += m_densityMatrix(r,s)*m_twoBodyMatrixElements(q,s)(p,r);
        }
        m_hartreeFockEnergy -= m_U(q,i) * rsSum * m_U(p,i);
    }
    m_hartreeFockEnergy += m_nucleusNucleusInteractionEnergy;
}

void RestrictedHartreeFock::printInitialInfo() {
    printf(" ================== Starting SCF iterations ================= \n");
    printf(" => Maximum iterations:    %-10d \n", (int) m_maximumIterations);
    printf(" => Convergence criterion: %-10g \n", m_convergenceCriterion);
    printf(" => Number of atoms:       %-10d \n", (int) m_system->getAtoms().size());
    for (Atom* atom : m_system->getAtoms()) {
        printf("    * %-20s (%5.3f, %5.3f, %5.3f) \n", atom->getInfo().c_str(),
                                                          atom->getPosition()(0),
                                                          atom->getPosition()(1),
                                                          atom->getPosition()(2));
    }
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

void RestrictedHartreeFock::computeDensityMatrix() {
    for(int p = 0; p < m_numberOfBasisFunctions; p++) {
        for(int q = 0; q < m_numberOfBasisFunctions; q++) {
            m_densityMatrix(q,p) = 0.0;
            for(int j = 0; j < m_numberOfBasisFunctions; j++) {
                m_densityMatrix(q,p) += 2.0 * m_U(q,j) * m_U(p,j);
            }
        }
    }
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
        m_twoBodyMatrixElements(p,q)(r,s) =         m_system->twoBodyElements(p,q,r,s)
                                            - 0.5 * m_system->twoBodyElements(p,q,s,r);
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
    printInitialInfo();
    for (int iteration = 1; iteration < m_maximumIterations+1; iteration++) {
        selfConsistentFieldIteration();

        if (iteration != 0) {
            m_convergenceTest = arma::abs(m_epsilon - m_epsilonOld).max();
            //m_convergenceTest = arma::sum(arma::abs(m_epsilon - m_epsilonOld)) / m_epsilon.n_elem;
            if(m_convergenceTest < m_convergenceCriterion) {
                m_reachedSelfConsistency    = true;
                m_iterationsUsed            = iteration;
                printIterationInfo(iteration);
                break;
            }
        }
        computeHartreeFockEnergy();
        m_epsilonOld = m_epsilon;
        printIterationInfo(iteration);
    }
    printFinalInfo();
    return m_reachedSelfConsistency;
}





















