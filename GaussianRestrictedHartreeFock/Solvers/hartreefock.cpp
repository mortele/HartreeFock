#include "hartreefock.h"
#include "Atoms/atom.h"

using std::cout;
using std::endl;
using arma::eye;
using arma::zeros;
using arma::mat;
using arma::vec;
using arma::field;


HartreeFock::HartreeFock(System* system) :
        m_system(system) {

    m_numberOfBasisFunctions    = m_system->getNumberOfBasisFunctions();
    m_numberOfElectrons = 0;
    for (Atom* atom : system->getAtoms()) {
        m_numberOfElectrons += atom->getNumberOfElectrons();
    }
    m_basis                     = m_system->getBasis();
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

double HartreeFock::solve(double convergenceCriterion, int maximumIterations) {
    m_convergenceCriterion      = convergenceCriterion;
    m_maximumIterations         = maximumIterations;
    m_reachedSelfConsistency    = false;

    if (! m_setupDone) {
        setup();
        m_setupDone = true;
    }
    if (! m_silent) printInitialInfo();
    for (m_iteration = 0; m_iteration < m_maximumIterations+1; m_iteration++) {
        selfConsistentFieldIteration();

        if (m_iteration != 0) {
            m_convergenceTest = convergenceTest();
            if(m_convergenceTest < m_convergenceCriterion) {
                m_reachedSelfConsistency    = true;
                m_iterationsUsed            = m_iteration;
                if (! m_silent) printIterationInfo(m_iteration);
                break;
            }
        }
        computeHartreeFockEnergy();
        m_electronicHartreeFockEnergy = m_hartreeFockEnergy - m_nucleusNucleusInteractionEnergy;
        storeEnergy();
        if (! m_silent) printIterationInfo(m_iteration);
    }
    if (! m_silent) printFinalInfo();
    return m_hartreeFockEnergy;
}

double HartreeFock::solveSilently(double convergenceCriterion, int maximumIterations) {
    m_silent = true;
    return solve(convergenceCriterion, maximumIterations);
}

void HartreeFock::setupOneBodyElements() {
    for (int p = 0; p < m_numberOfBasisFunctions; p++) {
        for (int q = 0; q < m_numberOfBasisFunctions; q++) {
            m_oneBodyMatrixElements(p,q) = m_system->oneBodyElements(p,q);
        }
    }
}

void HartreeFock::setupTwoBodyElements() {
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

void HartreeFock::setupOverlapMatrix() {
    for (int p = 0; p < m_numberOfBasisFunctions; p++) {
        for (int q = 0; q < m_numberOfBasisFunctions; q++) {
            m_overlapMatrix(p,q) = m_system->overlapIntegral(p,q);
        }
    }
}

void HartreeFock::diagonalizeOverlapMatrix() {
    vec s;
    mat A;
    arma::eig_sym(s, A, m_overlapMatrix);
    m_transformationMatrix = A * arma::diagmat(1.0 / sqrt(s));
}

void HartreeFock::printInitialInfo() {
    printf(" ================== Starting SCF iterations ================= \n");
    printf(" => Maximum iterations:    %-10d \n", (int) m_maximumIterations);
    printf(" => Convergence criterion: %-10g \n", m_convergenceCriterion);
    printf(" => Total basis size:      %-10d \n", (int) m_system->getBasis().size());
    printf(" => Number of atoms:       %-10d \n", (int) m_system->getAtoms().size());
    printf(" => Number of electrons:   %-10d \n", (int) m_system->getNumberOfSpinUpElectrons()+m_system->getNumberOfSpinDownElectrons());
    printf("   ------------------------------------------------------- \n");
    for (Atom* atom : m_system->getAtoms()) {
        printf("   | %-25s (%6.3f, %6.3f, %6.3f)  | \n", atom->getInfo().c_str(),
                                                          atom->getPosition()(0),
                                                          atom->getPosition()(1),
                                                          atom->getPosition()(2));
    }
    printf("   ------------------------------------------------------- \n\n");
    printf(" ============================================================ \n");
    printf(" %15s %20s %20s \n", "Iteration", "Energy", "Convergence");
    printf(" ------------------------------------------------------------ \n");
    fflush(stdout);
}

void HartreeFock::printIterationInfo(int iteration) {
    if (iteration != 0 && iteration % 20 == 0) {
        printf(" ------------------------------------------------------------ \n");
        printf(" %15s %20s %20s \n", "Iteration", "Energy", "Convergence");
        printf(" ------------------------------------------------------------ \n");
        fflush(stdout);
    }
    printf(" %15d %20.9g %20.9g \n", iteration, m_hartreeFockEnergy, m_convergenceTest);
}

void HartreeFock::printFinalInfo() {
    printf(" ============================================================ \n");
    if (m_reachedSelfConsistency) {
        printf("\n Self consistency SUCCESFULLY reached. \n\n");
        printf(" => Iterations used:         %30d   \n",  m_iterationsUsed);
        printf(" => Final convergence test:  %30.16g \n", m_convergenceTest);
        printf(" => Final electronic energy: %30.16g  \n", m_electronicHartreeFockEnergy);
        printf(" => Final energy (eV):       %30.16g \n", m_hartreeFockEnergy*27.21139);
        printf(" => Final energy:            %30.16g \n", m_hartreeFockEnergy);
    } else {
        printf("\n Self consistency -> NOT <- reached. \n\n");
        printf(" => Iterations used:         %30d    \n",   m_iterationsUsed);
        printf(" => Final convergence test:  %30.16g  \n",   m_convergenceTest);
        printf(" => Final electronic energy: %30.16g  \n", m_electronicHartreeFockEnergy);
        printf(" => Final energy (eV):       %30.16g \n", m_hartreeFockEnergy*27.21139);
        printf(" => Final energy:            %30.16g  \n", m_hartreeFockEnergy);
    }
    printf(" ============================================================ \n");
    cout << m_hartreeFockEnergy-m_electronicHartreeFockEnergy << endl;
}



