#pragma once
#include <armadillo>

class HartreeFock {
public:
    int         m_numberOfBasisFunctions;
    int         m_numberOfElectrons;
    int         m_iteration;
    int         m_maximumIterations;
    int         m_iterationsUsed;
    bool        m_dft = false;
    bool        m_reachedSelfConsistency;
    bool        m_silent = false;
    bool        m_setupDone = false;
    double      m_xcEnergy;
    double      m_oneElectronEnergy;
    double      m_twoElectronEnergy;
    double      m_convergenceCriterion;
    double      m_hartreeFockEnergy;
    double      m_electronicHartreeFockEnergy;
    double      m_nucleusNucleusInteractionEnergy;
    double      m_convergenceTest;
    class System*     m_system;
    arma::mat   m_oneBodyMatrixElements;
    arma::field<arma::mat>           m_twoBodyMatrixElements;
    std::vector<class ContractedGaussian*> m_basis;

    void setupOneBodyElements();
    void setupTwoBodyElements();
    void setupOverlapMatrix();
    void diagonalizeOverlapMatrix();

    void printInitialInfo();
    void printIterationInfo(int iteration);
    void printFinalInfo();

    virtual void    setup() = 0;
    virtual void    selfConsistentFieldIteration() = 0;
    virtual void    computeHartreeFockEnergy() = 0;
    virtual void    storeEnergy() = 0;
    virtual double  convergenceTest() = 0;

public:
    HartreeFock(class System* system);
    arma::mat   m_overlapMatrix;
    arma::mat   m_transformationMatrix;

    double solve(double convergenceCriterion=1e-14, int maximumIterations=50);
    double solveSilently(double convergenceCriterion=1e-14, int maximumIterations=50);
};
