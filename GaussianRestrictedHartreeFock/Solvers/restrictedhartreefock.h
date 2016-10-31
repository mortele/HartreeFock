#pragma once
#include "system.h"
#include "Orbitals/contractedgaussian.h"
#include <vector>
#include <armadillo>

class RestrictedHartreeFock {
private:
    int         m_numberOfBasisFunctions;
    int         m_numberOfElectrons;
    int         m_maximumIterations;
    int         m_iterationsUsed;
    bool        m_reachedSelfConsistency;
    bool        m_silent = false;
    double      m_convergenceCriterion;
    double      m_hartreeFockEnergy;
    double      m_nucleusNucleusInteractionEnergy;
    double      m_convergenceTest;
    System*     m_system;
    arma::vec   m_epsilon;
    arma::vec   m_epsilonOld;
    arma::mat   m_fockMatrix;
    arma::mat   m_fockMatrixTilde;
    arma::mat   m_coefficientMatrix;
    arma::mat   m_coefficientMatrixTilde;
    arma::mat   m_densityMatrix;
    arma::mat   m_overlapMatrix;
    arma::mat   m_transformationMatrix;
    arma::mat   m_oneBodyMatrixElements;
    arma::field<arma::mat>           m_twoBodyMatrixElements;
    std::vector<ContractedGaussian*> m_basis;

    void setup();
    void computeFockMatrix();
    void computeDensityMatrix();
    void setupOverlapMatrix();
    void setupTwoBodyMatrixElements();
    void setupOneBodyMatrixElements();
    void diagonalizeFockMatrix();
    void diagonalizeOverlapMatrix();
    void normalizeCoefficientMatrix();
    void selfConsistentFieldIteration();
    void computeHartreeFockEnergy();
    void printInitialInfo();
    void printIterationInfo(int iteration);
    void printFinalInfo();

public:
    RestrictedHartreeFock(System* system);
    double solve(double convergenceCriterion=1e-14, int maximumIterations=50);
    double solveSilently(double convergenceCriterion=1e-14, int maximumIterations=50);
};

