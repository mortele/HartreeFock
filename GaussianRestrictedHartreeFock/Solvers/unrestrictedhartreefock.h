#pragma once
#include "system.h"
#include "Solvers/hartreefock.h"

class UnrestrictedHartreeFock : public HartreeFock {
private:
    int         m_numberOfSpinUpElectrons;
    int         m_numberOfSpinDownElectrons;
    arma::vec   m_epsilonUp;
    arma::vec   m_epsilonDown;
    arma::vec   m_epsilonOldUp;
    arma::vec   m_epsilonOldDown;
    arma::mat   m_fockMatrixUp;
    arma::mat   m_fockMatrixDown;
    arma::mat   m_fockMatrixTildeUp;
    arma::mat   m_fockMatrixTildeDown;
    arma::mat   m_coefficientMatrixUp;
    arma::mat   m_coefficientMatrixDown;
    arma::mat   m_coefficientMatrixTildeUp;
    arma::mat   m_coefficientMatrixTildeDown;
    arma::mat   m_densityMatrixUp;
    arma::mat   m_densityMatrixDown;

    void setup();
    void computeFockMatrices();
    void computeDensityMatrices();
    void diagonalizeFockMatrices();
    void normalizeCoefficientMatrices();
    void selfConsistentFieldIteration();
    void computeHartreeFockEnergy();
    void storeEnergy();
    double convergenceTest();

public:
    UnrestrictedHartreeFock(System* system);
};

