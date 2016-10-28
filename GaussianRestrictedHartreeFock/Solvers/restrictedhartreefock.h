#pragma once
#include "system.h"
#include "Orbitals/contractedgaussian.h"
#include <vector>
#include <armadillo>

class RestrictedHartreeFock {
private:
    int                                 m_numberOfBasisFunctions;
    int                                 m_numberOfElectrons;
    System*                             m_system;
    std::vector<ContractedGaussian*>    m_basis;
    arma::vec                           m_epsilon;
    arma::vec                           m_epsilonOld;
    arma::mat                           m_fockMatrix;
    arma::mat                           m_U;
    arma::mat                           m_densityMatrix;

    void computeFockMatrix();

public:
    RestrictedHartreeFock(System* system);
    double solve(double convergenceCriteria=1.0e-8);
};

