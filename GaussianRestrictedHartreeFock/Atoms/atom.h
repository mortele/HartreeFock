#pragma once
#include "Orbitals/contractedgaussian.h"
#include <vector>
#include <armadillo>


class Atom {
protected:
    int                             m_numberOfElectrons;
    int                             m_numberOfOrbitals;
    arma::vec                       m_position;
    std::vector<ContractedGaussian> m_contractedGaussians;

public:
    Atom(arma::vec position, int numberOfOrbitals, int numberOfElectrons);

    std::vector<ContractedGaussian>& getContractedGaussians() { return m_contractedGaussians; }
    arma::vec getPosition() { return m_position; }
};

