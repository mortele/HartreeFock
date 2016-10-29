#pragma once
#include "Orbitals/contractedgaussian.h"
#include <vector>
#include <armadillo>
#include <string>


class Atom {
protected:
    int         m_numberOfElectrons;
    int         m_numberOfOrbitals;
    double      m_charge;
    arma::vec   m_position;
    std::vector<ContractedGaussian*>  m_contractedGaussians;

public:
    Atom(arma::vec position, int numberOfOrbitals, int numberOfElectrons, double charge);

    std::vector<ContractedGaussian*> getContractedGaussians() { return m_contractedGaussians; }
    arma::vec getPosition()     { return m_position; }
    int getNumberOfElectrons()  { return m_numberOfElectrons; }
    int getNumberOfOrbitals()   { return m_numberOfOrbitals; }
    double getCharge()          { return m_charge; }
    virtual std::string getInfo();
};

