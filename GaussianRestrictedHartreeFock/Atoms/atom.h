#pragma once
#include "Orbitals/contractedgaussian.h"
#include <vector>
#include <armadillo>
#include <string>


class Atom {
    friend class BasisSetParser;

protected:
    int         m_numberOfElectrons;
    int         m_numberOfOrbitals;
    double      m_charge;
    arma::vec   m_position;
    std::string                         m_info = "Unknown atom type";
    std::vector<ContractedGaussian*>    m_contractedGaussians;

public:
    Atom(arma::vec position, int numberOfElectrons, double charge);
    Atom(arma::vec position, int numberOfOrbitals, int numberOfElectrons, double charge);

    void setInfo(std::string info);
    void setNumberOfOrbitals(int numberOfOrbitals);

    void        setNumberOfElectrons(int n) { m_numberOfElectrons = n; }
    std::string getInfo()               { return m_info; }
    arma::vec   getPosition()           { return m_position; }
    int         getNumberOfElectrons()  { return m_numberOfElectrons; }
    int         getNumberOfOrbitals()   { return m_numberOfOrbitals; }
    double      getCharge()             { return m_charge; }
    std::vector<ContractedGaussian*> getContractedGaussians() { return m_contractedGaussians; }

    friend std::ostream& operator<<(std::ostream& stream, Atom& atom);
};

