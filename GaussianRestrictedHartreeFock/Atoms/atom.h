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

    ContractedGaussian* insertNewContracted();
    ContractedGaussian* create_S1(double a, double c);
    ContractedGaussian* create_S2(double a1, double a2, double c1, double c2);
    ContractedGaussian* create_S3(double a1, double a2, double a3, double c1, double c2, double c3);
    void create_P1(double a, double c);

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

