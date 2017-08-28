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

    // S orbitals.
    ContractedGaussian* create_S1(double a, double c);
    ContractedGaussian* create_S2(double a1, double a2, double c1, double c2);
    ContractedGaussian* create_S3(double a1, double a2, double a3, double c1, double c2, double c3);
    ContractedGaussian* create_S4(double a1, double a2, double a3, double a4, double c1, double c2, double c3, double c4);
    ContractedGaussian* create_S5(double a1, double a2, double a3, double a4, double a5, double c1, double c2, double c3, double c4, double c5);
    ContractedGaussian* create_S6(double a1, double a2, double a3, double a4, double a5, double a6, double c1, double c2, double c3, double c4, double c5, double c6);

    // P orbitals.
    ContractedGaussian* create_P1(double a, double c);
    ContractedGaussian* create_P2(double a1, double a2, double c1, double c2);
    ContractedGaussian* create_P3(double a1, double a2, double a3, double c1, double c2, double c3);

    // D orbitals.
    ContractedGaussian* create_D1(double a, double c);

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

#include "Atoms/hydrogen.h"
#include "Atoms/helium.h"
#include "Atoms/lithium.h"
#include "Atoms/beryllium.h"
#include "Atoms/boron.h"
#include "Atoms/carbon.h"
#include "Atoms/nitrogen.h"
#include "Atoms/oxygen.h"
#include "Atoms/fluorine.h"
#include "Atoms/neon.h"
