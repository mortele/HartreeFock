#pragma once
#include "Solvers/hartreefock.h"
#include "system.h"
#include "Orbitals/contractedgaussian.h"
#include <vector>
#include <armadillo>
#include <string>


class RestrictedDFT : public HartreeFock {
    friend class NumericalIntegrator;
    friend class ExchangeCorrelationFunctional;

public:
    bool        m_smoothing         = true;
    double      m_smoothingFactor   = 0.5;
    //double      m_xcEnergy;
    arma::vec   m_epsilon;
    arma::vec   m_epsilonOld;
    arma::mat   m_fockMatrix;
    arma::mat   m_fockMatrixTilde;
    arma::mat   m_coefficientMatrixTilde;
    arma::mat   m_coefficientMatrix;
    arma::mat   m_densityMatrix;
    arma::mat   m_xcMatrix;
    class NumericalIntegrator*              m_numericalIntegrator;
    class ExchangeCorrelationFunctional*    m_xcFunctional;

    double computeEnergyX();

    void setup();
    void computeXcMatrix();
    void computeFockMatrix();
    void computeDensityMatrix();
    void diagonalizeFockMatrix();
    void normalizeCoefficientMatrix();
    void selfConsistentFieldIteration();
    void computeHartreeFockEnergy();
    void storeEnergy();
    double twoBodyMatrixElements(int,int,int,int);
    double twoBodyMatrixElementsAntiSymmetric(int,int,int,int);
    double Vxc(int,int);
    double Exc(int,int);
    double convergenceTest();


public:
    RestrictedDFT(System* system);
    void setFunctional(std::string name);
};

