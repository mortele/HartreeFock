#pragma once
#include <numgrid.h>
#include <armadillo>



class NumericalIntegrator {
private:
    bool    m_xcFunctionalSet = false;
    bool    m_gridGenerated   = false;
    class   System*                         m_system;
    class   Grid*                           m_grid;
    class   ExchangeCorrelationFunctional*  m_xcFunctional;
    context_t*  m_context;
    arma::mat*  m_densityMatrix;

public:
    NumericalIntegrator(class System* system, arma::mat* densityMatrix);
    int generateBeckeGrid();
    double testIntegral();
    double integrateExchangeCorrelationPotential(double Ppq, class ContractedGaussian* Gp, class ContractedGaussian* Gq, int p, int q);
    double integrateExchangeCorrelationPotential(int,int);
    double integrateExchangeCorrelationEnergy();
    void setFunctional(class ExchangeCorrelationFunctional* functional);

    double testIntegral(const arma::mat& densityMatrix);
};

