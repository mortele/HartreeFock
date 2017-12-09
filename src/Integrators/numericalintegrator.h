#pragma once
#include <numgrid.h>
#include <armadillo>
#include <xc.h>



class NumericalIntegrator {
private:
    bool    m_xcFunctionalSet = false;
    bool    m_gridGenerated   = false;
    class   System*                         m_system;
    class   Grid*                           m_grid;
    class   ExchangeCorrelationFunctional*  m_xcFunctional;
    context_t*  m_context;
    arma::mat*  m_densityMatrix;
    xc_func_type m_c;
    xc_func_type m_x;



public:
    NumericalIntegrator(class System* system, arma::mat* densityMatrix);
    int generateBeckeGrid();
    double integrateDensity();
    double integrateExchangeCorrelationPotential(class ContractedGaussian* Gp, class ContractedGaussian* Gq, int p, int q);
    double integrateExchangeCorrelationPotential(int,int);
    double integrateExchangeCorrelationEnergy();
    void setFunctional(class ExchangeCorrelationFunctional* functional);

    double testIntegral(int p);
};

