#pragma once
#include <numgrid.h>
#include <armadillo>



class NumericalIntegrator {
private:
    bool                m_gridGenerated = false;
    class   System*     m_system;
    class   Grid*       m_grid;
    class ExchangeCorrelationFunctional* m_xcFunctional;
    context_t*  m_context;

public:
    NumericalIntegrator(class System* system);
    int generateBeckeGrid();
    double integrateDensity(const arma::mat& densityMatrix);
    void setFunctional(class ExchangeCorrelationFunctional* functional) { m_xcFunctional = functional; }

    double testIntegral(const arma::mat& densityMatrix);
};

