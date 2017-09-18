#pragma once
#include <numgrid.h>
#include <armadillo>



class NumericalIntegrator {
private:
    class   System*     m_system;
    class   Grid*       m_grid;
    context_t*  m_context;

public:
    NumericalIntegrator(class System* system);
    int generateBeckeGrid();
    double integrateDensity(const arma::mat& densityMatrix);

    double testIntegral(const arma::mat& densityMatrix);
};

