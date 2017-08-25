#pragma once
#include "Integrators/overlapintegrator.h"
#include <armadillo>

class IntegralTester {
private:
    double m_tollerance = 1e-13;
    arma::mat m_exactOverlap = arma::zeros(10,10);
    OverlapIntegrator* m_overlapIntegrator;

    void setupExactOverlapMatrix();

public:
    IntegralTester();
    bool runAllTests();
    bool runOverlapTests();

};
