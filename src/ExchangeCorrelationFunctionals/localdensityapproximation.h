#pragma once
#include "ExchangeCorrelationFunctionals/exchangecorrelationfunctional.h"
#include <armadillo>
#include <cmath>

class LocalDensityApproximation : public ExchangeCorrelationFunctional {
private:
    arma::mat* m_densityMatrix;

public:
    LocalDensityApproximation(class System* system, arma::mat* densityMatrix);
    double evaluateEnergy(double rho);
    double evaluatePotential(double rho);

    //double X(double x) { return x*x + m_b*x + m_c; }
    double  epsilonC(double rs, double rho);
    double dEpsilonC(double rs, double rho);
    double  epsilonX(double rs, double rho);
    double dEpsilonX(double rs, double rho);
    inline constexpr double Cx();
};

constexpr double LocalDensityApproximation::Cx() {
    // (3/4) (3/pi)^(1/3)
    return 0.73855876638202240;
}
