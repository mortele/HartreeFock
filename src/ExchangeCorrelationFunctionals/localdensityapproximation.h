#pragma once
#include "ExchangeCorrelationFunctionals/exchangecorrelationfunctional.h"
#include <armadillo>
#include <cmath>

class LocalDensityApproximation : public ExchangeCorrelationFunctional {
private:
    double m_A   ;
    double m_x0  ;
    double m_b   ;
    double m_c   ;
    double m_Q   ;
    double m_Xx0 ;
    double m_pi  ;
    arma::mat* m_densityMatrix;

public:
    LocalDensityApproximation(class System* system, arma::mat* densityMatrix);
    double evaluate(double x, double y, double z, int p, int q);

    double X(double x) { return x*x + m_b*x + m_c; }
    double epsilonC(double rs);
    double dEpsilonC(double rs);
    double epsilonX(double rs);
    double dEpsilonX(double rs);
    inline constexpr double Cx();
};

constexpr double LocalDensityApproximation::Cx() {
    // (3/4) (3/pi)^(1/3)
    return 0.73855876638202240;
}
