#pragma once
#include "ExchangeCorrelationFunctionals/exchangecorrelationfunctional.h"
#include <armadillo>
#include <cmath>

class LocalDensityApproximation : public ExchangeCorrelationFunctional {
private:
    arma::mat* m_densityMatrix;

public:
    LocalDensityApproximation(class System* system, arma::mat* densityMatrix);
    double evaluate(double x, double y, double z, int p, int q);
    inline constexpr double Cx();
};

constexpr double LocalDensityApproximation::Cx() {
    // (3/4) (3/pi)^(1/3)
    return 0.738558766382022405884230032680836267782320120689011037061;
}
