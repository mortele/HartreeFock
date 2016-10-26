#pragma once

class BoysFunction {
private:
    int m_taylorExpansionTerms      = 6;
    int m_numberOfIntegrationPoints = 1e5;

    double integrand(double x, double t, double n);
    double directIntegration(double x, double n);

    // From http://link.springer.com/article/10.1007/s10910-005-9023-3
    double analyticalIncompleteGammaFunction(double x, double n);
    double analyticalCompleteGammaFunction(double x, double n);

public:
    BoysFunction();
    double compute(double x, double n);
};

