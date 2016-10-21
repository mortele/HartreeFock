#pragma once

class BoysFunction {
private:
    int m_taylorExpansionTerms      = 6;
    int m_numberOfIntegrationPoints = 1e5;

    double integrand(double x, double t, double n);

public:
    BoysFunction();
    double directIntegration(double x, double n);


};

