#pragma once
#include <armadillo>

using arma::vec;

class GaussianPrimitive {
    friend class ContractedGaussian;
private:
    int      m_xExponent        = 0;
    int      m_yExponent        = 0;
    int      m_zExponent        = 0;
    double   m_exponent         = 0;
    vec      m_nucleusPosition  = arma::zeros<vec>(3);
    double   m_coefficient      = 0;

public:
    GaussianPrimitive(int i, int j, int k, double a, vec nucleusPosition);
    void setCoefficient(double coefficient);
    double evaluate(vec r);
};

