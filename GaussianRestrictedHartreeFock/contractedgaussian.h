#pragma once
#include <vector>
#include <armadillo>
#include "gaussianprimitive.h"

using std::vector;
using arma::vec;

class ContractedGaussian {
private:
    int                         m_numberOfPrimitives    = 0;
    vec                         m_nucleusPosition       = arma::zeros<vec>(3);
    vector<double>              m_coefficients;
    vector<GaussianPrimitive>   m_primitives;

public:
    ContractedGaussian(vec nucleusPosition);
    double evaluate(vec &r);
    void createNewPrimitive(int i, int j, int k, double a, double coefficient=1);
    void addPrimitive(GaussianPrimitive primitive, double coefficient=1);
};
