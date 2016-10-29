#pragma once
#include <vector>
#include <armadillo>
#include "Orbitals/gaussianprimitive.h"

class ContractedGaussian {
private:
    int         m_numberOfPrimitives    = 0;
    arma::vec   m_nucleusPosition       = arma::zeros<arma::vec>(3);
    std::vector<double>             m_coefficients;
    std::vector<GaussianPrimitive*> m_primitives;

public:
    double evaluate(arma::vec &r);
    void createNewPrimitive(int i, int j, int k, double a, double coefficient=1);
    void addPrimitive(GaussianPrimitive* primitive, double coefficient=1);
    std::vector<GaussianPrimitive*> getPrimitives() const { return m_primitives; }
};
