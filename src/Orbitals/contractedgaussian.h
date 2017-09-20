#pragma once
#include <vector>
#include <armadillo>
#include "Orbitals/gaussianprimitive.h"

class ContractedGaussian {
private:
    double      m_coefficient           = 1;
    int         m_numberOfPrimitives    = 0;
    arma::vec   m_nucleusPosition       = arma::zeros<arma::vec>(3);
    std::vector<double>             m_coefficients;
    std::vector<GaussianPrimitive*> m_primitives;

public:
    void setCoefficient(double coefficient) { m_coefficient = coefficient; }
    double getCoefficients(int i) { return m_coefficients.at(i); }
    double evaluate(arma::vec &r);
    double evaluate(double x, double y, double z);
    double operator()(double x, double y, double z);
    void createNewPrimitive(int i, int j, int k, double a, double coefficient=1);
    void addPrimitive(GaussianPrimitive* primitive, double coefficient=1);
    std::vector<GaussianPrimitive*> getPrimitives() const { return m_primitives; }
    arma::vec getNucleusPosition() { return m_nucleusPosition; }
    friend std::ostream& operator<<(std::ostream& stream, const ContractedGaussian& contracted);
};
