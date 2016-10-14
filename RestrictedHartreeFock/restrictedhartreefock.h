#pragma once
#include <armadillo>
#include <string>
#include "../Integrator/integraltable.h"

class RestrictedHartreeFock {
public:
    RestrictedHartreeFock(int nrOfParticles, int nrOfSpinOrbitals);
    void computeSolutionBySCF();
    void setAnalyticOneBodyElements(arma::vec oneBodyElements);
    void setAnalyticOneBodyElements(arma::mat oneBodyElements);
    bool setIntegralTable(std::string fileName);

private:
    int m_nrOfSpatialOrbitals;
    int m_nrOfSpinOrbitals;
    int m_nrOfParticles;
    int m_nrOfOccupiedOrbitals;
    double m_HartreeFockEnergy = 0;

    arma::mat m_DensityMatrix;
    arma::mat m_U;
    arma::mat m_FockMatrix;
    arma::vec m_eps;
    arma::mat m_oneBodyElements;

    IntegralTable* m_integralTable;

    //vector oneBodyElements
    //vector something twoBodyElements

    double getOneBodyMatrixElement(int p, int q);
    double getTwoBodyMatrixElement(int p, int q, int r, int s);

    double computeHartreeFockEnergy();

    void computeFockMatrix();
    void computeDensityMatrix();
    void diagonalizeFockMatrix();

    //Unit test functions?
    bool isFockMatrixHermitian();
};

