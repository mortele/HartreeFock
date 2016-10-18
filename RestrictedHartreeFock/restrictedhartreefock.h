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
    void setOverLapMatrix(arma::mat overLap);
    bool setIntegralTable(std::string fileName);

private:

    double m_convergencePrecision = 1e-15;
    int m_MAX_ITERS = 50;
    int m_nr_OfIters = 0;
    int m_nrOfSpatialOrbitals;
    int m_nrOfSpinOrbitals;
    int m_nrOfParticles;
    int m_nrOfOccupiedOrbitals;
    double m_HartreeFockEnergy = 0;
    bool m_overlap = false;
    bool m_reachedSelfConsistency = false;


    arma::mat m_DensityMatrix;
    arma::mat m_U;
    arma::mat m_FockMatrix;
    arma::vec m_eps;
    arma::vec m_eps_old;
    arma::mat m_oneBodyElements;
    arma::mat m_S;

    IntegralTable* m_integralTable;

    //vector oneBodyElements
    //vector something twoBodyElements

    double getOneBodyMatrixElement(int p, int q);
    double QRPS_AntiSym(int q, int r, int p, int s);

    double computeHartreeFockEnergy();

    void computeFockMatrix();
    void computeDensityMatrix();
    void diagonalizeFockMatrix();
    void printInfo();

    //Unit test functions?
    bool isFockMatrixHermitian();
};

