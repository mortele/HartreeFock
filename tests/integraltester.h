#pragma once
#include "Integrators/overlapintegrator.h"
#include "Integrators/kineticintegrator.h"
#include "Integrators/electronnucleusintegrator.h"
#include "Integrators/electronelectronintegrator.h"
#include "Orbitals/gaussianprimitive.h"
#include <vector>
#include <armadillo>

class IntegralTester {
private:
    double m_tollerance             = 1e-10;
    double m_tolleranceNumerical    = 1e-2;
    std::vector<GaussianPrimitive> m_primitives;

    arma::mat m_exactOverlap            = arma::zeros(10,10);
    arma::mat m_exactKinetic            = arma::zeros(10,10);
    arma::mat m_exactElectronNucleus    = arma::zeros(10,10);
    arma::mat m_exactElectronElectron   = arma::zeros(10);

    OverlapIntegrator*          m_overlapIntegrator;
    KineticIntegrator*          m_kineticIntegrator;
    ElectronNucleusIntegrator*  m_electronNucleusIntegrator;
    ElectronElectronIntegrator* m_electronElectronIntegrator;

    void setupExactOverlapMatrix();
    void setupExactKineticMatrix();
    void setupExactElectronNucleusMatrix();
    void setupExactElectronElectronVector();

public:
    IntegralTester();
    bool runAllTests(bool silent=false);
    bool runOverlapTests(bool silent=false);
    bool runKineticTests(bool silent=false);
    bool runElectronNucleusTests(bool silent=false);
    bool runElectronElectronTests(bool silent=false);
};
