#pragma once
#include <vector>
#include <armadillo>
#include "Orbitals/contractedgaussian.h"
#include "Orbitals/gaussianprimitive.h"
#include "Integrators/electronelectronintegrator.h"
#include "Integrators/electronnucleusintegrator.h"
#include "Integrators/kineticintegrator.h"
#include "Integrators/overlapintegrator.h"


class ContractedIntegrator {
private:
    OverlapIntegrator           m_overlapIntegrator;
    KineticIntegrator           m_kineticIntegrator;
    ElectronElectronIntegrator  m_electronElectronIntegrator;
    ElectronNucleusIntegrator   m_electronNucleusIntegrator;


public:
    double overlapIntegral(ContractedGaussian* contracted1,
                           ContractedGaussian* contracted2);
    double kineticIntegral(ContractedGaussian* contracted1,
                           ContractedGaussian* contracted2);
    double electronNucleusIntegral(ContractedGaussian* contracted1,
                                   ContractedGaussian* contracted2,
                                   arma::vec nucleusPosition);
    double electronElectronIntegral(ContractedGaussian* contracted1,
                                    ContractedGaussian* contracted2,
                                    ContractedGaussian* contracted3,
                                    ContractedGaussian* contracted4);
};
