#pragma once
#include <vector>
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
};
