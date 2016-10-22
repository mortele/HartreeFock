#pragma once
#include <vector>
#include "Orbitals/contractedgaussian.h"


class GaussianIntegrator {
private:
    std::vector<ContractedGaussian> m_contracteds;

public:
    GaussianIntegrator(std::vector<ContractedGaussian> contracteds);
};
