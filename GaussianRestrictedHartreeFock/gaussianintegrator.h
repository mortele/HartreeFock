#pragma once
#include <vector>
#include <contractedgaussian.h>


class GaussianIntegrator {
private:
    std::vector<ContractedGaussian> m_contracteds;

public:
    GaussianIntegrator(std::vector<ContractedGaussian> contracteds);
};
