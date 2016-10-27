#pragma once
#include "Atoms/atom.h"
#include "Orbitals/contractedgaussian.h"
#include "Integrators/contractedintegrator.h"
#include <vector>

class System {
private:
    int                                 m_numberOfAtoms;
    int                                 m_numberOfBasisFunctions;
    std::vector<Atom*>                  m_atoms;
    std::vector<ContractedGaussian*>    m_basis;
    ContractedIntegrator                m_integrator;

    void setupBasis();

public:
    System(int numberOfAtoms);
    void addAtom(Atom* atom);
    std::vector<Atom*>& getAtoms() { return m_atoms; }
    std::vector<ContractedGaussian*> getBasis() { return m_basis; }
    int getNumberOfBasisFunctions() { return m_numberOfBasisFunctions; }

    double overlapIntegral(int i, int j);
};