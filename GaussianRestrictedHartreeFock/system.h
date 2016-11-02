#pragma once
#include "Atoms/atom.h"
#include "Orbitals/contractedgaussian.h"
#include "Integrators/contractedintegrator.h"
#include <vector>

class System {
private:
    int     m_numberOfAtoms;
    int     m_numberOfBasisFunctions;
    int     m_numberOfElectrons;
    std::vector<Atom*>                  m_atoms;
    std::vector<ContractedGaussian*>    m_basis;
    ContractedIntegrator                m_integrator;

    void setupBasis();

public:
    System(int numberOfAtoms=2);
    void addAtom(Atom* atom);
    std::vector<Atom*>& getAtoms() { return m_atoms; }
    std::vector<ContractedGaussian*> getBasis() { return m_basis; }
    int getNumberOfBasisFunctions() { return m_numberOfBasisFunctions; }
    int getNumberOfSpinUpElectrons();
    int getNumberOfSpinDownElectrons();

    double overlapIntegral(int i, int j);
    double kineticIntegral(int i, int j);
    double electronNucleusIntegral(int i, int j, arma::vec nucleusPosition);
    double electronElectronIntegral(int i, int j, int k, int l);

    double oneBodyElements(int i, int j);
    double twoBodyElements(int i, int j, int k, int l);
    double nucleusNucleusInteractionEnergy();
};
