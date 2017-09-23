#pragma once
#include "Atoms/atom.h"
#include "Orbitals/contractedgaussian.h"
#include "Integrators/contractedintegrator.h"
#include <vector>
#include <string>

class System {
    friend class HartreeFock;
    friend class NumericalIntegrator;

public:
    int     m_numberOfAtoms;
    int     m_numberOfBasisFunctions;
    int     m_numberOfElectrons;
    std::vector<Atom*>                  m_atoms;
    std::vector<ContractedGaussian*>    m_basis;
    ContractedIntegrator                m_integrator;
    class HartreeFock*                  m_solver;
    void setupBasis();

protected:
    void setSolver(class HartreeFock* solver);
    class HartreeFock* getSolver() { return m_solver; }

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
