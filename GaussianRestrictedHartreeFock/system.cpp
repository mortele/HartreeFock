#include "system.h"

System::System(int numberOfAtoms) {
    m_atoms.reserve(numberOfAtoms);
    m_numberOfAtoms = numberOfAtoms;
    m_integrator = ContractedIntegrator();
}

void System::addAtom(Atom* atom) {
    m_atoms.push_back(atom);
    setupBasis();
}

double System::overlapIntegral(int i, int j) {
    return m_integrator.overlapIntegral(m_basis.at(i), m_basis.at(j));
}

double System::kineticIntegral(int i, int j) {
    return m_integrator.kineticIntegral(m_basis.at(i), m_basis.at(j));
}

double System::electronNucleusIntegral(int i, int j, arma::vec nucleusPosition) {
    return m_integrator.electronNucleusIntegral(m_basis.at(i), m_basis.at(j), nucleusPosition);
}

double System::electronElectronIntegral(int i, int j, int k, int l) {
    return m_integrator.electronElectronIntegral(m_basis.at(i), m_basis.at(j), m_basis.at(k), m_basis.at(l));
}

void System::setupBasis() {
    int numberOfBasisFunctions = 0;
    for (int i = 0; i < m_atoms.size(); i++) {
        numberOfBasisFunctions += m_atoms.at(i)->getNumberOfOrbitals();
    }
    m_basis.clear();
    m_basis.reserve(numberOfBasisFunctions);
    for (int i = 0; i < m_atoms.size(); i++) {
        for (int j = 0; j < m_atoms.at(i)->getContractedGaussians().size(); j++) {
            m_basis.push_back(m_atoms.at(i)->getContractedGaussians().at(j));
        }
    }
    m_numberOfBasisFunctions = m_basis.size();
    arma::vec x {0,1,2};
    int i=0; int j=0;
}
