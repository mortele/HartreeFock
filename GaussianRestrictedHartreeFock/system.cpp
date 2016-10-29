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

double System::oneBodyElements(int i, int j) {
    ContractedGaussian* contracted1 = m_basis.at(i);
    ContractedGaussian* contracted2 = m_basis.at(j);

    double kinetic = m_integrator.kineticIntegral(contracted1, contracted2);
    double nucleusCoulombInteraction = 0;
    for (Atom* atom : m_atoms) {
        nucleusCoulombInteraction -= atom->getCharge() *
                                     m_integrator.electronNucleusIntegral(contracted1,
                                                                          contracted2,
                                                                          atom->getPosition());
    }
    return kinetic + nucleusCoulombInteraction;
}

double System::twoBodyElements(int i, int j, int k, int l) {
    return electronElectronIntegral(i,j,k,l);
}

void System::setupBasis() {
    int numberOfBasisFunctions = 0;
    for (unsigned int i = 0; i < m_atoms.size(); i++) {
        numberOfBasisFunctions += m_atoms.at(i)->getNumberOfOrbitals();
    }
    m_basis.clear();
    m_basis.reserve(numberOfBasisFunctions);
    for (unsigned int i = 0; i < m_atoms.size(); i++) {
        for (unsigned int j = 0; j < m_atoms.at(i)->getContractedGaussians().size(); j++) {
            m_basis.push_back(m_atoms.at(i)->getContractedGaussians().at(j));
        }
    }
    m_numberOfBasisFunctions = m_basis.size();
}
