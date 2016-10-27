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
    ContractedGaussian* contracted1 = m_basis.at(i);
    ContractedGaussian* contracted2 = m_basis.at(j);
    return m_integrator.overlapIntegral(contracted1, contracted2);
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
            ContractedGaussian* contracted = m_atoms.at(i)->getContractedGaussians().at(j);
            m_basis.push_back(contracted);
        }
    }
    m_numberOfBasisFunctions = m_basis.size();
}
