#include "system.h"

System::System(int numberOfAtoms) {
    m_atoms.reserve(numberOfAtoms);
    m_numberOfAtoms = 0;
    m_integrator = ContractedIntegrator();
}

void System::addAtom(Atom* atom) {
    m_atoms.push_back(atom);
    m_numberOfAtoms += 1;
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

double System::nucleusNucleusInteractionEnergy() {
    double nucleusNucleusInteraction = 0;
    if (m_numberOfAtoms > 1) {
        for (int i = 0;     i < m_numberOfAtoms; i++)
        for (int j = i + 1; j < m_numberOfAtoms; j++) {
            Atom* atom1 = m_atoms.at(i);
            Atom* atom2 = m_atoms.at(j);
            arma::vec   r       = atom1->getPosition() - atom2->getPosition();
            double      rNorm   = arma::norm(r);
            nucleusNucleusInteraction += atom1->getCharge() *
                                         atom2->getCharge() *
                                         1.0 / rNorm;
        }
    }
    return nucleusNucleusInteraction;
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
