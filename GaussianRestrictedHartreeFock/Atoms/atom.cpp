#include "atom.h"


Atom::Atom(arma::vec    position,
           int          numberOfOrbitals,
           int          numberOfElectrons,
           double       charge) {

    m_charge = charge;
    m_position = position;
    m_contractedGaussians.clear();
    m_contractedGaussians.reserve(numberOfOrbitals);
    m_numberOfOrbitals  = numberOfOrbitals;
    m_numberOfElectrons = numberOfElectrons;
}

std::string Atom::getInfo() {
    return "Unknown atom type.";
}
