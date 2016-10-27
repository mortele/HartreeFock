#include "atom.h"


Atom::Atom(arma::vec position, int numberOfOrbitals, int numberOfElectrons) {
    m_position = position;
    m_contractedGaussians.clear();
    m_contractedGaussians.reserve(numberOfOrbitals);
    m_numberOfOrbitals  = numberOfOrbitals;
    m_numberOfElectrons = numberOfElectrons;
}
