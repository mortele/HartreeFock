#include "atom.h"
#include <iomanip>

using std::setprecision;

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

void Atom::setInfo(std::string info) {
    m_info = info;
}

std::ostream& operator<<(std::ostream& stream, Atom& atom) {
    stream << atom.getInfo() << " at (" << setprecision(4)
                                        << atom.m_position(0) << ", "
                                        << atom.m_position(1) << ", "
                                        << atom.m_position(2) << ")";
    return stream;
}
