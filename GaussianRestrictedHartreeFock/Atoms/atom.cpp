#include "atom.h"
#include <iomanip>

using std::setprecision;


ContractedGaussian* Atom::insertNewContracted() {
    ContractedGaussian* contracted = new ContractedGaussian();
    m_contractedGaussians.push_back(contracted);
    return contracted;
}

ContractedGaussian* Atom::create_S1(double a, double c) {
    ContractedGaussian* contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,0,0,a,m_position,c), c);
    return contracted;
}

ContractedGaussian* Atom::create_S2(double a1, double a2, double c1, double c2) {
    ContractedGaussian* contracted = create_S1(a1,c1);
    contracted->addPrimitive(new GaussianPrimitive(0,0,0,a2,m_position,c2), c2);
    return contracted;
}

ContractedGaussian* Atom::create_S3(double a1, double a2, double a3, double c1, double c2, double c3) {
    ContractedGaussian* contracted = create_S2(a1,a2,c1,c2);
    contracted->addPrimitive(new GaussianPrimitive(0,0,0,a3,m_position,c3), c3);
    return contracted;
}

ContractedGaussian* Atom::create_P2(double a1, double a2, double c1, double c2) {
    ContractedGaussian* contracted = create_P1(a1, c1);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a2,m_position,c2), c2);
    return contracted;
}

ContractedGaussian* Atom::create_P1(double a, double c) {
    ContractedGaussian* contracted;
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a,m_position,c), c);
    return contracted;
}

Atom::Atom(arma::vec position, int numberOfElectrons, double charge) {
    m_position = position;
    m_numberOfElectrons = numberOfElectrons;
    m_charge = charge;
}

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

void Atom::setNumberOfOrbitals(int numberOfOrbitals) {
    m_numberOfOrbitals = numberOfOrbitals;
    m_contractedGaussians.clear();
    m_contractedGaussians.reserve(numberOfOrbitals);
}

std::ostream& operator<<(std::ostream& stream, Atom& atom) {
    stream << atom.getInfo() << " at (" << setprecision(4)
                                        << atom.m_position(0) << ", "
                                        << atom.m_position(1) << ", "
                                        << atom.m_position(2) << ")";
    return stream;
}
