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

ContractedGaussian*Atom::create_S4(double a1, double a2, double a3, double a4, double c1, double c2, double c3, double c4) {
    ContractedGaussian* contracted = create_S3(a1,a2,a3,c1,c2,c3);
    contracted->addPrimitive(new GaussianPrimitive(0,0,0,a4,m_position,c4), c4);
    return contracted;
}

ContractedGaussian*Atom::create_S5(double a1, double a2, double a3, double a4, double a5, double c1, double c2, double c3, double c4, double c5) {
    ContractedGaussian* contracted = create_S4(a1,a2,a3,a4,c1,c2,c3,c4);
    contracted->addPrimitive(new GaussianPrimitive(0,0,0,a5,m_position,c5), c5);
    return contracted;
}

ContractedGaussian*Atom::create_S6(double a1, double a2, double a3, double a4, double a5, double a6, double c1, double c2, double c3, double c4, double c5, double c6) {
    ContractedGaussian* contracted = create_S5(a1,a2,a3,a4,a5,c1,c2,c3,c4,c5);
    contracted->addPrimitive(new GaussianPrimitive(0,0,0,a6,m_position,c6), c6);
    return contracted;
}


ContractedGaussian* Atom::create_P3(double a1, double a2, double a3, double c1, double c2, double c3) {
    ContractedGaussian* contracted;
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a3,m_position,c3), c3);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a3,m_position,c3), c3);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a3,m_position,c3), c3);
    return contracted;
}

ContractedGaussian*Atom::create_P6(double a1, double a2, double a3, double a4, double a5, double a6, double c1, double c2, double c3, double c4, double c5, double c6) {
    ContractedGaussian* contracted;
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a3,m_position,c3), c3);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a4,m_position,c4), c4);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a5,m_position,c5), c5);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a6,m_position,c6), c6);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a3,m_position,c3), c3);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a4,m_position,c4), c4);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a5,m_position,c5), c5);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a6,m_position,c6), c6);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a2,m_position,c2), c2);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a3,m_position,c3), c3);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a4,m_position,c4), c4);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a5,m_position,c5), c5);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a6,m_position,c6), c6);
    return contracted;
}

ContractedGaussian* Atom::create_D1(double a, double c) {
    ContractedGaussian* contracted;
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(2,0,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,2,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,0,2,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,1,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,0,1,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,1,1,a,m_position,c), c);
    return contracted;
}

ContractedGaussian*Atom::create_F1(double a, double c) {
    ContractedGaussian* contracted;
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(3,0,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,3,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,0,3,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(2,1,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(2,0,1,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,2,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,2,1,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,2,0,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,0,2,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,1,2,a,m_position,c), c);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,1,1,a,m_position,c), c);
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

ContractedGaussian* Atom::create_P2(double a1, double a2, double c1, double c2) {
    ContractedGaussian* contracted;
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(1,0,0,a2,m_position,c2), c2);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(0,1,0,a2,m_position,c2), c2);
    contracted = insertNewContracted();
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a1,m_position,c1), c1);
    contracted->addPrimitive(new GaussianPrimitive(0,0,1,a2,m_position,c2), c2);
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
