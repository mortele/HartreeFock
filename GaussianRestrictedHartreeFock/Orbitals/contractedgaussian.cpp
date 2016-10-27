#include "Orbitals/contractedgaussian.h"

using arma::vec;
using std::vector;
using std::cout;
using std::endl;


//ContractedGaussian::ContractedGaussian(const ContractedGaussian& copy) {
//    cout << "copy constructor contracted" << endl;
//    m_numberOfPrimitives = copy.m_numberOfPrimitives;
//    m_nucleusPosition = arma::zeros<vec>(3);
//    for (int i=0; i<3; i++) {
//        m_nucleusPosition.at(i) = copy.m_nucleusPosition.at(i);
//    }
//    m_coefficients.clear();
//    m_coefficients.reserve(copy.m_coefficients.size());
//    for (int i=0; i<copy.m_coefficients.size(); i++) {
//        m_coefficients.at(i) = copy.m_coefficients.at(i);
//    }
//    m_primitives.clear();
//    m_primitives.reserve(copy.m_primitives.size());
//    for (int i=0; i<copy.m_primitives.size(); i++) {
//        m_primitives.at(i) = copy.m_primitives.at(i);
//    }
//}
//
//ContractedGaussian::ContractedGaussian(ContractedGaussian&& move) {
//    cout << "move constructor contracted" << endl;
//    m_numberOfPrimitives = move.m_numberOfPrimitives;
//    m_nucleusPosition = arma::zeros<vec>(3);
//    for (int i=0; i<3; i++) {
//        m_nucleusPosition.at(i) = move.m_nucleusPosition.at(i);
//    }
//    m_coefficients.clear();
//    m_coefficients.reserve(move.m_coefficients.size());
//    for (int i=0; i<move.m_coefficients.size(); i++) {
//        m_coefficients.at(i) = move.m_coefficients.at(i);
//    }
//    m_primitives.clear();
//    m_primitives.reserve(move.m_primitives.size());
//    for (int i=0; i<move.m_primitives.size(); i++) {
//        m_primitives.at(i) = move.m_primitives.at(i);
//    }
//}

double ContractedGaussian::evaluate(vec& r) {
    double functionValue = 0;
    for (int i=0; i<m_numberOfPrimitives; i++) {
        double coefficient = m_coefficients.at(i);
        functionValue += coefficient * m_primitives.at(i)->evaluate(r);
    }
    return functionValue;
}

void ContractedGaussian::createNewPrimitive(int     i,
                                            int     j,
                                            int     k,
                                            double  a,
                                            double  coefficient) {
    GaussianPrimitive* p = new GaussianPrimitive(i,j,k,a,m_nucleusPosition);
    m_primitives.push_back(p);
    m_coefficients.push_back(coefficient);
    p->setCoefficient(coefficient);
    m_numberOfPrimitives += 1;
}

void ContractedGaussian::addPrimitive(GaussianPrimitive* primitive,
                                      double             coefficient) {
    m_nucleusPosition = primitive->nucleusPosition();
    m_primitives.push_back(primitive);
    m_coefficients.push_back(coefficient);
    primitive->setCoefficient(coefficient);
    m_numberOfPrimitives += 1;
}
