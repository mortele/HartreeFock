#include "Orbitals/contractedgaussian.h"

using arma::vec;
using std::vector;
using std::cout;
using std::endl;



double ContractedGaussian::evaluate(vec& r) {
    double functionValue = 0;
    for (int i=0; i<m_numberOfPrimitives; i++) {
        double coefficient = m_coefficients.at(i);
        functionValue += coefficient * m_primitives.at(i)->evaluate(r);
    }
    return m_coefficient * functionValue;
}

double ContractedGaussian::evaluate(double x, double y, double z) {
    double functionValue = 0;
    for (int i=0; i<m_numberOfPrimitives; i++) {
        double coefficient = m_coefficients.at(i);
        functionValue += coefficient * m_primitives.at(i)->evaluate(x,y,z);
    }
    return m_coefficient * functionValue;
}

double ContractedGaussian::operator()(double x, double y, double z) {
    return evaluate(x,y,z);
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

std::ostream& operator<<(std::ostream& stream, const ContractedGaussian& contracted) {

    cout << "nucleus: " << "(" << contracted.m_nucleusPosition.at(0) << ","
                               << contracted.m_nucleusPosition.at(1) << ","
                               << contracted.m_nucleusPosition.at(2) << ")" << endl;

    int i = 0;
    for (GaussianPrimitive* primitive : contracted.m_primitives) {
        stream << "   * primitive" << i << ": " << *primitive << endl;
        i++;
    }
    return stream;
}
