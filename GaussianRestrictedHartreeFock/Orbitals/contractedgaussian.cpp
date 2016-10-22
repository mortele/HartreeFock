#include "Orbitals/contractedgaussian.h"

using arma::vec;
using std::vector;

ContractedGaussian::ContractedGaussian(vec nucleusPosition) :
        m_numberOfPrimitives(0),
        m_nucleusPosition(nucleusPosition) {
}

double ContractedGaussian::evaluate(vec& r) {
    double functionValue = 0;
    for (int i=0; i<m_numberOfPrimitives; i++) {
        GaussianPrimitive*  primitive   = &m_primitives.at(i);
        double*             coefficient = &m_coefficients.at(i);
        functionValue += *coefficient * primitive->evaluate(r);
    }
    return functionValue;
}

void ContractedGaussian::createNewPrimitive(int     i,
                                            int     j,
                                            int     k,
                                            double  a,
                                            double  coefficient) {
    GaussianPrimitive* p = new GaussianPrimitive(i,j,k,a,m_nucleusPosition);
    m_primitives.push_back(*p);
    m_coefficients.push_back(coefficient);
    p->setCoefficient(coefficient);
    m_numberOfPrimitives += 1;
}

void ContractedGaussian::addPrimitive(GaussianPrimitive primitive,
                                      double            coefficient) {
    m_primitives.push_back(primitive);
    m_coefficients.push_back(coefficient);
    primitive.setCoefficient(coefficient);
    m_numberOfPrimitives += 1;
}
