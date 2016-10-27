#include "Integrators/contractedintegrator.h"

using std::cout;
using std::endl;
using arma::vec;
using arma::zeros;

double ContractedIntegrator::overlapIntegral(ContractedGaussian* contracted1,
                                             ContractedGaussian* contracted2) {

    double integral = 0;
    for (int i = 0; i < contracted1->getPrimitives().size(); i++) {
        for (int j = 0; j < contracted2->getPrimitives().size(); j++) {
            integral += contracted1->getPrimitives().at(i)->getConstantTerm() *
                        contracted2->getPrimitives().at(j)->getConstantTerm() *
                        m_overlapIntegrator.computeIntegral(contracted1->getPrimitives().at(i),
                                                            contracted2->getPrimitives().at(j));
        }
    }
    return integral;
}
