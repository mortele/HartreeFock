#include "Integrators/contractedintegrator.h"

using std::cout;
using std::endl;
using arma::vec;
using arma::zeros;



double ContractedIntegrator::overlapIntegral(ContractedGaussian* contracted1,
                                             ContractedGaussian* contracted2) {

    double integral = 0;
    for (unsigned int i = 0; i < contracted1->getPrimitives().size(); i++) {
        for (unsigned int j = 0; j < contracted2->getPrimitives().size(); j++) {

            integral += contracted1->getPrimitives().at(i)->getConstantTerm() *
                        contracted2->getPrimitives().at(j)->getConstantTerm() *
                        m_overlapIntegrator.computeIntegral(contracted1->getPrimitives().at(i),
                                                            contracted2->getPrimitives().at(j));
        }
    }
    return integral;
}

double ContractedIntegrator::kineticIntegral(ContractedGaussian* contracted1,
                                             ContractedGaussian* contracted2) {

    double integral = 0;
    for (unsigned int i = 0; i < contracted1->getPrimitives().size(); i++) {
        for (unsigned int j = 0; j < contracted2->getPrimitives().size(); j++) {

            integral += contracted1->getPrimitives().at(i)->getConstantTerm() *
                        contracted2->getPrimitives().at(j)->getConstantTerm() *
                        m_kineticIntegrator.computeIntegral(contracted1->getPrimitives().at(i),
                                                            contracted2->getPrimitives().at(j));
        }
    }
    return integral;
}

double ContractedIntegrator::electronNucleusIntegral(ContractedGaussian* contracted1,
                                                     ContractedGaussian* contracted2,
                                                     arma::vec nucleusPosition) {
    m_electronNucleusIntegrator.setNucleusPosition(nucleusPosition);
    double integral = 0;
    for (unsigned int i = 0; i < contracted1->getPrimitives().size(); i++) {
        for (unsigned int j = 0; j < contracted2->getPrimitives().size(); j++) {

            integral += contracted1->getPrimitives().at(i)->getConstantTerm() *
                        contracted2->getPrimitives().at(j)->getConstantTerm() *
                        m_electronNucleusIntegrator.computeIntegral(contracted1->getPrimitives().at(i),
                                                                    contracted2->getPrimitives().at(j));
        }
    }
    return integral;
}

double ContractedIntegrator::electronElectronIntegral(ContractedGaussian* contracted1,
                                                      ContractedGaussian* contracted2,
                                                      ContractedGaussian* contracted3,
                                                      ContractedGaussian* contracted4) {
    double integral = 0;
    for (unsigned int i = 0; i < contracted1->getPrimitives().size(); i++)
    for (unsigned int j = 0; j < contracted2->getPrimitives().size(); j++)
    for (unsigned int k = 0; k < contracted3->getPrimitives().size(); k++)
    for (unsigned int l = 0; l < contracted4->getPrimitives().size(); l++) {
        integral += contracted1->getPrimitives().at(i)->getConstantTerm() *
                    contracted2->getPrimitives().at(j)->getConstantTerm() *
                    contracted3->getPrimitives().at(k)->getConstantTerm() *
                    contracted4->getPrimitives().at(l)->getConstantTerm() *
                    m_electronElectronIntegrator.computeIntegral(contracted1->getPrimitives().at(i),
                                                                 contracted2->getPrimitives().at(j),
                                                                 contracted3->getPrimitives().at(k),
                                                                 contracted4->getPrimitives().at(l));
    }
    return integral;
}









