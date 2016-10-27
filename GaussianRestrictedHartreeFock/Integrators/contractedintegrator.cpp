#include "Integrators/contractedintegrator.h"

using std::cout;
using std::endl;
using arma::vec;
using arma::zeros;

double ContractedIntegrator::overlapIntegral(ContractedGaussian* contracted1,
                                             ContractedGaussian* contracted2) {

    vec x = zeros<vec>(3);
    cout << contracted1->getPrimitives().at(0).evaluate(x) << endl;
    std::vector<GaussianPrimitive> primitives2 = contracted2->getPrimitives();
/*
    double integral = 0;
    for (int i = 0; i < contracted1.getPrimitives().size(); i++) {
        for (int j = 0; j < contracted2.getPrimitives().size(); j++) {
            auto primitive1 = contracted1.getPrimitives().at(i);
            auto primitive2 = contracted2.getPrimitives().at(j);
            integral += primitive1.getConstantTerm() *
                        primitive2.getConstantTerm() *
                        m_overlapIntegrator.computeIntegral(primitive1, primitive2);
        }
    }
    return integral;*/
}
