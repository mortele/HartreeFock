#include "Integrators/electronnucleusintegrator.h"
#include <cmath>

using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;

ElectronNucleusIntegrator::ElectronNucleusIntegrator() :
        m_nucleusPosition(zeros<vec>(3)),
        m_hermiteGaussian(HermiteGaussian()),
        m_hermiteGaussianIntegral(HermiteGaussianIntegral()) {

}

void ElectronNucleusIntegrator::setNucleusPosition(arma::vec nucleusPosition) {
    m_nucleusPosition = nucleusPosition;
}

double ElectronNucleusIntegrator::computeIntegral(GaussianPrimitive* primitive1,
                                                  GaussianPrimitive* primitive2) {

    m_hermiteGaussianIntegral.setupCoefficients(primitive1,
                                                primitive2,
                                                m_nucleusPosition);
    m_hermiteGaussian.setupCoefficients(primitive1,
                                        primitive2);
    const double p      = primitive1->exponent()  + primitive2->exponent();

    const int    x1     = primitive1->xExponent();
    const int    y1     = primitive1->yExponent();
    const int    z1     = primitive1->zExponent();

    const int    x2     = primitive2->xExponent();
    const int    y2     = primitive2->yExponent();
    const int    z2     = primitive2->zExponent();

    const int    tLimit = x1 + x2 + 1;
    const int    uLimit = y1 + y2 + 1;
    const int    vLimit = z1 + z2 + 1;

    double integral = 0;
    for (int t = 0; t < tLimit; t++) {
        for (int u = 0; u < uLimit; u++) {
            for (int v = 0; v < vLimit; v++) {
                double Eproduct = 1;
                Eproduct *= m_hermiteGaussian.getCoefficientDimension(x1,x2,t,0);
                Eproduct *= m_hermiteGaussian.getCoefficientDimension(y1,y2,u,1);
                Eproduct *= m_hermiteGaussian.getCoefficientDimension(z1,z2,v,2);
                integral += Eproduct * m_hermiteGaussianIntegral.getCoefficient(0,t,u,v);
            }
        }
    }
    return integral * 2*M_PI / p;
}

double ElectronNucleusIntegrator::computeIntegral(GaussianPrimitive* primitive1,
                                                  GaussianPrimitive* primitive2,
                                                  arma::vec          nucleusPosition) {
    setNucleusPosition(nucleusPosition);
    return computeIntegral(primitive1, primitive2);
}
