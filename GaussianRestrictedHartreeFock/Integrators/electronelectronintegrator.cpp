#include "electronelectronintegrator.h"
#include <cmath>

using std::sqrt;
using std::cout;
using std::endl;
using arma::vec;

ElectronElectronIntegrator::ElectronElectronIntegrator() :
    m_2sqrtPiToThe5(2*sqrt(M_PI*M_PI*M_PI*M_PI*M_PI)),
    m_hermiteGaussian12(HermiteGaussian()),
    m_hermiteGaussian34(HermiteGaussian()),
    m_hermiteGaussianIntegral(HermiteGaussianIntegral()) {

}

void ElectronElectronIntegrator::setupHermiteGaussianIntegral(GaussianPrimitive& primitive1,
                                                              GaussianPrimitive& primitive2,
                                                              GaussianPrimitive& primitive3,
                                                              GaussianPrimitive& primitive4) {

    int t = primitive1.xExponent() + primitive2.xExponent() +
            primitive3.xExponent() + primitive4.xExponent();
    int u = primitive1.yExponent() + primitive2.yExponent() +
            primitive3.yExponent() + primitive4.yExponent();
    int v = primitive1.zExponent() + primitive2.zExponent() +
            primitive3.zExponent() + primitive4.zExponent();
    double alpha = primitive1.exponent();
    double beta  = primitive2.exponent();
    double gamma = primitive3.exponent();
    double delta = primitive4.exponent();

    double p1 = (alpha + beta);
    double p2 = (gamma + delta);
    double p  = p1 * p2 / (p1 + p2);

    vec nucleus1    = primitive1.nucleusPosition();
    vec nucleus2    = primitive2.nucleusPosition();
    vec nucleus3    = primitive3.nucleusPosition();
    vec nucleus4    = primitive4.nucleusPosition();
    vec P1          = (alpha * nucleus1 + beta  * nucleus2) / (alpha + beta);
    vec P2          = (gamma * nucleus3 + delta * nucleus4) / (gamma + delta);
    vec PC          = P1 - P2;

    m_hermiteGaussianIntegral.setupCoefficients(t, u, v, p, PC);
}

double ElectronElectronIntegrator::computeIntegral(GaussianPrimitive& primitive1,
                                                   GaussianPrimitive& primitive2,
                                                   GaussianPrimitive& primitive3,
                                                   GaussianPrimitive& primitive4) {
    m_hermiteGaussian12.setupCoefficients(primitive1, primitive2);
    m_hermiteGaussian34.setupCoefficients(primitive3, primitive4);
    setupHermiteGaussianIntegral(primitive1, primitive2, primitive3, primitive4);

    const int x1 = primitive1.xExponent();
    const int y1 = primitive1.yExponent();
    const int z1 = primitive1.zExponent();

    const int x2 = primitive2.xExponent();
    const int y2 = primitive2.yExponent();
    const int z2 = primitive2.zExponent();

    const int x3 = primitive3.xExponent();
    const int y3 = primitive3.yExponent();
    const int z3 = primitive3.zExponent();

    const int x4 = primitive4.xExponent();
    const int y4 = primitive4.yExponent();
    const int z4 = primitive4.zExponent();

    int tuvLimits[] {x1 + x2 + 1,   y1 + y2 + 1,   z1 + z2 + 1,
                     x3 + x4 + 1,   y3 + y4 + 1,   z3 + z4 + 1};

    double integral = 0;

    for (int t = 0; t < tuvLimits[0]; t++) {
        for (int u = 0; u < tuvLimits[1]; u++) {
            for (int v = 0; v < tuvLimits[2]; v++) {
                for (int t_ = 0; t_ < tuvLimits[3]; t_++) {
                    for (int u_ = 0; u_ < tuvLimits[4]; u_++) {
                        for (int v_ = 0; v_ < tuvLimits[5]; v_++) {

                            double Eproduct = 1;
                            Eproduct *= m_hermiteGaussian12.getCoefficientDimension(x1,x2,t,0);
                            Eproduct *= m_hermiteGaussian12.getCoefficientDimension(y1,y2,u,1);
                            Eproduct *= m_hermiteGaussian12.getCoefficientDimension(z1,z2,v,2);

                            Eproduct *= m_hermiteGaussian34.getCoefficientDimension(x3,x4,t_,0);
                            Eproduct *= m_hermiteGaussian34.getCoefficientDimension(y3,y4,u_,1);
                            Eproduct *= m_hermiteGaussian34.getCoefficientDimension(z3,z4,v_,2);


                            double R = m_hermiteGaussianIntegral.getCoefficient(0, t+t_, u+u_, v+v_);
                            double sign = ((t_ + u_ + v_) % 2) == 0 ? 1 : -1;

                            integral += Eproduct * R * sign;
                        }
                    }
                }
            }
        }
    }

    double p1       = primitive1.exponent() + primitive2.exponent();
    double p2       = primitive3.exponent() + primitive4.exponent();
    return m_2sqrtPiToThe5 * integral / (p1*p2*sqrt(p1+p2));
}






















