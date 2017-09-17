#include "grid.h"
#include <cmath>
#include "system.h"
#include "Orbitals/gaussianprimitive.h"
#include "Orbitals/contractedgaussian.h"

using arma::vec;
using arma::mat;
using arma::linspace;
using arma::zeros;
using std::cout;
using std::endl;
using std::sqrt;
using std::cos;
using std::sin;
using std::acos;
using std::log;


Grid::Grid(System* system) {
    m_system = system;
}

void Grid::createSimpleOneAtomGrid(int radialPoints, int angularPoints) {

    double lowestPrimitiveExponent = 10000;
    for (ContractedGaussian* contracted : m_system->getBasis()) {
        for (GaussianPrimitive* primitive : contracted->getPrimitives()) {
            const double exponent = primitive->exponent();
            if (exponent < lowestPrimitiveExponent) {
                lowestPrimitiveExponent = exponent;
            }
        }
    }
    const double cutoffValue    = 1e-5;
    const double maximumRadius  = sqrt(log(cutoffValue)/(-lowestPrimitiveExponent));
    const double pi             = acos(-1.0);

    mat m_points  = zeros<mat>(radialPoints*angularPoints*angularPoints, 3);
    vec m_weights = zeros<vec>(radialPoints*angularPoints*angularPoints);
    vec radius  = linspace(cutoffValue, maximumRadius,  radialPoints);
    vec theta   = linspace(0,           pi,             angularPoints);
    vec phi     = linspace(0,           2*pi,           angularPoints);

    int index = 0;
    for (int i = 0; i < radialPoints; i++) {
        for (int j = 0; j  < angularPoints; j ++) {
            for (int k = 0; k < angularPoints; k++) {
                const double sinTheta       = sin(theta[j]);
                const double cosTheta       = cos(theta[j]);
                const double sinPhi         = sin(phi[k]);
                const double cosPhi         = cos(phi[k]);
                const double r              = radius[i];
                const double volumeElement  = sin(theta[j])*r*r;
                m_weights(index)      = volumeElement;
                m_points (index,0)    = r * sinTheta * cosPhi;
                m_points (index,1)    = r * sinTheta * sinPhi;
                m_points (index,2)    = r * cosTheta;
                index++;
            }
        }
    }
    vec nucleus = m_system->getAtoms().at(0)->getPosition();
    if ( (fabs(nucleus(0)) > 1e-10) ||
         (fabs(nucleus(1)) > 1e-10) ||
         (fabs(nucleus(2)) > 1e-10) ) {
        for (int i = 0; i < radialPoints*angularPoints*angularPoints; i++) {
            m_points(i,0) -= nucleus(0);
            m_points(i,1) -= nucleus(1);
            m_points(i,2) -= nucleus(2);
        }
    }
}
