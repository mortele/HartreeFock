#include "numericalintegrator.h"
#include "system.h"
#include "Integrators/grid.h"
#include <cmath>
#include <numgrid.h>
#include "Atoms/atom.h"
#include "Orbitals/contractedgaussian.h"
#include "Orbitals/gaussianprimitive.h"


using arma::vec;
using arma::mat;
using std::sqrt;
using std::exp;
using std::cout;
using std::endl;


NumericalIntegrator::NumericalIntegrator(System* system) {
    m_system = system;
    m_grid   = new Grid(system);
}

double NumericalIntegrator::testIntegral() {
    m_grid->createSimpleOneAtomGrid(100,150,5);

    const vec& w = m_grid->getWeights();
    const mat& p = m_grid->getPoints();

    double integral = 0;
    for (int i = 0; i < w.n_elem; i++) {
        const double x = p(i,0);
        const double y = p(i,1);
        const double z = p(i,2);
        const double r = sqrt(x*x + y*y + z*z);
        integral += exp(-r)*w(i);
    }
    return integral;
}

double NumericalIntegrator::testBecke() {
    double radialPrecision = 1e-12;
    int maximumRadialPoints = 302;
    int minimumRadialPoints = 86;

    int numberOfAtoms = m_system->getAtoms().size();
    double atomCoordinates[numberOfAtoms*3];
    int atomSize[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++) {
        const vec& position = m_system->getAtoms().at(i)->getPosition();
        atomCoordinates[i*3+0] = position(0);
        atomCoordinates[i*3+1] = position(1);
        atomCoordinates[i*3+2] = position(2);
        atomSize[i] = m_system->getAtoms().at(i)->getNumberOfElectrons();
    }

    int numberOfBasisFunctions = m_system->getBasis().size();
    int basisCenters[numberOfBasisFunctions];
    int basisAngularMomentum[numberOfBasisFunctions];
    int basisNumberOfPrimitives[numberOfBasisFunctions];
    int totalNumberOfPrimitives = 0;

    for (int i = 0; i < numberOfBasisFunctions; i++) {
        GaussianPrimitive* primitive = m_system->getBasis().at(i)->getPrimitives().at(0);
        basisAngularMomentum[i] = primitive->xExponent() +
                                  primitive->yExponent() +
                                  primitive->zExponent();

        const vec& center = primitive->nucleusPosition();
        for (int j = 0; j < numberOfAtoms; j++) {
            const vec& position = m_system->getAtoms().at(j)->getPosition();
            if ((fabs(position(0) - center(0)) < 1e-5) &&
                (fabs(position(1) - center(1)) < 1e-5) &&
                (fabs(position(2) - center(2)) < 1e-5)) {
                basisCenters[i] = j;
            }
        }
        basisNumberOfPrimitives[i] = m_system->getBasis().at(i)->getPrimitives().size();
        totalNumberOfPrimitives += basisNumberOfPrimitives[i];
    }

    double primitiveExponents[totalNumberOfPrimitives];
    int index = 0;
    for (int i = 0; i < numberOfBasisFunctions; i++) {
        ContractedGaussian* contracted = m_system->getBasis().at(i);
        for (int j = 0; j < contracted->getPrimitives().size(); j++) {
            GaussianPrimitive* primitive = contracted->getPrimitives().at(j);
            primitiveExponents[index] = primitive->exponent();
            index++;
        }
    }

    context_t* m_context = numgrid_new_context();
    int error = numgrid_generate_grid(m_context,
                                      radialPrecision,
                                      minimumRadialPoints,
                                      maximumRadialPoints,
                                      numberOfAtoms,
                                      atomCoordinates,
                                      basisCenters,
                                      0,
                                      NULL,
                                      NULL,
                                      numberOfBasisFunctions,
                                      basisCenters,
                                      basisAngularMomentum,
                                      basisNumberOfPrimitives,
                                      primitiveExponents);









}
