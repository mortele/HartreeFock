#include "numericalintegrator.h"
#include "system.h"
#include "Integrators/grid.h"
#include <cmath>
#include "Atoms/atom.h"
#include "Orbitals/contractedgaussian.h"
#include "Orbitals/gaussianprimitive.h"
#include "Solvers/hartreefock.h"


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

double NumericalIntegrator::testIntegral(const arma::mat& densityMatrix) {
    m_grid->createSimpleOneAtomGrid(200,30);

    const vec& weights = m_grid->getWeights();
    const mat& points = m_grid->getPoints();

    std::vector<ContractedGaussian*> basis = m_system->getBasis();
    int basisSize = basis.size();

    double integral = 0;
    for (int i = 0; i < weights.n_elem; i++) {
        const double x = points(i,0);
        const double y = points(i,1);
        const double z = points(i,2);
        const double w = weights(i);
        double tmp = 0;
        for (int p = 0; p < basisSize; p++) {
            ContractedGaussian* pPhi = basis.at(p);
            for (int q = 0; q < basisSize; q++) {
                ContractedGaussian* qPhi = basis.at(q);
                tmp += densityMatrix(p,q) * pPhi->evaluate(x,y,z) * qPhi->evaluate(x,y,z);
            }
        }
        integral += w * tmp;
    }
    return integral;
}

int NumericalIntegrator::generateBeckeGrid() {
    //double radialPrecision = 1e-12;
    //int maximumRadialPoints = 302;
    //int minimumRadialPoints = 86;
    double radialPrecision = 1e-20;
    int maximumRadialPoints = 2000;
    int minimumRadialPoints = 1000;

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
        if (numberOfAtoms==1) {
            basisCenters[i] = 1;
        } else {
            for (int j = 0; j < numberOfAtoms; j++) {
                const vec& position = m_system->getAtoms().at(j)->getPosition();
                if ((fabs(position(0) - center(0)) < 1e-5) &&
                    (fabs(position(1) - center(1)) < 1e-5) &&
                    (fabs(position(2) - center(2)) < 1e-5)) {
                    basisCenters[i] = j+1;
                }
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
    int     numberOfOuterAtoms = 0;
    int*    outerAtomSize = NULL;
    double* outerAtomCoordinates = NULL;

    m_context = numgrid_new_context();
    int error = numgrid_generate_grid(m_context,                   //  context,
                                      radialPrecision,             //  radial_precision,
                                      minimumRadialPoints,         //  min_num_angular_points,
                                      maximumRadialPoints,         //  max_num_angular_points,
                                      numberOfAtoms,               //  num_centers,
                                      atomCoordinates,             //  center_coordinates,
                                      atomSize,                    //  center_elements,
                                      numberOfOuterAtoms,          //  num_outer_centers,
                                      outerAtomCoordinates,        //  outer_center_coordinates,
                                      outerAtomSize,               //  outer_center_elements,
                                      numberOfBasisFunctions,      //  num_shells,
                                      basisCenters,                //  shell_centers,
                                      basisAngularMomentum,        //  shell_l_quantum_numbers,
                                      basisNumberOfPrimitives,     //  shell_num_primitives,
                                      primitiveExponents);         //  primitive_exponents);
    return error;
}

double NumericalIntegrator::integrateDensity(const mat& densityMatrix) {
    generateBeckeGrid();
    int             numberOfGridPoints  = numgrid_get_num_points(m_context);
    const double*   grid                = numgrid_get_grid      (m_context);

    std::vector<ContractedGaussian*> basis = m_system->getBasis();
    int basisSize = basis.size();

    cout << "GRID SIZE: " << numberOfGridPoints << endl;
    double integral = 0;
    for (int i = 0; i < 4*numberOfGridPoints; i+=4) {
        const double x = grid[i+0];
        const double y = grid[i+1];
        const double z = grid[i+2];
        const double w = grid[i+3];

        /*if (i%40==0) {
            printf("%15.10f %15.10f %15.10f %15.8g \n",x,y,z,w);
        }*/

        double tmp = 0;
        for (int p = 0; p < basisSize; p++) {
            ContractedGaussian* pPhi = basis.at(p);
            for (int q = 0; q < basisSize; q++) {
                ContractedGaussian* qPhi = basis.at(q);
                tmp += densityMatrix(p,q) * pPhi->evaluate(x,y,z) * qPhi->evaluate(x,y,z);
            }
        }
        integral += w * tmp;
    }
    return integral;
}









