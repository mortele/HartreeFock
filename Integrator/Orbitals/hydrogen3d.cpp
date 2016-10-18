#include "hydrogen3d.h"
#include <cmath>

using std::sqrt;
using std::exp;

double Hydrogen3D::computeWavefunction(double*  coordinates,
                                       int*     quantumNumbers) {
    const double r        = coordinates[0];
    const double theta    = coordinates[1];   // [0, pi]
    const double phi      = coordinates[2];   // [0, 2pi]

    const int n = quantumNumbers[0];
    const int l = quantumNumbers[1];
    const int m = quantumNumbers[2];

    const double n1 = 2./n;
    const double n2 = factorial(n-l-1);
    const double n3 = 2*n;
    const double n4 = factorial(n+l);
    const double normalization = sqrt(n1*n1*n1 * n2 / (n3*n4*n4*n4));
    const double exponential = exp(-r/((double) n));


}

double*Hydrogen3D::getCoordinateScales() {
    double* scales = new double[3];
    scales[0] = m_rMax;
    scales[1] = m_thetaMax;
    scales[2] = m_phiMax;
    return scales;
}
