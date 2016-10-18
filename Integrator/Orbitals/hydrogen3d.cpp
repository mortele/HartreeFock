#include "hydrogen3d.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::sqrt;
using std::exp;
using std::cos;
using std::sin;

double Hydrogen3D::solidHarmonics(int    l,
                                  int    m,
                                  double r,
                                  double theta,
                                  double phi) {
    // http://www.phy.ohiou.edu/~elster/phys5071/extras/MHJ_Ch11.pdf  pp.280
    if (l==0) {
        return 1;
    } else if(l==1) {
        if (m==-1) {
            const double z = r*cos(theta);
            return z;
        } else if (m==0) {
            const double y = r*sin(theta)*sin(phi);
            return y;
        } else if (m==1) {
            const double x = r*sin(theta)*cos(phi);
            return x;
        }
    } else if (l==2) {
        const double sqrt3 = 1.732050807568877293527446341505872366;
        if (m==-2) {
            const double x = r*sin(theta)*cos(phi);
            const double y = r*sin(theta)*sin(phi);
            return sqrt3*x*y;
        } else if (m==-1) {
            const double y = r*sin(theta)*sin(phi);
            const double z = r*cos(theta);
            return sqrt3*y*z;
        } else if (m==0) {
            const double z = r*cos(theta);
            return 0.5*(3*z*z-r*r);
        } else if (m==1) {
            const double x = r*sin(theta)*cos(phi);
            const double z = r*cos(theta);
            return sqrt3*x*z;
        } else if (m==2) {
            const double x = r*sin(theta)*cos(phi);
            const double y = r*sin(theta)*sin(phi);
            return 0.5*sqrt3*(x*x-y*y);
        }
    } else {
        cout << "Unkown solid harmonics for (l,m)=" << l << "," << m << ")." << endl;
        return 0;
    }
}

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
    const double normalization  = sqrt(n1*n1*n1 * n2 / (n3*n4*n4*n4));
    const double exponential    = exp(-r/((double) n));
    const double solid          = solidHarmonics(l,m,r,theta,phi);
    const double laguerre       = Orbital::associatedLaguerrePolynomial(2*r/((double) n),
                                                                        n-l-1,
                                                                        2*l-1);
    return normalization * exponential * solid * laguerre;
}

double Hydrogen3D::integrandOne(double* allCoordinates, int* allQuantumNumbers) {
    double r                    = allCoordinates[0];
    double theta                = allCoordinates[1];
    double integrationMeasure   = r*r*sin(theta);

    double waveFunction1 = computeWavefunction(allCoordinates, allQuantumNumbers);
    double waveFunction2 = computeWavefunction(allCoordinates, allQuantumNumbers+2);
    double integrand     = integrationMeasure * waveFunction1 * waveFunction2;
    return integrand;
}

double*Hydrogen3D::getCoordinateScales() {
    double* scales = new double[3];
    scales[0] = m_rMax;
    scales[1] = m_thetaMax;
    scales[2] = m_phiMax;
    return scales;
}
