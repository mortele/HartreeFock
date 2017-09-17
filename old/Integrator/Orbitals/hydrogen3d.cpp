#include "hydrogen3d.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::sqrt;
using std::exp;
using std::cos;
using std::sin;
using std::pow;

double Hydrogen3D::solidHarmonics(int    l,
                                  int    m,
                                  double r,
                                  double theta,
                                  double phi) {

    //cout << "l,m,r,theta,phi=" << l << "," << m << "," << r << "," << theta << "," << phi << endl;
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

double Hydrogen3D::normalization(int n, int l, int m) {
    if (n==1 && l==0 && m==0) {
        return 0.5641718413;
    } else if (n==2 && l==0 && m==0) {
        return 0.09972601571;
    } else if (n==2 && l==1 && m==-1) {
        return 0.09976101671;
    } else if (n==2 && l==1 && m==0) {
        return 0.09972977957;
    } else if (n==2 && l==1 && m==1) {
        return 0.09972259164;
    } else if (n==3 && l==0 && m==0) {
        return 0.03619545492;
    } else if (n==3 && l==1 && m==-1) {
        return 0.009675197503;
    } else if (n==3 && l==1 && m==0) {
        return 0.009670842141;
    } else if (n==3 && l==1 && m==1) {
        return 0.009672302328;
    } else if (n==3 && l==2 && m==-2) {
        return 0.005687606665;
    } else if (n==3 && l==2 && m==-1) {
        return 0.00568748345;
    } else if (n==3 && l==2 && m==0) {
        return 0.005687948022;
    } else if (n==3 && l==2 && m==1) {
        return 0.005687113449;
    } else if (n==3 && l==2 && m==2) {
        return 0.005687235524;
    } else {
        return 1;
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

    const double norm           = normalization(n,l,m);
    const double exponential    = exp(-r/((double) n));
    const double solid          = solidHarmonics(l,m,r,theta,phi);
    const double laguerre       = Orbital::associatedLaguerrePolynomial(2*r/((double) n),
                                                                        n-l-1,
                                                                        2*l-1);
    return norm * exponential * solid * laguerre;
}

double Hydrogen3D::integrandOne(double* allCoordinates, int* allQuantumNumbers) {
    double r                    = allCoordinates[0];
    double theta                = allCoordinates[1];
    double integrationMeasure   = r*r*sin(theta);

    double waveFunction1 = computeWavefunction(allCoordinates, allQuantumNumbers);
    double waveFunction2 = computeWavefunction(allCoordinates, allQuantumNumbers+3);
    double integrand     = integrationMeasure * waveFunction1 * waveFunction2;
    return integrand;
}

double Hydrogen3D::integrandTwo(double* allCoordinates, int* allQuantumNumbers) {
    const double r1                   = allCoordinates[0];
    const double theta1               = allCoordinates[1];
    const double r2                   = allCoordinates[3];
    const double theta2               = allCoordinates[4];
    const double integrationMeasure   = r1*r1*sin(theta1) * r2*r2*sin(theta2);
    const double r12                  = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos(theta2-theta1));
    const double oneOverR12           = r12 < 1e-12 ? 0 : 1./r12;

    const double waveFunction1 = computeWavefunction(allCoordinates,   allQuantumNumbers);
    const double waveFunction2 = computeWavefunction(allCoordinates+3, allQuantumNumbers+3);
    const double waveFunction3 = computeWavefunction(allCoordinates,   allQuantumNumbers+6);
    const double waveFunction4 = computeWavefunction(allCoordinates+3, allQuantumNumbers+9);
    const double integrand     = integrationMeasure * waveFunction1 * waveFunction2 *
                                 oneOverR12         * waveFunction3 * waveFunction4;
    return integrand;
}

double*Hydrogen3D::getCoordinateScales() {
    double* scales = new double[3];
    scales[0] = m_rMax;
    scales[1] = m_thetaMax;
    scales[2] = m_phiMax;
    return scales;
}

void Hydrogen3D::updateCoordinateScales(int* allQuantumNumbers,
                                        int  numberOfQuantumNumbers) {
    int nMax = 0;
    int lMax = 0;
    for (int i=0; i<numberOfQuantumNumbers; i+=3) {
        if (allQuantumNumbers[i] >= nMax) {
            nMax = allQuantumNumbers[i];
            if (allQuantumNumbers[i+1] > lMax) {
                lMax = allQuantumNumbers[i+1];
            }
        }
    }
    if (nMax==1) {
        m_rMax = 12;
    } else if (nMax==2) {
        m_rMax = 16;
    } else if (nMax==3) {
        m_rMax = 35;
    } else {
        m_rMax = 12;
    }
}







