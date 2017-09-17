#include "restrictedhartreefock.h"
#include <armadillo>
#include <cassert>
#include <iostream>
#include <cmath>
#include <iomanip>

using std::cout;
using std::endl;

RestrictedHartreeFock::RestrictedHartreeFock(int nrOfParticles, int nrOfSpinOrbitals) {

    assert(nrOfParticles > 0);
    assert(nrOfParticles % 2 == 0);

    m_nrOfParticles        = nrOfParticles;
    m_nrOfSpinOrbitals     = nrOfSpinOrbitals;
    m_nrOfSpatialOrbitals  = m_nrOfSpinOrbitals/2;
    m_nrOfOccupiedOrbitals = m_nrOfParticles/2;

    //std::cout << m_nrOfSpatialOrbitals << std::endl;
    //std::cout << m_nrOfOccupiedOrbitals << std::endl;

    m_U               = arma::eye(m_nrOfSpatialOrbitals,m_nrOfSpatialOrbitals); //Intitial guess U0 = Identity matrix
    m_DensityMatrix   = arma::zeros<arma::mat>(m_nrOfSpatialOrbitals,m_nrOfSpatialOrbitals);
    m_FockMatrix      = arma::zeros<arma::mat>(m_nrOfSpatialOrbitals,m_nrOfSpatialOrbitals);
    m_oneBodyElements = arma::zeros<arma::mat>(m_nrOfSpatialOrbitals,m_nrOfSpatialOrbitals);
    m_eps             = arma::zeros<arma::vec>(m_nrOfSpatialOrbitals);
    m_eps_old         = arma::zeros<arma::vec>(m_nrOfSpatialOrbitals);

    m_integralTable = new IntegralTable();

}

void RestrictedHartreeFock::setAnalyticOneBodyElements(arma::vec oneBodyElements) {
    for(int i = 0; i < m_nrOfSpatialOrbitals; i++) {
        m_oneBodyElements(i,i) = oneBodyElements(i);
    }
}

void RestrictedHartreeFock::setAnalyticOneBodyElements(arma::mat oneBodyElements) {
    m_oneBodyElements = oneBodyElements;
}

bool RestrictedHartreeFock::setIntegralTable(std::string fileName) {
    return m_integralTable->readTableFromFile(fileName);
}

void RestrictedHartreeFock::setOverLapMatrix(arma::mat overLap) {
    m_overlap = true;
    m_S = overLap;
}

void RestrictedHartreeFock::setMaximumIterations(int maximum) {
    m_MAX_ITERS = maximum;
}

void RestrictedHartreeFock::computeSolutionBySCF() {

    computeDensityMatrix();
    computeFockMatrix();
    diagonalizeFockMatrix();
    m_eps_old = m_eps;
    m_nr_OfIters += 1;

    double energy = computeHartreeFockEnergy();

    for(int k = 1; k < m_MAX_ITERS; k++) {

        computeDensityMatrix();
        computeFockMatrix();
        diagonalizeFockMatrix(); //computes m_eps and m_u
        m_nr_OfIters += 1;

        energy = computeHartreeFockEnergy();
        //cout << std::setprecision(16) << energy << endl;

        if(arma::abs(m_eps-m_eps_old).max() < m_convergencePrecision) {
            m_eps_old = m_eps;
            m_reachedSelfConsistency = true;
            break;
        }

        if (k != m_MAX_ITERS-1) {
            m_eps_old = m_eps;
        }

    }

    m_HartreeFockEnergy = energy;
    this->printInfo();

}

void RestrictedHartreeFock::printInfo() {

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << "Number of particles: " << m_nrOfParticles << endl;
    cout << "Number of spatial orbitals: " << m_nrOfSpatialOrbitals << endl;
    cout << endl;
    cout << "  -- SCF-iteration conditions -- " << endl;
    cout << "Convergence precision : 10^" << std::log10(m_convergencePrecision) << endl;
    cout << "Maximum SCF-iterations: " << m_MAX_ITERS << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    if(m_reachedSelfConsistency) {
        cout << "Self Consistency Reached: True"  << endl;
        cout << "Number of SCF-iterations before convergence: " << m_nr_OfIters << endl;
        cout << "E_rhf: " << std::setprecision(-std::log10(m_convergencePrecision)) << m_HartreeFockEnergy << endl;
    } else {
        cout << "Self Consistency Reached: False" << endl;
        cout << "SCF did not converge to the given precision" << endl;
        cout << "|eps-stuff|: " << arma::abs(m_eps-m_eps_old).max() << endl;
        cout << "E_rhf: " << std::setprecision(-std::log10(m_convergencePrecision)) << m_HartreeFockEnergy << endl;
    }

    cout << endl;
    /*
    for(int i = 0; i < m_nrOfOccupiedOrbitals; i++) {
        cout << m_U.col(i) << endl;
    }
    */

}


void RestrictedHartreeFock::diagonalizeFockMatrix() {

    /* Note that eig_sym returns eigenvalues in ascending order
     * Should also then have eigenvectors in corresponing columns?
     */

    if(m_overlap) {
        m_U = m_S*m_U;
    }

    arma::eig_sym(m_eps, m_U, m_FockMatrix);

}

double RestrictedHartreeFock::getOneBodyMatrixElement(int p, int q) {
    return m_oneBodyElements(p,q);
}

double RestrictedHartreeFock::QRPS_AntiSym(int q, int r, int p, int s) {
    return m_integralTable->getIntegral(q,r,p,s) - 0.5*m_integralTable->getIntegral(q,r,s,p);
}

void RestrictedHartreeFock::computeDensityMatrix() {

    for(int r = 0; r < m_nrOfSpatialOrbitals; r++) {
        for(int s = 0; s < m_nrOfSpatialOrbitals; s++) {
            m_DensityMatrix(s,r) = 0.0;
            for(int j = 0; j < m_nrOfOccupiedOrbitals; j++) {
                m_DensityMatrix(s,r) += 2.0*m_U(s,j)*m_U(r,j); //If U complex we must have conj(m_U(r,j))
            }
        }
    }

}

void RestrictedHartreeFock::computeFockMatrix() {

    for(int q = 0; q < m_nrOfSpatialOrbitals; q++) {
        for(int p = 0; p < m_nrOfSpatialOrbitals; p++) {

            m_FockMatrix(q,p) = 0;

            m_FockMatrix(q,p) = getOneBodyMatrixElement(q,p);

            double tmp = 0.0;

            for(int r = 0; r < m_nrOfSpatialOrbitals; r++) {
                for(int s = 0; s < m_nrOfSpatialOrbitals; s++) {
                    tmp += m_DensityMatrix(s,r)*QRPS_AntiSym(q,r,p,s);
                }
            }
            m_FockMatrix(q,p) += tmp;
        }
    }
}

double RestrictedHartreeFock::computeHartreeFockEnergy() {

    double ERHF = 0.0;

    for(int i = 0; i < m_nrOfOccupiedOrbitals; i++) {
        ERHF += m_eps(i);
    }

    ERHF *= 2.0;

    for(int i = 0; i < m_nrOfOccupiedOrbitals; i++) {
        for(int p = 0; p < m_nrOfSpatialOrbitals; p++) {
            for(int q = 0; q < m_nrOfSpatialOrbitals; q++) {

                double tmp = 0.0;

                for(int s = 0; s < m_nrOfSpatialOrbitals; s++) {
                    for(int r = 0; r < m_nrOfSpatialOrbitals; r++) {
                        tmp += m_DensityMatrix(s,r)*QRPS_AntiSym(q,r,p,s);
                    }
                }

                ERHF -= m_U(q,i)*tmp*m_U(p,i); //Assume U real, else conjugate(m_U[q,i])
            }
        }
    }

    m_HartreeFockEnergy = ERHF;
    return ERHF;
}
