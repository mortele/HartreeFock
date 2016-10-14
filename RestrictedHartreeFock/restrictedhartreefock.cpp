#include "restrictedhartreefock.h"
#include <armadillo>
#include <cassert>
#include <iostream>

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

void RestrictedHartreeFock::computeSolutionBySCF() {
    //double convergencePrecision = 1e-8;
    computeDensityMatrix();
    computeFockMatrix();
    diagonalizeFockMatrix();
    double energy = computeHartreeFockEnergy();
    cout << m_DensityMatrix << std::endl;
    cout << energy << endl;
}


void RestrictedHartreeFock::diagonalizeFockMatrix() {

    /* Note that eig_sym returns eigenvalues in ascending order
     * Should also then have eigenvectors in corresponing columns?
     */
    arma::eig_sym(m_eps, m_U, m_FockMatrix);
}

double RestrictedHartreeFock::getOneBodyMatrixElement(int p, int q) {
    return m_oneBodyElements(p,q);
}

double RestrictedHartreeFock::getTwoBodyMatrixElement(int p, int q, int r, int s) {
    //<pq|w|rs>
    return m_integralTable->getIntegral(p,q,r,s);
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

            if(q == p) {
                m_FockMatrix(q,p) = getOneBodyMatrixElement(p,p);
            }

            for(int r = 0; r < m_nrOfSpatialOrbitals; r++) {
                for(int s = 0; s < m_nrOfSpatialOrbitals; s++) {
                    m_FockMatrix(q,p) += m_DensityMatrix(s,r)*getTwoBodyMatrixElement(q,r,p,s);
                }
            }
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
                        tmp += m_DensityMatrix(s,r)*getTwoBodyMatrixElement(q,r,p,s);
                    }
                }

                ERHF -= m_U(q,i)*tmp*m_U(p,i); //Assume U real, else conjugate(m_U[q,i])
            }
        }
    }
    m_HartreeFockEnergy = ERHF;
    return ERHF;
}
