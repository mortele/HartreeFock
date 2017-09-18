#include "unrestrictedhartreefock.h"
#include "Orbitals/contractedgaussian.h"
#include "Orbitals/gaussianprimitive.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <boost/format.hpp>
#include "system.h"

using std::cout;
using std::endl;
using std::printf;
using std::sqrt;
using arma::eye;
using arma::zeros;
using arma::mat;
using arma::vec;
using arma::field;
using arma::randu;


UnrestrictedHartreeFock::UnrestrictedHartreeFock(System* system) :
        HartreeFock(system) {

    m_numberOfSpinUpElectrons   = m_system->getNumberOfSpinUpElectrons();
    m_numberOfSpinDownElectrons = m_system->getNumberOfSpinDownElectrons();

    m_epsilonUp                   = zeros(m_numberOfSpinUpElectrons);
    m_epsilonOldUp                = zeros(m_numberOfSpinUpElectrons);
    m_fockMatrixUp                = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTildeUp           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixUp         = zeros(m_numberOfBasisFunctions, m_numberOfSpinUpElectrons);
    m_coefficientMatrixTildeUp    = zeros(m_numberOfBasisFunctions, m_numberOfSpinUpElectrons);
    //m_densityMatrixUp             = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_densityMatrixUp             = m_coefficientMatrixUp   * m_coefficientMatrixUp.t();
    m_densityMatrixUp(0,1) = 0.1;

    m_epsilonDown                 = zeros(m_numberOfSpinDownElectrons);
    m_epsilonOldDown              = zeros(m_numberOfSpinDownElectrons);
    m_fockMatrixDown              = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_fockMatrixTildeDown         = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_coefficientMatrixDown       = zeros(m_numberOfBasisFunctions, m_numberOfSpinDownElectrons);
    m_coefficientMatrixTildeDown  = zeros(m_numberOfBasisFunctions, m_numberOfSpinDownElectrons);
    //m_densityMatrixDown           = zeros(m_numberOfBasisFunctions, m_numberOfBasisFunctions);
    m_densityMatrixDown           = m_coefficientMatrixDown * m_coefficientMatrixDown.t();
}

void UnrestrictedHartreeFock::setup() {
    setupOverlapMatrix();
    diagonalizeOverlapMatrix();
    setupOneBodyElements();
    setupTwoBodyElements();
    computeDensityMatrices();
    m_nucleusNucleusInteractionEnergy = m_system->nucleusNucleusInteractionEnergy();
}

void UnrestrictedHartreeFock::computeFockMatrices() {
    for(int p = 0; p < m_numberOfBasisFunctions; p++)
    for(int q = 0; q < m_numberOfBasisFunctions; q++) {
    //for(int q = p; q < m_numberOfBasisFunctions; q++) {

        m_fockMatrixUp(p,q)   = m_oneBodyMatrixElements(p,q);
        m_fockMatrixDown(p,q) = m_oneBodyMatrixElements(p,q);

        for(int r = 0; r < m_numberOfBasisFunctions; r++)
        for(int s = 0; s < m_numberOfBasisFunctions; s++) {
            const double prqs = m_twoBodyMatrixElements(p,r)(q,s);
            const double prsq = m_twoBodyMatrixElements(p,r)(s,q);
            m_fockMatrixUp(p,q)   += (prqs-prsq) * m_densityMatrixUp(s,r)
                                   +  prqs       * m_densityMatrixDown(s,r);
            m_fockMatrixDown(p,q) += (prqs-prsq) * m_densityMatrixDown(s,r)
                                   +  prqs       * m_densityMatrixUp(s,r);
        }
    }
}

void UnrestrictedHartreeFock::computeDensityMatrices() {
    if (m_smoothing) {
        double a = m_smoothingFactor;
        mat densityMatrixUpTmp   = m_coefficientMatrixUp   * m_coefficientMatrixUp.t();;
        mat densityMatrixDownTmp = m_coefficientMatrixDown * m_coefficientMatrixDown.t();;
        m_densityMatrixUp   = a*m_densityMatrixUp    + (1-a)*densityMatrixUpTmp;
        m_densityMatrixDown = a*m_densityMatrixDown  + (1-a)*densityMatrixDownTmp;
    } else {
        m_densityMatrixUp   = m_coefficientMatrixUp   * m_coefficientMatrixUp.t();
        m_densityMatrixDown = m_coefficientMatrixDown * m_coefficientMatrixDown.t();
    }
}

void UnrestrictedHartreeFock::diagonalizeFockMatrices() {
    const mat& A = m_transformationMatrix;
    m_fockMatrixTildeUp   = A.t() * m_fockMatrixUp   * A;
    m_fockMatrixTildeDown = A.t() * m_fockMatrixDown * A;
    arma::eig_sym(m_epsilonUp,   m_coefficientMatrixTildeUp,   m_fockMatrixTildeUp);
    arma::eig_sym(m_epsilonDown, m_coefficientMatrixTildeDown, m_fockMatrixTildeDown);
    m_coefficientMatrixUp   = A * m_coefficientMatrixTildeUp.submat(0,0,m_numberOfBasisFunctions-1, m_numberOfSpinUpElectrons-1);
    m_coefficientMatrixDown = A * m_coefficientMatrixTildeDown.submat(0,0,m_numberOfBasisFunctions-1, m_numberOfSpinDownElectrons-1);
    normalizeCoefficientMatrices();
}

void UnrestrictedHartreeFock::normalizeCoefficientMatrices() {
    for (int k = 0; k < m_numberOfSpinUpElectrons; k++) {
        double normalizationFactor = 0;
        for (int p = 0; p < m_numberOfBasisFunctions; p++) {
            for (int q = 0; q < m_numberOfBasisFunctions; q++) {
                normalizationFactor += m_coefficientMatrixUp(p,k) *
                                       m_coefficientMatrixUp(q,k) *
                                       m_overlapMatrix(p,q);
            }
        }
        normalizationFactor = sqrt(normalizationFactor);
        normalizationFactor = (normalizationFactor == 0.0 ? 1. : normalizationFactor);
        m_coefficientMatrixUp.col(k) = m_coefficientMatrixUp.col(k) /
                                       normalizationFactor;
    }
    for (int k = 0; k < m_numberOfSpinDownElectrons; k++) {
        double normalizationFactor = 0;
        for (int p = 0; p < m_numberOfBasisFunctions; p++) {
            for (int q = 0; q < m_numberOfBasisFunctions; q++) {
                normalizationFactor += m_coefficientMatrixDown(p,k) *
                                       m_coefficientMatrixDown(q,k) *
                                       m_overlapMatrix(p,q);
            }
        }
        normalizationFactor = sqrt(normalizationFactor);
        normalizationFactor = (normalizationFactor == 0.0 ? 1. : normalizationFactor);
        m_coefficientMatrixDown.col(k) = m_coefficientMatrixDown.col(k) /
                                       normalizationFactor;
    }
}

void UnrestrictedHartreeFock::selfConsistentFieldIteration() {
    computeFockMatrices();
    diagonalizeFockMatrices();
    computeDensityMatrices();
}

void UnrestrictedHartreeFock::computeHartreeFockEnergy() {
    m_hartreeFockEnergy = 0;

    for (int p = 0; p < m_numberOfBasisFunctions; p++)
    for (int q = 0; q < m_numberOfBasisFunctions; q++) {
        m_hartreeFockEnergy += m_oneBodyMatrixElements(p,q)
                               * (m_densityMatrixUp(p,q) + m_densityMatrixDown(p,q));
        m_hartreeFockEnergy += m_fockMatrixUp(p,q)   * m_densityMatrixUp(p,q);
        m_hartreeFockEnergy += m_fockMatrixDown(p,q) * m_densityMatrixDown(p,q);
    }
    m_hartreeFockEnergy *= 0.5;
    m_hartreeFockEnergy += m_nucleusNucleusInteractionEnergy;
}

void UnrestrictedHartreeFock::storeEnergy() {
    m_epsilonOldUp   = m_epsilonUp;
    m_epsilonOldDown = m_epsilonDown;
}

double UnrestrictedHartreeFock::convergenceTest() {
    double sumUp    = arma::sum(arma::abs(m_epsilonUp   - m_epsilonOldUp))   / m_epsilonUp.n_elem;
    double sumDown  = arma::sum(arma::abs(m_epsilonDown - m_epsilonOldDown)) / m_epsilonDown.n_elem;
    return 0.5 * (sumUp + sumDown);
}

std::string UnrestrictedHartreeFock::getDateTime() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d-%H.%M.%S");
    std::string dateTime = ss.str();
    return dateTime;
}

std::string UnrestrictedHartreeFock::dumpBasisToFile(std::string fileNameIn) {
    std::string fileName;
    if (fileNameIn=="") {
        fileName = "../HartreeFockBases/basis-" + getDateTime();
    } else {
        fileName = "../HartreeFockBases/" + fileNameIn;
    }

    std::vector<ContractedGaussian*> basis = m_system->getBasis();

    std::ofstream outFile;
    outFile.open(fileName, std::ios::out);

    int basisSize = m_numberOfBasisFunctions;
    int electrons = m_numberOfElectrons;
    int spinUp    = m_numberOfSpinUpElectrons;
    int spinDown  = m_numberOfSpinDownElectrons;
    int atoms     = m_system->getAtoms().size();
    outFile << boost::format("%d %d %d %d %d\n") % electrons % spinUp % spinDown % basisSize % atoms;

    for (Atom* atom : m_system->getAtoms()) {
        vec position = atom->getPosition();
        int Z        = atom->getCharge();
        outFile << boost::format("%d %.15f %.15f %.15f\n") % Z % position(0) % position(1) % position(2);
    }

    for (ContractedGaussian* contracted : basis) {
        int primitives = contracted->getPrimitives().size();
        vec nucleus = contracted->getNucleusPosition();

        outFile << boost::format("%d\n") % primitives;
        outFile << boost::format("%.15f %.15f %.15f\n") % nucleus(0) % nucleus(1) % nucleus(2);

        for (GaussianPrimitive* primitive : contracted->getPrimitives()) {
            int x = primitive->xExponent();
            int y = primitive->yExponent();
            int z = primitive->zExponent();
            double exponent = primitive->exponent();
            double constant = primitive->getCoefficient();
            outFile << boost::format("%d %d %d %.15f %.15f\n") % x % y % z % exponent % constant;
        }
    }
    for (int electron = 0; electron < m_numberOfSpinUpElectrons; electron++) {
        for (int basisFunction = 0; basisFunction < m_numberOfBasisFunctions; basisFunction++) {
            outFile << std::setprecision(5) << m_coefficientMatrixUp(basisFunction, electron) << " ";
        }
        outFile << "\n";
    }
    for (int electron = 0; electron < m_numberOfSpinDownElectrons; electron++) {
        for (int basisFunction = 0; basisFunction < m_numberOfBasisFunctions; basisFunction++) {
            outFile << std::setprecision(5) << m_coefficientMatrixDown(basisFunction, electron) << " ";
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Dumped HF basis to file: " << fileName << std::endl;
    return fileName;
}






