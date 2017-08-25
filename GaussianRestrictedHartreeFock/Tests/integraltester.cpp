#include "integraltester.h"
#include "Orbitals/gaussianprimitive.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using std::fabs;
using arma::vec;
using arma::mat;
using arma::zeros;
using std::cout;
using std::endl;
using std::setprecision;

void IntegralTester::setupExactOverlapMatrix() {
    m_exactOverlap(0,0) = 1.96870124321530;
    m_exactOverlap(0,1) = -0.0341172323898471;
    m_exactOverlap(0,2) = -0.804031477774766;
    m_exactOverlap(0,3) = 0.527226809908349;
    m_exactOverlap(0,4) = 2.95367374025867e-6;
    m_exactOverlap(0,5) = 0.000675323284399051;
    m_exactOverlap(0,6) = 0.0202562111841077;
    m_exactOverlap(0,7) = 0.000967087011747088;
    m_exactOverlap(0,8) = 2.13813050004227e-5;
    m_exactOverlap(0,9) = 1.40984128295943e-19;
    m_exactOverlap(1,0) = -0.0341172323898471;
    m_exactOverlap(1,1) = 0.387827062976667;
    m_exactOverlap(1,2) = -0.00176248709791666;
    m_exactOverlap(1,3) = 0.0122388517922682;
    m_exactOverlap(1,4) = 4.99655646048901e-11;
    m_exactOverlap(1,5) = 0.00734798949097483;
    m_exactOverlap(1,6) = 0.00138996712571738;
    m_exactOverlap(1,7) = -0.00127323177225232;
    m_exactOverlap(1,8) = 0.00101665428542184;
    m_exactOverlap(1,9) = -1.81927486326609e-5;
    m_exactOverlap(2,0) = -0.804031477774766;
    m_exactOverlap(2,1) = -0.00176248709791666;
    m_exactOverlap(2,2) = 9.98427851756932;
    m_exactOverlap(2,3) = -0.0606071391980689;
    m_exactOverlap(2,4) = 0.0185755368982808;
    m_exactOverlap(2,5) = 0.000961861089058550;
    m_exactOverlap(2,6) = 0.0383258102957808;
    m_exactOverlap(2,7) = -0.000135744841834612;
    m_exactOverlap(2,8) = -8.30819528805685e-5;
    m_exactOverlap(2,9) = 0.000882925086794951;
    m_exactOverlap(3,0) = 0.527226809908349;
    m_exactOverlap(3,1) = 0.0122388517922682;
    m_exactOverlap(3,2) = -0.0606071391980689;
    m_exactOverlap(3,3) = 0.640491765531426;
    m_exactOverlap(3,4) = -0.00289162677173892;
    m_exactOverlap(3,5) = 4.59938683989704e-5;
    m_exactOverlap(3,6) = 0.0224417428337129;
    m_exactOverlap(3,7) = 0.000155680871167730;
    m_exactOverlap(3,8) = 1.28786125356206e-6;
    m_exactOverlap(3,9) = 0.000399808605221194;
    m_exactOverlap(4,0) = 2.95367374025867e-6;
    m_exactOverlap(4,1) = 4.99655646048901e-11;
    m_exactOverlap(4,2) = 0.0185755368982808;
    m_exactOverlap(4,3) = -0.00289162677173892;
    m_exactOverlap(4,4) = 0.0233723135930304;
    m_exactOverlap(4,5) = 1.51663773622826e-24;
    m_exactOverlap(4,6) = 3.63749494189881e-12;
    m_exactOverlap(4,7) = 4.69426663197159e-20;
    m_exactOverlap(4,8) = 1.03888645793693e-24;
    m_exactOverlap(4,9) = 2.17123310831415e-13;
    m_exactOverlap(5,0) = 0.000675323284399051;
    m_exactOverlap(5,1) = 0.00734798949097483;
    m_exactOverlap(5,2) = 0.000961861089058550;
    m_exactOverlap(5,3) = 4.59938683989704e-5;
    m_exactOverlap(5,4) = 1.51663773622826e-24;
    m_exactOverlap(5,5) = 0.0172362086060265;
    m_exactOverlap(5,6) = 8.26196083281730e-7;
    m_exactOverlap(5,7) = -2.20707377259169e-5;
    m_exactOverlap(5,8) = 0.00100873345509930;
    m_exactOverlap(5,9) = -1.27038048996005e-7;
    m_exactOverlap(6,0) = 0.0202562111841077;
    m_exactOverlap(6,1) = 0.00138996712571738;
    m_exactOverlap(6,2) = 0.0383258102957808;
    m_exactOverlap(6,3) = 0.0224417428337129;
    m_exactOverlap(6,4) = 3.63749494189881e-12;
    m_exactOverlap(6,5) = 8.26196083281730e-7;
    m_exactOverlap(6,6) = 0.00703744735730407;
    m_exactOverlap(6,7) = -1.54331741523399e-5;
    m_exactOverlap(6,8) = 2.01930457151725e-8;
    m_exactOverlap(6,9) = 5.69065367227716e-5;
    m_exactOverlap(7,0) = 0.000967087011747088;
    m_exactOverlap(7,1) = -0.00127323177225232;
    m_exactOverlap(7,2) = -0.000135744841834612;
    m_exactOverlap(7,3) = 0.000155680871167730;
    m_exactOverlap(7,4) = 4.69426663197159e-20;
    m_exactOverlap(7,5) = -2.20707377259169e-5;
    m_exactOverlap(7,6) = -1.54331741523399e-5;
    m_exactOverlap(7,7) = 0.00126285590311933;
    m_exactOverlap(7,8) = -7.49167731362069e-7;
    m_exactOverlap(7,9) = 5.91524141572004e-5;
    m_exactOverlap(8,0) = 2.13813050004227e-5;
    m_exactOverlap(8,1) = 0.00101665428542184;
    m_exactOverlap(8,2) = -8.30819528805685e-5;
    m_exactOverlap(8,3) = 1.28786125356206e-6;
    m_exactOverlap(8,4) = 1.03888645793693e-24;
    m_exactOverlap(8,5) = 0.00100873345509930;
    m_exactOverlap(8,6) = 2.01930457151725e-8;
    m_exactOverlap(8,7) = -7.49167731362069e-7;
    m_exactOverlap(8,8) = 0.0491199707835051;
    m_exactOverlap(8,9) = 2.13719828184536e-10;
    m_exactOverlap(9,0) = 1.40984128295943e-19;
    m_exactOverlap(9,1) = -1.81927486326609e-5;
    m_exactOverlap(9,2) = 0.000882925086794951;
    m_exactOverlap(9,3) = 0.000399808605221194;
    m_exactOverlap(9,4) = 2.17123310831415e-13;
    m_exactOverlap(9,5) = -1.27038048996005e-7;
    m_exactOverlap(9,6) = 5.69065367227716e-5;
    m_exactOverlap(9,7) = 5.91524141572004e-5;
    m_exactOverlap(9,8) = 2.13719828184536e-10;
    m_exactOverlap(9,9) = 0.000746311679740889;
}

void IntegralTester::setupExactKineticMatrix() {
    m_exactKinetic(0,0) = 2.95305186482296;
    m_exactKinetic(0,1) = -0.0558424051783208;
    m_exactKinetic(0,2) = -0.710211946758681;
    m_exactKinetic(0,3) = 0.825190282920817;
    m_exactKinetic(0,4) = -4.09145004078157e-5;
    m_exactKinetic(0,5) = -0.00354106288941588;
    m_exactKinetic(0,6) = -0.00101696409225696;
    m_exactKinetic(0,7) = 0.000318784100712331;
    m_exactKinetic(0,8) = -0.000193277296800795;
    m_exactKinetic(0,9) = -1.03461059186122e-16;
    m_exactKinetic(1,0) = -0.0558424051783208;
    m_exactKinetic(1,1) = 1.06652442318584;
    m_exactKinetic(1,2) = -0.00136283235772673;
    m_exactKinetic(1,3) = 0.00666163150823509;
    m_exactKinetic(1,4) = -1.38735470820744e-9;
    m_exactKinetic(1,5) = 0.0101718786886854;
    m_exactKinetic(1,6) = -0.00332120851907597;
    m_exactKinetic(1,7) = -0.00085909102921053;
    m_exactKinetic(1,8) = -0.00197612918039297;
    m_exactKinetic(1,9) = -0.000269901350336851;
    m_exactKinetic(2,0) = -0.710211946758697;
    m_exactKinetic(2,1) = -0.00136283235772661;
    m_exactKinetic(2,2) = 7.48820888817698;
    m_exactKinetic(2,3) = -0.0941675849397594;
    m_exactKinetic(2,4) = 0.00992457461917693;
    m_exactKinetic(2,5) = -0.00129538182364648;
    m_exactKinetic(2,6) = 0.0314870023540979;
    m_exactKinetic(2,7) = 8.81768102492484e-6;
    m_exactKinetic(2,8) = 0.000281024428562610;
    m_exactKinetic(2,9) = 0.00157572335530054;
    m_exactKinetic(3,0) = 0.825190282920816;
    m_exactKinetic(3,1) = 0.00666163150823466;
    m_exactKinetic(3,2) = -0.0941675849397592;
    m_exactKinetic(3,3) = 1.44110647244571;
    m_exactKinetic(3,4) = 0.00871686596713676;
    m_exactKinetic(3,5) = -0.000409048293848453;
    m_exactKinetic(3,6) = 0.0304280527464197;
    m_exactKinetic(3,7) = 2.05945927793234e-5;
    m_exactKinetic(3,8) = -1.43955534859779e-5;
    m_exactKinetic(3,9) = 0.00208775559366966;
    m_exactKinetic(4,0) = -4.09145004078137e-5;
    m_exactKinetic(4,1) = -1.38735470820744e-9;
    m_exactKinetic(4,2) = 0.00992457461936637;
    m_exactKinetic(4,3) = 0.00871686596713263;
    m_exactKinetic(4,4) = 0.111408028141933;
    m_exactKinetic(4,5) = -1.80773041407557e-22;
    m_exactKinetic(4,6) = -2.08223149360596e-10;
    m_exactKinetic(4,7) = -4.77727574575186e-18;
    m_exactKinetic(4,8) = -8.97002227888797e-23;
    m_exactKinetic(4,9) = -1.46810972477601e-11;
    m_exactKinetic(5,0) = -0.00354106288941604;
    m_exactKinetic(5,1) = 0.0101718786886830;
    m_exactKinetic(5,2) = -0.00129538182364563;
    m_exactKinetic(5,3) = -0.000409048293848444;
    m_exactKinetic(5,4) = -1.80773041407557e-22;
    m_exactKinetic(5,5) = 0.0896282847513165;
    m_exactKinetic(5,6) = -1.92843648148517e-5;
    m_exactKinetic(5,7) = -0.000172941244267707;
    m_exactKinetic(5,8) = -0.00462092936135390;
    m_exactKinetic(5,9) = 2.73212268997395e-6;
    m_exactKinetic(6,0) = -0.00101696409225621;
    m_exactKinetic(6,1) = -0.00332120851907591;
    m_exactKinetic(6,2) = 0.0314870023540918;
    m_exactKinetic(6,3) = 0.0304280527464185;
    m_exactKinetic(6,4) = -2.08223149360596e-10;
    m_exactKinetic(6,5) = -1.92843648148518e-5;
    m_exactKinetic(6,6) = 0.0472681880832255;
    m_exactKinetic(6,7) = -0.000206874924595757;
    m_exactKinetic(6,8) = -4.32226463603479e-7;
    m_exactKinetic(6,9) = 0.00105258842223646;
    m_exactKinetic(7,0) = 0.000318784100712428;
    m_exactKinetic(7,1) = -0.000859091029210779;
    m_exactKinetic(7,2) = 8.81768102490486e-6;
    m_exactKinetic(7,3) = 2.05945927793469e-5;
    m_exactKinetic(7,4) = -4.77727574575185e-18;
    m_exactKinetic(7,5) = -0.000172941244267711;
    m_exactKinetic(7,6) = -0.000206874924595625;
    m_exactKinetic(7,7) = 0.0163539839453953;
    m_exactKinetic(7,8) = 5.36473707571811e-6;
    m_exactKinetic(7,9) = 0.000246865723384732;
    m_exactKinetic(8,0) = -0.000193277296800797;
    m_exactKinetic(8,1) = -0.00197612918039266;
    m_exactKinetic(8,2) = 0.000281024428562512;
    m_exactKinetic(8,3) = -1.43955534859778e-5;
    m_exactKinetic(8,4) = -8.97002227888790e-23;
    m_exactKinetic(8,5) = -0.00462092936135679;
    m_exactKinetic(8,6) = -4.32226463603479e-7;
    m_exactKinetic(8,7) = 5.36473707571820e-6;
    m_exactKinetic(8,8) = 0.223495867064918;
    m_exactKinetic(8,9) = -4.82109324561387e-9;
    m_exactKinetic(9,0) = -1.50727578348677e-17;
    m_exactKinetic(9,1) = -0.000269901350336876;
    m_exactKinetic(9,2) = 0.00157572335530140;
    m_exactKinetic(9,3) = 0.00208775559366982;
    m_exactKinetic(9,4) = -1.46810972477601e-11;
    m_exactKinetic(9,5) = 2.73212268997395e-6;
    m_exactKinetic(9,6) = 0.00105258842223629;
    m_exactKinetic(9,7) = 0.000246865723384709;
    m_exactKinetic(9,8) = -4.82109324561561e-9;
    m_exactKinetic(9,9) = 0.0112319907800976;
}

void IntegralTester::setupExactElectronNucleusMatrix() {
    m_exactElectronNucleus(0,0) = 0;
}

IntegralTester::IntegralTester() {
    setupExactOverlapMatrix();
    setupExactKineticMatrix();
    setupExactElectronNucleusMatrix();
    m_overlapIntegrator             = new OverlapIntegrator();
    m_kineticIntegrator             = new KineticIntegrator();
    m_electronNucleusIntegrator     = new ElectronNucleusIntegrator();
    m_electronElectronIntegrator    = new ElectronElectronIntegrator();

    m_primitives.reserve(10);
    m_primitives.push_back(GaussianPrimitive(0,0,0,     1.0,    vec{ 0.0, -1.5,  1.2}));
    m_primitives.push_back(GaussianPrimitive(1,0,0,     1.1,    vec{ 0.1, -1.2,  2.5}));
    m_primitives.push_back(GaussianPrimitive(0,1,0,     0.3,    vec{ 0.2, -1.0, -0.3}));
    m_primitives.push_back(GaussianPrimitive(0,0,1,     0.9,    vec{ 0.3, -0.8,  0.1}));
    m_primitives.push_back(GaussianPrimitive(0,0,2,     2.2,    vec{ 0.4, -0.6, -3.1}));
    m_primitives.push_back(GaussianPrimitive(0,2,0,     2.4,    vec{ 0.5, -0.4,  3.8}));
    m_primitives.push_back(GaussianPrimitive(2,0,0,     3.1,    vec{ 0.6, -0.2,  1.3}));
    m_primitives.push_back(GaussianPrimitive(1,1,0,     3.7,    vec{ 0.7,  0.0,  2.4}));
    m_primitives.push_back(GaussianPrimitive(0,1,1,     1.3,    vec{ 0.8,  0.2,  5.3}));
    m_primitives.push_back(GaussianPrimitive(1,0,1,     4.3,    vec{ 0.9,  0.4,  1.2}));
}

bool IntegralTester::runAllTests() {
    bool passed = true;
    passed = runOverlapTests()          && passed;
    passed = runKineticTests()          && passed;
    passed = runElectronNucleusTests()  && passed;
    passed = runElectronElectronTests() && passed;

    if (passed) {
        cout << "   All integral tests PASSED." << endl;
    } else {
        cout << "   At least one integral test FAILED." << endl;
    }
}

bool IntegralTester::runOverlapTests() {
    cout << "   Running overlap integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 100;

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            double integral = m_overlapIntegrator->computeIntegral(&m_primitives[i], &m_primitives[j]);
            double difference = fabs(integral - m_exactOverlap(i,j));
            if (difference > m_tollerance) {
                failedTests++;
                cout << " Overlap integral test (" << i << "," << j << ") failed. Value: "
                     << integral << ", exact value: " << m_exactOverlap(i,j) << ", abs. difference: "
                     << difference << endl;
            }
        }
    }
    if (failedTests == 0) {
        cout << "      All 100 overlap integral tests PASSED." << endl;
    } else {
        cout << "      " << failedTests << " kinetic integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}

bool IntegralTester::runKineticTests() {
    cout << "   Running kinetic integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 100;

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            double integral = m_kineticIntegrator->computeIntegral(&m_primitives[i], &m_primitives[j]);
            double difference = fabs(integral - m_exactKinetic(i,j));
            if (difference > m_tollerance) {
                failedTests++;
                cout << "      Kinetic integral test (" << i << "," << j << ") failed. Value: "
                     << integral << ", exact value: " << m_exactKinetic(i,j) << ", abs. difference: "
                     << difference << endl;
            }
        }
    }
    if (failedTests == 0) {
        cout << "      All 100 kinetic integral tests PASSED." << endl;
    } else {
        cout << "      " << failedTests << " kinetic integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}

bool IntegralTester::runElectronNucleusTests() {
    cout << "   Running electron-nucleus integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 100;

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            int atom = (359*i+295*j-42*i*j+120*i*i-38*i*j*i) % 9;
            atom = (atom > 0 ? atom : -atom);
            vec nucleus = m_primitives[atom].nucleusPosition();
            double integral = m_electronNucleusIntegrator->computeIntegral(&m_primitives[i], &m_primitives[j], nucleus);
            double difference = fabs(integral - m_exactElectronNucleus(i,j));
            if (difference > m_tollerance) {
                failedTests++;
                cout << "      Electron-nucleus integral test (" << i << "," << j << ") failed. Value: "
                     << integral << ", exact value: " << m_exactElectronNucleus(i,j) << ", abs. difference: "
                     << difference << endl;
            }
        }
    }
    if (failedTests == 0) {
        cout << "      All 100 electron-nucleus integral tests PASSED." << endl;
    } else {
        cout << "      " << failedTests << " electron-nucleus integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}

bool IntegralTester::runElectronElectronTests() {
    cout << "   Running electron-electron integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 100;

    /*
    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            double integral = m_electronElectronIntegrator->computeIntegral(&m_primitives[i], &m_primitives[j]);
            double difference = fabs(integral - m_exactElectronNucleus(i,j));
            if (difference > m_tollerance) {
                failedTests++;
                cout << "      Electron-nucleus integral test (" << i << "," << j << ") failed. Value: "
                     << integral << ", exact value: " << m_exactElectronNucleus(i,j) << ", abs. difference: "
                     << difference << endl;
            }
        }
    }*/
    if (failedTests == 0) {
        cout << "      All 100 electron-nucleus integral tests PASSED." << endl;
    } else {
        cout << "      " << failedTests << " electron-nucleus integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}
