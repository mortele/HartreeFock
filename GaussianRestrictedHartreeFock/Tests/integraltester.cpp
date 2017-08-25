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

IntegralTester::IntegralTester() {
    setupExactOverlapMatrix();
    m_overlapIntegrator = new OverlapIntegrator();
}

bool IntegralTester::runAllTests() {
    runOverlapTests();
}

bool IntegralTester::runOverlapTests() {
    GaussianPrimitive primitive0 = GaussianPrimitive(0,0,0,     1.0,    vec{ 0.0, -1.5,  1.2});
    GaussianPrimitive primitive1 = GaussianPrimitive(1,0,0,     1.1,    vec{ 0.1, -1.2,  2.5});
    GaussianPrimitive primitive2 = GaussianPrimitive(0,1,0,     0.3,    vec{ 0.2, -1.0, -0.3});
    GaussianPrimitive primitive3 = GaussianPrimitive(0,0,1,     0.9,    vec{ 0.3, -0.8,  0.1});
    GaussianPrimitive primitive4 = GaussianPrimitive(0,0,2,     2.2,    vec{ 0.4, -0.6, -3.1});
    GaussianPrimitive primitive5 = GaussianPrimitive(0,2,0,     2.4,    vec{ 0.5, -0.4,  3.8});
    GaussianPrimitive primitive6 = GaussianPrimitive(2,0,0,     3.1,    vec{ 0.6, -0.2,  1.3});
    GaussianPrimitive primitive7 = GaussianPrimitive(1,1,0,     3.7,    vec{ 0.7,  0.0,  2.4});
    GaussianPrimitive primitive8 = GaussianPrimitive(0,1,1,     1.3,    vec{ 0.8,  0.2,  5.3});
    GaussianPrimitive primitive9 = GaussianPrimitive(1,0,1,     4.3,    vec{ 0.9,  0.4,  1.2});

    GaussianPrimitive* primitives = new GaussianPrimitive[10];
    primitives[0] = primitive0;
    primitives[1] = primitive1;
    primitives[2] = primitive2;
    primitives[3] = primitive3;
    primitives[4] = primitive4;
    primitives[5] = primitive5;
    primitives[6] = primitive6;
    primitives[7] = primitive7;
    primitives[8] = primitive8;
    primitives[9] = primitive9;

    int failedTests = 0;
    int totalTests  = 100;

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            double integral = m_overlapIntegrator->computeIntegral(&primitives[i], &primitives[j]);
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
        cout << "All 100 tests passed." << endl;
    } else {
        cout << failedTests << " tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}
