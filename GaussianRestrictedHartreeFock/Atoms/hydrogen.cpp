#include "hydrogen.h"
#include <iostream>

using std::cout;
using std::endl;

void Hydrogen::basis_321G() {
    /*   h   3-21G
     *   *
     *      2  s
     *        5.4471780              0.1562850
     *        0.8245470              0.9046910
     *      1  s
     *        0.1831920              1.0000000
     */
    m_info = "Hydrogen : 3-21G";
    setNumberOfOrbitals(2);

    create_S2(5.4471780, 0.8245470, 0.1562850, 0.9046910);
    create_S1(0.1831920, 1.0000000);
}

void Hydrogen::basis_321ppG() {
    /*   h   3-21++G
     *   *
     *      2  s
     *         5.4471780              0.1562850
     *         0.8245470              0.9046910
     *      1  s
     *          0.1831920             1.0000000
     *      1  s
     *          0.0360000             1.0000000
     */
    m_info = "Hydrogen : 3-21++G";
    setNumberOfOrbitals(3);

    create_S2(5.4471780, 0.8245470, 0.1562850, 0.9046910);
    create_S1(0.1831920, 1.0000000);
    create_S1(0.0360000, 1.0000000);
}

void Hydrogen::basis_631G() {
    /*  h   6-31G
    *   *
    *     3  s
    *      18.7311370              0.03349460
    *       2.8253937              0.23472695
    *       0.6401217              0.81375733
    *     1  s
    *       0.1612778              1.0000000
    */
    m_info = "Hydrogen : 6-31G";
    setNumberOfOrbitals(2);

    create_S3(18.7311370, 2.8253937, 0.6401217, 0.03349460, 0.23472695, 0.81375733);
    create_S1(0.1612778, 1.0000000);
}

void Hydrogen::basis_631Gss() {
    /*  h   6-31G**
    *   *
    *     3  s
    *      18.7311370              0.03349460
    *       2.8253937              0.23472695
    *       0.6401217              0.81375733
    *     1  s
    *       0.1612778              1.0000000
    *     1  p
    *       1.1000000              1.0000000
    */
    m_info = "Hydrogen : 6-31G**";
    setNumberOfOrbitals(3);

    create_S3(18.7311370, 2.8253937, 0.6401217, 0.03349460, 0.23472695, 0.81375733);
    create_S1(0.1612778, 1.0000000);
    create_P1(1.1000000, 1.0000000);
}

void Hydrogen::basis_631ppGss() {
    /*  h   6-31++G**
    *   *
    *     3  s
    *      18.7311370              0.03349460
    *       2.8253937              0.23472695
    *       0.6401217              0.81375733
    *     1  s
    *       0.1612778              1.0000000
    *     1  s
    *       0.0360000              1.0000000
    *     1  p
    *       1.1000000              1.0000000
    */
    m_info = "Hydrogen : 6-31++G**";
    setNumberOfOrbitals(4);

    create_S3(18.7311370, 2.8253937, 0.6401217, 0.03349460, 0.23472695, 0.81375733);
    create_S1(0.1612778, 1.0000000);
    create_S1(0.0360000, 1.0000000);
    create_P1(1.1000000, 1.0000000);
}

void Hydrogen::basis_6311ppGss() {
    /*  h   6-311++G**
    *   *
    *     3  s
    *      33.8650000              0.0254938
    *       5.0947900              0.1903730
    *       1.1587900              0.8521610
    *     1  s
    *       0.3258400              1.0000000
    *     1  s
    *       0.1027410              1.0000000
    *     1  s
    *       0.0360000              1.0000000
    *     1  p
    *       0.7500000              1.0000000
    */
    m_info = "Hydrogen : 6-311++G**";
    setNumberOfOrbitals(5);

    create_S3(33.8650000, 5.0947900, 1.1587900, 0.0254938, 0.1903730, 0.8521610);
    create_S1(0.3258400, 1.0000000);
    create_S1(0.1027410, 1.0000000);
    create_S1(0.0360000, 1.0000000);
    create_P1(0.7500000, 1.0000000);
}

void Hydrogen::basis_6311ppG2d2p() {
    /*  h   6-311++G(2d,2p)
    *   *
    *     3  s
    *      33.8650000              0.0254938
    *       5.0947900              0.1903730
    *       1.1587900              0.8521610
    *     1  s
    *       0.3258400              1.0000000
    *     1  s
    *       0.1027410              1.0000000
    *     1  s
    *       0.0360000              1.0000000
    *     1  p
    *       1.5000000              1.0000000
    *     1  p
    *       0.3750000              1.0000000
    */
    m_info = "Hydrogen : 6-311++G(2d,2p)";
    setNumberOfOrbitals(6);

    create_S3(33.8650000, 5.0947900, 1.1587900, 0.0254938, 0.1903730, 0.8521610);
    create_S1(0.3258400, 1.0000000);
    create_S1(0.1027410, 1.0000000);
    create_S1(0.0360000, 1.0000000);
    create_P1(1.5000000, 1.0000000);
    create_P1(0.3750000, 1.0000000);
}

void Hydrogen::basis_augccpVQZ() {
    m_info = "Hydrogen : aug-cc-pVQZ";
    setNumberOfOrbitals(14);
    create_S3(82.64,12.41,2.824,0.002006,0.015343,0.075579);
    create_S1(0.7977,1.0);
    create_S1(0.2581,1.0);
    create_S1(0.08989,1.0);
    create_S1(0.02363,1.0);
    create_P1(2.292,1.0);
    create_P1(0.838,1.0);
    create_P1(0.292,1.0);
    create_P1(0.0848,1.0);
    create_D1(2.062,1.0);
    create_D1(0.662,1.0);
    create_D1(0.19,1.0);
    create_F1(1.397,1.0);
    create_F1(0.36,1.0);
}



Hydrogen::Hydrogen(std::string basisName, arma::vec position) :
    Atom(position, 1, 1.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if (basisName == "3-21++G") {
        basis_321ppG();
    } else if (basisName == "6-31G") {
        basis_631G();
    } else if (basisName == "6-31G**") {
        basis_631Gss();
    } else if (basisName == "6-31++G**") {
        basis_631ppGss();
    } else if (basisName == "6-311++G**") {
        basis_6311ppGss();
    } else if (basisName == "6-311++G(2d,2p)") {
        basis_6311ppG2d2p();
    } else if (basisName == "aug-cc-pVQZ") {
        basis_augccpVQZ();
    } else {
        cout << "Unknown basis: " << basisName << endl;
        cout << "Currently known basis sets for Hydrogen: " << endl;
        cout << " * 3-21G"              << endl;
        cout << " * 3-21++G"            << endl;
        cout << " * 6-31G"              << endl;
        cout << " * 6-31G**"            << endl;
        cout << " * 6-31++G**"          << endl;
        cout << " * 6-311++G**"         << endl;
        cout << " * 6-311++G(2d,2p)"    << endl;
        cout << " * aug-cc-pVQZ"        << endl;
        exit(1);
    }
}

