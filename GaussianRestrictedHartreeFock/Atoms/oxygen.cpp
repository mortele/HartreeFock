#include "oxygen.h"
#include <iostream>

using std::cout;
using std::endl;

void Oxygen::basis_321G() {
   /*   o   3-21G
    *   *
    *       3  s
    *       322.0370000              0.0592394
    *        48.4308000              0.3515000
    *        10.4206000              0.7076580
    *       2  s
    *         7.4029400             -0.4044530
    *         1.5762000              1.2215600
    *       1  s
    *         0.3736840              1.0000000
    *       2  p
    *         7.4029400              0.2445860
    *         1.5762000              0.8539550
    *       1  p
    *         0.3736840              1.0000000
    */
    m_info = "Oxygen : 3-21G";
    setNumberOfOrbitals(5);

    create_S3(322.0370000, 48.4308000, 10.4206000, 0.0592394, 0.3515000, 0.7076580);
    create_S2(7.4029400, 1.5762000, -0.4044530, 1.2215600);
    create_S1(0.3736840, 1.0000000);
    create_P2(7.4029400, 1.5762000, 0.2445860, 0.8539550);
    create_P1(0.3736840, 1.0000000);
}

Oxygen::Oxygen(std::string basisName, arma::vec position) :
        Atom(position, 8, 8.0) {
    if (basisName == "3-21G") {
        basis_321G();
    /*} else if (basisName == "3-21++G") {
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
        basis_6311ppG2d2p();*/
    } else {
        cout << "Unknown basis: " << basisName << endl;
        cout << "Currently known basis sets for Oxygen: " << endl;
        cout << " * 3-21G"              << endl;
        /*cout << " * 3-21++G"            << endl;
        cout << " * 6-31G"              << endl;
        cout << " * 6-31G**"            << endl;
        cout << " * 6-31++G**"          << endl;
        cout << " * 6-311++G**"         << endl;
        cout << " * 6-311++G(2d,2p)"    << endl;*/
    }
}
