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
    m_info = "Oxygen   : 3-21G";
    setNumberOfOrbitals(5);

    create_S3(322.0370000, 48.4308000, 10.4206000, 0.0592394, 0.3515000, 0.7076580);
    create_S2(7.4029400, 1.5762000, -0.4044530, 1.2215600);
    create_S1(0.3736840, 1.0000000);
    create_P2(7.4029400, 1.5762000, 0.2445860, 0.8539550);
    create_P1(0.3736840, 1.0000000);
}

void Oxygen::basis_631pGss() {
   /* o   6-31+G**
    * *
    *    6  s
    *   5484.6717000              0.0018311
    *    825.2349500              0.0139501
    *    188.0469600              0.0684451
    *     52.9645000              0.2327143
    *     16.8975700              0.4701930
    *      5.7996353              0.3585209
    *    3  s
    *     15.5396160             -0.1107775
    *      3.5999336             -0.1480263
    *      1.0137618              1.1307670
    *    1  s
    *      0.2700058              1.0000000
    *    1  s
    *      0.0845000              1.0000000
    *    3  p
    *     15.5396160              0.0708743
    *      3.5999336              0.3397528
    *      1.0137618              0.7271586
    *    1  p
    *      0.2700058              1.0000000
    *    1  p
    *      0.0845000              1.0000000
    *    1  d
    *      0.8000000              1.0000000
    */
    m_info = "Oxygen   : 6-31+G**";
    setNumberOfOrbitals(8);
    create_S6(5484.6717000, 825.2349500, 188.0469600, 52.9645000, 16.8975700, 5.7996353, 0.0018311, 0.0139501, 0.0684451, 0.2327143, 0.4701930, 0.3585209);
    create_S3(15.5396160, 3.5999336, 1.0137618, -0.1107775, -0.1480263, 1.1307670);
    create_S1(0.2700058, 1.0000000);
    create_S1(0.0845000, 1.0000000);
    create_P3(15.5396160, 3.5999336, 1.0137618, 0.0708743, 0.3397528, 0.7271586);
    create_P1(0.2700058, 1.0000000);
    create_P1(0.0845000, 1.0000000);
    create_D1(0.8000000, 1.0000000);
}

void Oxygen::basis_6311ppGss() {
   /* o   6-311++G**
    * *
    *     6  s
    *    8588.5000000              0.00189515
    *    1297.2300000              0.0143859
    *     299.2960000              0.0707320
    *      87.3771000              0.2400010
    *      25.6789000              0.5947970
    *       3.7400400              0.2808020
    *     3  s
    *      42.1175000              0.1138890
    *       9.6283700              0.9208110
    *       2.8533200             -0.00327447
    *     1  s
    *       0.9056610              1.0000000
    *     1  s
    *       0.2556110              1.0000000
    *     1  s
    *       0.0845000              1.0000000
    *     3  p
    *      42.1175000              0.0365114
    *       9.6283700              0.2371530
    *       2.8533200              0.8197020
    *     1  p
    *       0.9056610              1.0000000
    *     1  p
    *       0.2556110              1.0000000
    *     1  p
    *       0.0845000              1.0000000
    *     1  d
    *       1.2920000              1.0000000
    */
    m_info = "Oxygen   : 6-311++G**";
    setNumberOfOrbitals(10);

    create_S6(8588.5000000, 1297.2300000, 299.2960000, 87.3771000, 25.6789000, 3.7400400, 0.00189515, 0.0143859, 0.0707320, 0.2400010, 0.5947970, 0.2808020);
    create_S3(42.1175000, 9.6283700, 2.8533200,  0.1138890, 0.9208110, -0.00327447);
    create_S1(0.9056610, 1.0000000);
    create_S1(0.2556110, 1.0000000);
    create_S1(0.0845000, 1.0000000);
    create_P3(42.1175000, 9.6283700, 2.8533200, 0.0365114, 0.2371530, 0.8197020);
    create_P1(0.9056610, 1.0000000);
    create_P1(0.2556110, 1.0000000);
    create_P1(0.0845000, 1.0000000);
    create_D1(1.2920000, 1.0000000);
}

Oxygen::Oxygen(std::string basisName, arma::vec position) :
        Atom(position, 8, 8.0) {
    if (basisName == "3-21G") {
        basis_321G();
    } else if (basisName == "6-31+G**") {
        basis_631pGss();
    } else if (basisName == "6-311++G**") {
        basis_6311ppGss();
    } else {
        cout << "Unknown basis: " << basisName << endl;
        cout << "Currently known basis sets for Oxygen: " << endl;
        cout << " * 3-21G"              << endl;
        cout << " * 6-311++G**"         << endl;
    }
}
