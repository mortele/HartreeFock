#include "helium.h"
#include <iostream>

using std::cout;
using std::endl;

Helium::Helium(std::string basisName, arma::vec position) :
    Atom(position, 2, 2.0) {
    if (basisName == "3-21G") {
        basis_321G();
    } else if (basisName == "6-311+G**") {
        basis_6311pGss();
    } else if (basisName == "6-311G(2df,2pd)") {
        basis_6311G2df2pd();
    }
    else {
        cout << "Unknown basis: " << basisName << endl;
        cout << "Currently known basis sets for Helium: " << endl;
        cout << " * 3-21G"              << endl;
        cout << " * 6-311+G**"          << endl;
        cout << " * 6-311G(2df,2pd)"    << endl;
    }
}

void Helium::basis_321G() {
    /*  he   3-21G
     *  *
     *      2  s
     *       13.6267000              0.1752300
     *        1.9993500              0.8934830
     *      1  s
     *        0.3829930              1.0000000
     */
    m_info = "Helium : 3-21G";
    setNumberOfOrbitals(2);

    create_S2(13.6267, 1.99935, 0.1752300, 0.8934830);
    create_S1(0.3829930, 1.0000000);
}

void Helium::basis_6311pGss() {
    /*
     * he   6-311+G**
     * *
     *     3  s
     *      98.1243000              0.0287452
     *      14.7689000              0.2080610
     *       3.3188300              0.8376350
     *     1  s
     *       0.8740470              1.0000000
     *     1  s
     *       0.2445640              1.0000000
     *     1  p
     *       0.7500000              1.0000000
     *
     */
    m_info = "Helium : 6-311+G**";
    setNumberOfOrbitals(4);
    create_S3(98.1243000, 14.7689000, 3.3188300, 0.0287452, 0.2080610, 0.8376350);
    create_S1(0.8740470, 1.0000000);
    create_S1(0.2445640, 1.0000000);
    create_P1(0.7500000, 1.0000000);
}

void Helium::basis_6311G2df2pd() {
    /*
     * he   6-311G(2df,2pd)
     * *
     *     3  s
     *      98.1243000              0.0287452
     *      14.7689000              0.2080610
     *       3.3188300              0.8376350
     *     1  s
     *       0.8740470              1.0000000
     *     1  s
     *       0.2445640              1.0000000
     *     1  p
     *       1.5000000              1.0000000
     *     1  p
     *       0.3750000              1.0000000
     *     1  d
     *       2.0000000              1.0000000
     */
    m_info = "Helium : 6-311G(2df,2pd)";
    setNumberOfOrbitals(6);
    create_S3(98.1243000, 14.7689000, 3.3188300, 0.0287452, 0.2080610 , 0.8376350);
    create_S1(0.8740470, 1.0000000);
    create_S1(0.2445640, 1.0000000);
    create_P1(1.5000000, 1.0000000);
    create_P1(0.3750000, 1.0000000);
    create_D1(2.0000000, 1.0000000);
}
