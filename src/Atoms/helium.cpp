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
    } else if (basisName == "STO-2G") {
        basis_STO2G();
    } else if (basisName == "STO-3G") {
        basis_STO3G();
    } else if (basisName == "STO-6G") {
        basis_STO6G();
    } else if (basisName == "toy") {
        basis_toy();
    } else if (basisName == "thijssen") {
        basis_thijssen();
    } else {
        cout << "Unknown basis: " << basisName << endl;
        cout << "Currently known basis sets for Helium: " << endl;
        cout << " * 3-21G"              << endl;
        cout << " * 6-311+G**"          << endl;
        cout << " * 6-311G(2df,2pd)"    << endl;
        cout << " * STO-6G"             << endl;
        cout << " * Thijssen"           << endl;
        exit(1);
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

void Helium::basis_STO2G() {
    m_info = "Helium : STO-2G";
    setNumberOfOrbitals(1);
    create_S2(2.4328790, 0.4330510, 0.4301280, 0.6789140);
}

void Helium::basis_STO3G() {
    m_info = "Helium : STO-3G";
    setNumberOfOrbitals(1);
    create_S3(6.36242139, 1.15892300, 0.31364979, 0.15432897, 0.53532814,0.44463454);
}

void Helium::basis_STO6G() {
    m_info = "Helium : STO-6G";
    setNumberOfOrbitals(1);
    create_S6(65.98456824,12.09819836,3.384639924,1.162715163,0.451516322,0.185959356,0.00916359628,0.04936149294,0.1685383049,0.3705627997,0.4164915298,0.1303340841);
}

void Helium::basis_toy() {
    m_info = "Toy example";
    setNumberOfOrbitals(2);
    //create_S1(3.8, 1.0);
    create_S1(0.5, 1.0);
    create_S1(3.8, 1.0);
    //create_S2(98.1243000, 3.3188300, 0.0287452, 0.8376350);

}

void Helium::basis_thijssen() {
    m_info = "Thijssen example";
    setNumberOfOrbitals(4);
    create_S1(0.298073,  1.0);
    create_S1(1.242567,  1.0);
    create_S1(5.782948,  1.0);
    create_S1(38.474970, 1.0);
}








