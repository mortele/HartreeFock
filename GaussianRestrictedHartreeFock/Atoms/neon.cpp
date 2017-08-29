#include "neon.h"

void Neon::basis_321G() {
    m_info = "Neon : 3-21G";
    setNumberOfOrbitals(5);
    create_S3(515.724,77.6538,16.8136,0.058143,0.347951,0.710714);
    create_S2(12.483,2.66451,-0.409922,1.22431);
    create_S1(0.60625,1.0);
    create_P2(12.483,2.66451,0.24746,0.851743);
    create_P1(0.60625,1.0);
}

void Neon::basis_631ppGss() {
    m_info = "Neon : 6-311++G**";
    setNumberOfOrbitals(10);
    create_S6(13995.7,2117.1,490.425,143.833,41.9265,6.15684,0.00183276,0.0138827,0.0680687,0.231328,0.58589,0.305883);
    create_S3(69.1211,15.835,4.67326,0.119149,0.917375,-0.00405839);
    create_S1(1.45756,1.0);
    create_S1(0.397057,1.0);
    create_S1(0.13,1.0);
    create_P3(69.1211,15.835,4.67326,0.0356574,0.239477,0.818461);
    create_P1(1.45756,1.0);
    create_P1(0.397057,1.0);
    create_P1(0.13,1.0);
    create_D1(2.304,1.0);
}

Neon::Neon(std::string basisName, arma::vec position) :
    Atom(position, 10, 10.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if(basisName == "6-311++G**") {
        basis_631ppGss();
    } else {
        std::cout << "Could not find basis " << basisName << " for Neon." << std::endl;
        exit(1);
    }
}
