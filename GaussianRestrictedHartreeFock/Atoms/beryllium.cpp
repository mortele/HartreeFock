#include "beryllium.h"

void Beryllium::basis_321G() {
    m_info = "Beryllium : 3-21G";
    setNumberOfOrbitals(5);
    create_S3(71.8876,10.7289,2.22205,0.0644263,0.366096,0.695934);
    create_S2(1.29548,0.268881,-0.421064,1.22407);
    create_S1(0.07735,1.0);
    create_P2(1.29548,0.268881,0.205132,0.882528);
    create_P1(0.07735,1.0);
}

void Beryllium::basis_631ppGss() {
    m_info = "Beryllium : 6-311++G**";
    setNumberOfOrbitals(10);
    create_S6(1682.8,251.715,57.4116,16.5171,4.85364,0.626863,0.00228574,0.0175938,0.0863315,0.281835,0.640594,0.144467);
    create_S3(8.30938,1.74075,0.485816,0.108621,0.927301,-0.00297169);
    create_S1(0.163613,1.0);
    create_S1(0.0567285,1.0);
    create_S1(0.0207,1.0);
    create_P3(8.30938,1.74075,0.485816,0.0361344,0.216958,0.841839);
    create_P1(0.163613,1.0);
    create_P1(0.0567285,1.0);
    create_P1(0.0207,1.0);
    create_D1(0.255,1.0);
}

Beryllium::Beryllium(std::string basisName, arma::vec position) :
    Atom(position, 4, 4.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if(basisName == "6-311++G**") {
        basis_631ppGss();
    } else {
        std::cout << "Could not find basis " << basisName << " for Beryllium." << std::endl;
        exit(1);
    }

}
