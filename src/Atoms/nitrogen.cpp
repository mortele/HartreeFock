#include "nitrogen.h"

void Nitrogen::basis_321G() {
    m_info = "Nitrogen : 3-21G";
    setNumberOfOrbitals(5);
    create_S3(242.766,36.4851,7.81449,0.0598657,0.352955,0.706513);
    create_S2(5.42522,1.14915,-0.413301,1.22442);
    create_S1(0.283205,1.0);
    create_P2(5.42522,1.14915,0.237972,0.858953);
    create_P1(0.283205,1.0);
}

void Nitrogen::basis_631ppGss() {
    m_info = "Nitrogen : 6-311++G**";
    setNumberOfOrbitals(10);
    create_S6(6293.48,949.044,218.776,63.6916,18.8282,2.72023,0.00196979,0.0149613,0.0735006,0.248937,0.60246,0.256202);
    create_S3(30.6331,7.02614,2.11205,0.111906,0.921666,-0.00256919);
    create_S1(0.684009,1.0);
    create_S1(0.200878,1.0);
    create_S1(0.0639,1.0);
    create_P3(30.6331,7.02614,2.11205,0.0383119,0.237403,0.817592);
    create_P1(0.684009,1.0);
    create_P1(0.200878,1.0);
    create_P1(0.0639,1.0);
    create_D1(0.913,1.0);
}

Nitrogen::Nitrogen(std::string basisName, arma::vec position) :
    Atom(position, 7, 7.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if(basisName == "6-311++G**") {
        basis_631ppGss();
    } else {
        std::cout << "Could not find basis " << basisName << " for Nitrogen." << std::endl;
        exit(1);
    }

}
