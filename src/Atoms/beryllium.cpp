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

void Beryllium::basis_STO6G() {
    m_info = "Beryllium : STO-6G";
    setNumberOfOrbitals(3);
    create_S6(312.8704937,57.36446253,16.0485094,5.513096119,2.140896553,0.8817394283,0.00916359628,0.04936149294,0.1685383049,0.3705627997,0.4164915298,0.1303340841);
    create_S6(13.63324744,2.698375464,0.8386530829,0.3226600698,0.1401314882,0.0642325139,-0.01325278809,-0.04699171014,-0.03378537151,0.2502417861,0.5951172526,0.2407061763);
    create_P6(13.63324744,2.698375464,0.8386530829,0.3226600698,0.1401314882,0.0642325139,0.0037596966,0.0376793698,0.1738967435,0.4180364347,0.4258595477,0.1017082955);
}

Beryllium::Beryllium(std::string basisName, arma::vec position) :
    Atom(position, 4, 4.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if(basisName == "6-311++G**") {
        basis_631ppGss();
    } else if (basisName == "STO-6G") {
        basis_STO6G();
    } else {
        std::cout << "Could not find basis " << basisName << " for Beryllium." << std::endl;
        exit(1);
    }

}
