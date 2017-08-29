#include "boron.h"

void Boron::basis_321G() {
    m_info = "Boron : 3-21G";
    setNumberOfOrbitals(5);
    create_S3(116.434,17.4314,3.68016,0.0629605,0.363304,0.697255);
    create_S2(2.28187,0.465248,-0.368662,1.19944);
    create_S1(0.124328,1.0);
    create_P2(2.28187,0.465248,0.231152,0.866764);
    create_P1(0.124328,1.0);
}

void Boron::basis_6311ppGss() {
    m_info = "Boron : 6-311++G**";
    setNumberOfOrbitals(10);
    create_S6(2858.89,428.14,97.5282,27.9693,8.21577,1.11278,0.00215375,0.0165823,0.082187,0.276618,0.629316,0.17377);
    create_S3(13.2415,3.00166,0.912856,0.117443,0.918002,-0.00265105);
    create_S1(0.315454,1.0);
    create_S1(0.0988563,1.0);
    create_S1(0.0315,1.0);
    create_P3(13.2415,3.00166,0.912856,0.04181,0.236575,0.816214);
    create_P1(0.315454,1.0);
    create_P1(0.0988563,1.0);
    create_P1(0.0315,1.0);
    create_D1(0.401,1.0);
}

Boron::Boron(std::string basisName, arma::vec position) :
    Atom(position, 5, 5.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if(basisName == "6-311++G**") {
        basis_6311ppGss();
    } else {
        std::cout << "Could not find basis " << basisName << " for Boron." << std::endl;
        exit(1);
    }
}
