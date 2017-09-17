#include "fluorine.h"


void Fluorine::basis_321G() {
    m_info = "Fluorine : 3-21G";
    setNumberOfOrbitals(5);
    create_S3(413.801,62.2446,13.434,0.0585483,0.349308,0.709632);
    create_S2(9.77759,2.08617,-0.407327,1.22314);
    create_S1(0.482383,1.0);
    create_P2(9.77759,2.08617,0.24668,0.852321);
    create_P1(0.482383,1.0);
}

void Fluorine::basis_631ppGss() {
    m_info = "Fluorine : 6-311++G**";
    setNumberOfOrbitals(10);
    create_S6(11427.1,1722.35,395.746,115.139,33.6026,4.91901,0.00180093,0.0137419,0.0681334,0.233325,0.589086,0.299505);
    create_S3(55.4441,12.6323,3.71756,0.114536,0.920512,-0.00337804);
    create_S1(1.16545,1.0);
    create_S1(0.321892,1.0);
    create_S1(0.1076,1.0);
    create_P3(55.4441,12.6323,3.71756,0.0354609,0.237451,0.820458);
    create_P1(1.16545,1.0);
    create_P1(0.321892,1.0);
    create_P1(0.1076,1.0);
    create_D1(1.75,1.0);
}

Fluorine::Fluorine(std::string basisName, arma::vec position) :
    Atom(position, 9, 9.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if(basisName == "6-311++G**") {
        basis_631ppGss();
    } else {
        std::cout << "Could not find basis " << basisName << " for Flourine." << std::endl;
        exit(1);
    }
}
