#include "carbon.h"

void Carbon::basis_321G() {
    m_info = "Carbon : 3-21G";
    setNumberOfOrbitals(5);
    create_S3(172.256,25.9109,5.53335,0.0617669,0.358794,0.700713);
    create_S2(3.66498,0.770545,-0.395897,1.21584);
    create_S1(0.195857,1.0);
    create_P2(3.66498,0.770545,0.23646,0.860619);
    create_P1(0.195857,1.0);
}

void Carbon::basis_6311ppGss() {
    m_info = "Carbon : 6-311++G**";
    setNumberOfOrbitals(10);
    create_S6(4563.24,682.024,154.973,44.4553,13.029,1.82773,0.00196665,0.0152306,0.0761269,0.260801,0.616462,0.221006);
    create_S3(20.9642,4.80331,1.45933,0.11466,0.919999,-0.00303068);
    create_S1(0.483456,1.0);
    create_S1(0.145585,1.0);
    create_S1(0.0438,1.0);
    create_P3(20.9642,4.80331,1.45933,0.0402487,0.237594,0.815854);
    create_P1(0.483456,1.0);
    create_P1(0.145585,1.0);
    create_P1(0.0438,1.0);
    create_D1(0.626,1.0);
}

Carbon::Carbon(std::string basisName, arma::vec position) :
    Atom(position, 6, 6.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if(basisName == "6-311++G**") {
        basis_6311ppGss();
    } else {
        std::cout << "Could not find basis " << basisName << " for Carbon." << std::endl;
        exit(1);
    }
}
