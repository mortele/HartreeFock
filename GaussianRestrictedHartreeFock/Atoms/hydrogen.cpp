#include "hydrogen.h"
#include <iostream>

using std::cout;
using std::endl;

void Hydrogen::basis_321G() {
    /*   h   3-21G
     *   *
     *      2  s
     *        5.4471780              0.1562850
     *        0.8245470              0.9046910
     *      1  s
     *        0.1831920              1.0000000
     */
    m_info = "Hydrogen : 3-21G";
    setNumberOfOrbitals(2);

    // s orbital (2 primitives)
    GaussianPrimitive* primitive1 = new GaussianPrimitive(0, 0, 0, 5.4471780, m_position, 0.1562850);
    GaussianPrimitive* primitive2 = new GaussianPrimitive(0, 0, 0, 0.8245470, m_position, 0.9046910);

    // s orbital (1 primitive)
    GaussianPrimitive* primitive3 = new GaussianPrimitive(0, 0, 0, 0.1831920, m_position, 1.0000000);

    ContractedGaussian* contracted1 = new ContractedGaussian();
    contracted1->addPrimitive(primitive1, primitive1->getCoefficient());
    contracted1->addPrimitive(primitive2, primitive2->getCoefficient());

    ContractedGaussian* contracted2 = new ContractedGaussian();
    contracted2->addPrimitive(primitive3, primitive3->getCoefficient());

    m_contractedGaussians.push_back(contracted1);
    m_contractedGaussians.push_back(contracted2);
}

void Hydrogen::basis_321Gpp() {
    /*   h   3-21++G
     *   *
     *      2  s
     *         5.4471780              0.1562850
     *         0.8245470              0.9046910
     *      1  s
     *          0.1831920             1.0000000
     *      1  s
     *          0.0360000             1.0000000
     */
    m_info = "Hydrogen : 3-21++G";
    setNumberOfOrbitals(3);

    // s orbital (2 primitives)
    GaussianPrimitive* primitive1 = new GaussianPrimitive(0, 0, 0, 5.4471780, m_position, 0.1562850);
    GaussianPrimitive* primitive2 = new GaussianPrimitive(0, 0, 0, 0.8245470, m_position, 0.9046910);

    // s orbital (1 primitive)
    GaussianPrimitive* primitive3 = new GaussianPrimitive(0, 0, 0, 0.1831920, m_position, 1.0000000);

    // s orbital (1 primitive)
    GaussianPrimitive* primitive4 = new GaussianPrimitive(0, 0, 0, 0.0360000, m_position, 1.0000000);

    ContractedGaussian* contracted1 = new ContractedGaussian();
    contracted1->addPrimitive(primitive1, primitive1->getCoefficient());
    contracted1->addPrimitive(primitive2, primitive2->getCoefficient());

    ContractedGaussian* contracted2 = new ContractedGaussian();
    contracted2->addPrimitive(primitive3, primitive3->getCoefficient());

    ContractedGaussian* contracted3 = new ContractedGaussian();
    contracted3->addPrimitive(primitive4, primitive4->getCoefficient());

    m_contractedGaussians.push_back(contracted1);
    m_contractedGaussians.push_back(contracted2);
    m_contractedGaussians.push_back(contracted3);
}

void Hydrogen::basis_631G() {
    /*  h   6-31G
    *   *
    *     3  s
    *      18.7311370              0.03349460
    *       2.8253937              0.23472695
    *       0.6401217              0.81375733
    *     1  s
    *       0.1612778              1.0000000
    */
    m_info = "Hydrogen : 6-31G";
    setNumberOfOrbitals(2);

    // s orbital (3 primitives)
    GaussianPrimitive* primitive1 = new GaussianPrimitive(0, 0, 0, 18.7311370, m_position, 0.03349460);
    GaussianPrimitive* primitive2 = new GaussianPrimitive(0, 0, 0,  2.8253937, m_position, 0.23472695);
    GaussianPrimitive* primitive3 = new GaussianPrimitive(0, 0, 0,  0.6401217, m_position, 0.81375733);

    // s orbital (1 primitive)
    GaussianPrimitive* primitive4 = new GaussianPrimitive(0, 0, 0, 0.1612778, m_position, 1.0000000);

    ContractedGaussian* contracted1 = new ContractedGaussian();
    contracted1->addPrimitive(primitive1, primitive1->getCoefficient());
    contracted1->addPrimitive(primitive2, primitive2->getCoefficient());
    contracted1->addPrimitive(primitive3, primitive3->getCoefficient());

    ContractedGaussian* contracted2 = new ContractedGaussian();
    contracted2->addPrimitive(primitive4, primitive4->getCoefficient());

    m_contractedGaussians.push_back(contracted1);
    m_contractedGaussians.push_back(contracted2);
}

void Hydrogen::basis_631Gss() {
    /*  h   6-31G**
    *   *
    *     3  s
    *      18.7311370              0.03349460
    *       2.8253937              0.23472695
    *       0.6401217              0.81375733
    *     1  s
    *       0.1612778              1.0000000
    *     1  p
    *       1.1000000              1.0000000
    */
    m_info = "Hydrogen : 6-31G**";
    setNumberOfOrbitals(5);

    // s orbital (3 primitives) (1 contracted)
    GaussianPrimitive*  primitive1 = new GaussianPrimitive(0, 0, 0, 18.7311370, m_position, 0.03349460);
    GaussianPrimitive*  primitive2 = new GaussianPrimitive(0, 0, 0,  2.8253937, m_position, 0.23472695);
    GaussianPrimitive*  primitive3 = new GaussianPrimitive(0, 0, 0,  0.6401217, m_position, 0.81375733);
    ContractedGaussian* contracted = new ContractedGaussian();
    contracted->addPrimitive(primitive1, primitive1->getCoefficient());
    contracted->addPrimitive(primitive2, primitive2->getCoefficient());
    contracted->addPrimitive(primitive3, primitive3->getCoefficient());
    m_contractedGaussians.push_back(contracted);

    // s orbital (1 primitive) (1 contracted)
    primitive1 = new GaussianPrimitive(0, 0, 0, 0.1612778, m_position, 1.0000000);
    contracted = new ContractedGaussian();
    contracted->addPrimitive(primitive1, primitive1->getCoefficient());
    m_contractedGaussians.push_back(contracted);

    // p orbital (1 primitive) (3 contracteds)
    primitive1 = new GaussianPrimitive(1, 0, 0, 1.1000000, m_position, 1.0000000);
    primitive2 = new GaussianPrimitive(0, 1, 0, 1.1000000, m_position, 1.0000000);
    primitive3 = new GaussianPrimitive(0, 0, 1, 1.1000000, m_position, 1.0000000);
    contracted = new ContractedGaussian();
    contracted->addPrimitive(primitive1, primitive1->getCoefficient());
    m_contractedGaussians.push_back(contracted);

    contracted = new ContractedGaussian();
    contracted->addPrimitive(primitive2, primitive2->getCoefficient());
    m_contractedGaussians.push_back(contracted);

    contracted = new ContractedGaussian();
    contracted->addPrimitive(primitive3, primitive3->getCoefficient());
    m_contractedGaussians.push_back(contracted);
}

Hydrogen::Hydrogen(std::string basisName, arma::vec position) :
    Atom(position, 1, 1.0) {

    if (basisName == "3-21G") {
        basis_321G();
    } else if (basisName == "3-21G++") {
        basis_321Gpp();
    } else if (basisName == "6-31G") {
        basis_631G();
    } else if (basisName == "6-31G**") {
        basis_631Gss();
    } else {
        cout << "Unknown basis: " << basisName << endl;
        cout << "Currently known basis sets for Hydrogen: " << endl;
        cout << " * 3-21G"   << endl;
        cout << " * 3-21G++" << endl;
        cout << " * 6-31G"   << endl;
        cout << " * 6-31G**" << endl;
    }
}

