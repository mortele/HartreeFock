#include "hydrogen_631gss.h"

Hydrogen_631Gss::Hydrogen_631Gss(arma::vec position) :
        //             orbitals, electrons, charge
        Atom(position, 5,        1,         1) {

    m_info = "Hydrogen : 6-31G**";

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



/*
 * $basis
 * *
 * h   6-31G**
 * *
 *     3  s
 *      18.7311370              0.03349460
 *       2.8253937              0.23472695
 *       0.6401217              0.81375733
 *     1  s
 *       0.1612778              1.0000000
 *     1  p
 *       1.1000000              1.0000000
 * *
 * $end
 */
