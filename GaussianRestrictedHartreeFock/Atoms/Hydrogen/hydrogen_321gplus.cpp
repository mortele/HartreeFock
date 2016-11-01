#include "Atoms/Hydrogen/hydrogen_321Gplus.h"
#include "Orbitals/gaussianprimitive.h"
#include "Orbitals/contractedgaussian.h"


Hydrogen_321Gplus::Hydrogen_321Gplus(arma::vec position) :
        //             orbitals, electrons, charge
        Atom(position, 2,        1,         1) {

    m_info = "Hydrogen : 3-21++G";

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


/* $basis
 * *
 * h   3-21++G
 * *
 *     2  s
 *       5.4471780              0.1562850
 *       0.8245470              0.9046910
 *     1  s
 *       0.1831920              1.0000000
 *     1  s
 *       0.0360000              1.0000000
 * *
 * $end
 */
