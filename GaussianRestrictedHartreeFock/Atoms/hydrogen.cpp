#include "hydrogen.h"
#include "Orbitals/gaussianprimitive.h"
#include "Orbitals/contractedgaussian.h"


Hydrogen::Hydrogen(arma::vec position) :
        //             orbitals, electrons
        Atom(position, 2,        1        ) {

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


/*   $basis
 *   *
 *   h   3-21G
 *   *
 *       2  s
 *         5.4471780              0.1562850
 *         0.8245470              0.9046910
 *       1  s
 *         0.1831920              1.0000000
 *   *
 *   $end
 */
