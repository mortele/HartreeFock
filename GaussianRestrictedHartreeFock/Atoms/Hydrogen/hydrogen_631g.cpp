#include "hydrogen_631g.h"
#include "Orbitals/gaussianprimitive.h"
#include "Orbitals/contractedgaussian.h"


Hydrogen_631G::Hydrogen_631G(arma::vec position) :
        //             orbitals, electrons
        Atom(position, 2,        1        ) {

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


/* $basis
 * *
 * h   6-31G
 * *
 *     3  s
 *      18.7311370              0.03349460
 *       2.8253937              0.23472695
 *       0.6401217              0.81375733
 *     1  s
 *       0.1612778              1.0000000
 * *
 * $end
 */
