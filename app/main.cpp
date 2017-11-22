//#define ARMA_MAT_PREALLOC 16
//#define ARMA_EXTRA_DEBUG
//#define ARMA_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <armadillo>
#include <xc.h>
#include "examples.h"

#include "system.h"
#include "Solvers/restricteddft.h"
#include "Solvers/restrictedhartreefock.h"
#include "Orbitals/gaussianprimitive.h"
#include "Orbitals/contractedgaussian.h"
#include "Atoms/hydrogen.h"
#include "Atoms/helium.h"
#include "Atoms/beryllium.h"
#include "Atoms/carbon.h"
#include "Integrators/numericalintegrator.h"
#include "ExchangeCorrelationFunctionals/localdensityapproximation.h"

using std::cout;
using std::endl;
using std::setprecision;


int main(int, char**) {
    Examples::firstExample();
    return 0;


    /*
    System* system = new System();
    system->addAtom(new Helium("toy", arma::vec{0,0,0}));
    RestrictedDFT* DFT = new RestrictedDFT(system);
    DFT->setFunctional("LDA");
    DFT->solve(1e-5,1);
    LocalDensityApproximation LDA(system, &DFT->m_densityMatrix);
    double rho_ = 1.52708108890903;
    double ans = -1.29877680133385;
    double ans2 = -0.11477000207871;
    double ans3 = -1.41354680341256;
    double mine =  rho_ * LDA.evaluateEnergy(rho_);
    //cout << setprecision(14) << mine << endl;
    //cout << fabs(mine-ans3) << endl;

    rho_ = 0.98333768675191;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ContractedGaussian* Gi = system->getBasis().at(i);
            ContractedGaussian* Gj = system->getBasis().at(j);
            double x =  0.0;
            double y =  0.1;
            double z = -0.1;
            double iPhi = Gi->evaluate(x,y,z);
            double jPhi = Gj->evaluate(x,y,z);

            cout << setprecision(14) << iPhi * jPhi * LDA.evaluatePotential(rho_) << " ";
        }
        cout << endl;
    }

    xc_func_type c;
    xc_func_type x;
    double rho[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    double sigma[5] = {0,0,0,0,0};
    double ex[5];
    double ec[5];
    double exc[5];
    int i, vmajor, vminor, vmicro, func_id = 1;

    xc_func_init(&x, XC_LDA_X,          XC_UNPOLARIZED);
    xc_func_init(&c, XC_LDA_C_VWN_4,  XC_UNPOLARIZED);

    xc_lda_exc(&x, 5, rho, ex);
    xc_lda_exc(&c, 5, rho, ec);

    double exc_homemade[5];
    for (int i = 0; i < 5; i++) {
        exc[i] = ex[i] + ec[i];
        //exc_homemade[i] = LDA.evaluatePotential(rho[i]);
        exc_homemade[i] = LDA.evaluateEnergy(rho[i]);
    }

    for(i=0; i<5; i+=1) {
        printf("%20.15f %20.15f %20.15f %20.15g\n", rho[i], exc[i], exc_homemade[i], fabs(exc[i]-exc_homemade[i]));
    }

    xc_func_end(&x);
    xc_func_end(&c);



    return 0;

    */
    //Examples::SingleAtom(4,"STO-6G",4,1e4,1e-8,"cuspTest");
    //Examples::H2();
    //return 0;
    System*         system = new System();

    //Helium*         helium      = new Helium    ("3-21G", arma::vec{0,0,0});
    //Helium*         helium      = new Helium    ("6-311+G**", arma::vec{0,0,0});
    //Helium*         helium      = new Helium    ("6-311G(2df,2pd)", arma::vec{0,0,0});
    Helium*         helium      = new Helium    ("toy", arma::vec{0,0,0});
    Hydrogen*       hydrogen1   = new Hydrogen  ("3-21G", arma::vec{0,0,0});
    Hydrogen*       hydrogen2   = new Hydrogen  ("3-21G", arma::vec{0,0,1.4});
    Beryllium*      beryllium   = new Beryllium ("3-21G", arma::vec{0,0,0});
    //Beryllium*      beryllium   = new Beryllium ("6-311++G**", arma::vec{0,0,0});
    Carbon*         carbon      = new Carbon    ("3-21G", arma::vec{0,0,0});

    //system->addAtom(hydrogen1);
    //system->addAtom(hydrogen2);
    system->addAtom(helium);
    //system->addAtom(beryllium);
    //system->addAtom(carbon);

    RestrictedHartreeFock*  rhf     = new RestrictedHartreeFock(system);
    RestrictedDFT*          rdft    = new RestrictedDFT(system);
    rdft->setFunctional("LDA");
    rdft->solve(1e-8,1);
    for (int i = 0; i < system->getNumberOfBasisFunctions(); i++) {
        cout << setprecision(15) << fabs(rdft->m_numericalIntegrator->testIntegral(i)-1) << endl;
    }
    //cout << setprecision(15) << "" << rdft->m_coefficientMatrix(0,0) << endl;
    //cout << setprecision(15) << "" << rdft->m_coefficientMatrix(1,0) << endl;

    // HF
    //rdft->m_coefficientMatrix(0,0) = 0.300859;
    //rdft->m_coefficientMatrix(1,0) = 0.811650;

    // DFT
    //rdft->m_coefficientMatrix(0,0) = 0.295500;
    //rdft->m_coefficientMatrix(1,0) = 0.815618;

    // 3-21G
    //rdft->m_coefficientMatrix(0,0) = 0.44176629832108;
    //rdft->m_coefficientMatrix(1,0) = 0.67192440307208;

    // Toy basis @psi4
    rdft->m_coefficientMatrix(0,0) = -0.29548391472761;
    rdft->m_coefficientMatrix(1,0) = -0.81563000823352;

    // 3-21G test psi4
    //rdft->m_coefficientMatrix(0,0) =  0.44741355811568;
    //rdft->m_coefficientMatrix(1,0) =  0.66682729446697;

    rdft->m_smoothing = false;
    rdft->computeDensityMatrix();
    arma::mat& S = rdft->m_overlapMatrix;
    arma::mat& P = rdft->m_densityMatrix;
    arma::mat& C = rdft->m_coefficientMatrix;
    //cout << C << endl;
    arma::mat& A = rdft->m_transformationMatrix;
    //rdft->m_densityMatrix = 2*C*C.t();
    GaussianPrimitive* pr1 = system->getBasis().at(0)->getPrimitives().at(0);
    GaussianPrimitive* pr2 = system->getBasis().at(1)->getPrimitives().at(0);
    //cout << *pr1 << endl;
    //cout << *pr2 << endl;

    //// !!!!!!!
    //// !!!!!!!
    //// https://github.com/psi4/psi4numpy/blob/master/Tutorials/04_Density_Functional_Theory/4b_LDA_kernel.ipynb
    //// !!!!!!!
    //// !!!!!!!

    return 0;

}

