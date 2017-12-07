#include "examples.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <boost/timer.hpp>
#include "system.h"
#include "Solvers/restrictedhartreefock.h"
#include "Solvers/unrestrictedhartreefock.h"
#include "Atoms/atom.h" // Includes all the specific atoms currently available.



using arma::vec;
using arma::zeros;
using std::cout;
using std::endl;
using std::cos;
using std::sin;



void Examples::firstExample() {
    System He;
    He.addAtom(new Helium("3-21G", vec{0,0,0}));
    RestrictedHartreeFock solver(&He);
    solver.solve();
}

void Examples::secondExample() {
    vec O  { 0.000, 0.000, 0.000};
    vec H1 {-1.430, 1.108, 0.000};
    vec H2 { 1.430, 1.108, 0.000};

    System H2O;
    H2O.addAtom(new Oxygen  ("6-311++G**", O));
    H2O.addAtom(new Hydrogen("6-311++G**", H1));
    H2O.addAtom(new Hydrogen("6-311++G**", H2));

    UnrestrictedHartreeFock solver(&H2O);
    solver.solve();
}

void Examples::methane() {
    // CH4
    // ========================================================================
    vec CH4_nucleus1  {0,            0,             0};
    vec CH4_nucleus2  { 1.67381799,  0.        , -1.18356805};
    vec CH4_nucleus3  {-1.67381799,  0.        , -1.18356805};
    vec CH4_nucleus4  { 0.        ,  1.67381799,  1.18356805};
    vec CH4_nucleus5  { 0.        , -1.67381799,  1.18356805};
    System* system_CH4_big   = new System(5);
    system_CH4_big->  addAtom(new Carbon  ("6-311++G**", CH4_nucleus1));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus2));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus3));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus4));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus5));
    UnrestrictedHartreeFock un_solver_CH4_big  (system_CH4_big);

    double un_CH4_big    = un_solver_CH4_big.  solve(1e-10,1e4);
    un_solver_CH4_big.dumpBasisToFile("methane");
}

void Examples::diberyllium() {
    vec n1  { 0.000, 0.000, 0.000};
    vec n2  { 4.630, 0.000, 0.000};

    System Be2;
    Be2.addAtom(new Beryllium("6-311++G**", n1));
    Be2.addAtom(new Beryllium("6-311++G**", n2));

    UnrestrictedHartreeFock solver(&Be2);
    solver.solve();
}

void Examples::Hm() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};

    System* system = new System(1);
    Hydrogen* Hm = new Hydrogen("6-311++G(2d,2p)", nucleus1);
    Hm->setNumberOfElectrons(2);
    system->addAtom(Hm);

    //UnrestrictedHartreeFock solver(system);
    RestrictedHartreeFock solver(system);
    solver.solve(1e-5, 1e4);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;


    // -13.25093 eV  UHF H-
    // -14.348   eV      H- (http://nist.gov/data/PDFfiles/jpcrd68.pdf)
    // -13.6     eV      H
}

void Examples::He() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};

    System* system = new System(1);
    system->addAtom(new Helium("6-311+G**", nucleus1));

    //UnrestrictedHartreeFock solver(system);
    RestrictedHartreeFock solver(system);
    solver.solve(1e-5, 1e4);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}

void Examples::HeHp() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};
    vec nucleus2  {0, 0, 1.4632};

    System* system = new System(2);
    system->addAtom(new Helium("6-311+G**", nucleus1));
    Hydrogen* Hm = new Hydrogen("6-311++G(2d,2p)", nucleus2);
    Hm->setNumberOfElectrons(0);
    system->addAtom(Hm);

    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);
    solver.solve(1e-5, 1e4);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}


void Examples::H2() {
    boost::timer t;
    vec nucleus1  {0, 0, 0};
    vec nucleus2  {0, 0, 1.4};

    System* system = new System(2);
    //system->addAtom(new Hydrogen("3-21G", nucleus1));
    //system->addAtom(new Hydrogen("3-21G", nucleus2));
    system->addAtom(new Hydrogen("6-311++G**", nucleus1));
    system->addAtom(new Hydrogen("6-311++G**", nucleus2));
    //system->addAtom(new Hydrogen("cc-pVQZ", nucleus1));
    //system->addAtom(new Hydrogen("cc-pVQZ", nucleus2));
    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);
    solver.solve(1e-10, 1e4);
    solver.dumpBasisToFile();
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}


void Examples::H20() {
    boost::timer t;

    vec nucleus1 { 0.000, 0.000, 0.000};
    vec nucleus2 {-1.430, 1.108, 0.000};
    vec nucleus3 { 1.430, 1.108, 0.000};

    System* system = new System(3);
    system->addAtom(new Oxygen  ("6-311++G**", nucleus1));
    system->addAtom(new Hydrogen("6-311++G**", nucleus2));
    system->addAtom(new Hydrogen("6-311++G**", nucleus3));

    UnrestrictedHartreeFock solver(system);
    //RestrictedHartreeFock solver(system);

    solver.solve(1e-12, 1e3);
    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}

void Examples::ValidationTableEnergy() {
    boost::timer t;

    // H2
    // ========================================================================
    vec H2_nucleus1  {0, 0, 0};
    vec H2_nucleus2  {0, 0, 1.4};
    System* system_H2_small = new System(2);
    System* system_H2_big   = new System(2);
    system_H2_small->addAtom(new Hydrogen("3-21G",      H2_nucleus1));
    system_H2_small->addAtom(new Hydrogen("3-21G",      H2_nucleus2));
    system_H2_big->  addAtom(new Hydrogen("6-311++G**", H2_nucleus1));
    system_H2_big->  addAtom(new Hydrogen("6-311++G**", H2_nucleus2));
    UnrestrictedHartreeFock un_solver_H2_small(system_H2_small);
    RestrictedHartreeFock   re_solver_H2_small(system_H2_small);
    UnrestrictedHartreeFock un_solver_H2_big  (system_H2_big);
    RestrictedHartreeFock   re_solver_H2_big  (system_H2_big);

    double un_H2_small  = un_solver_H2_small.solveSilently(1e-14,1e4);
    double re_H2_small  = re_solver_H2_small.solveSilently(1e-14,1e4);
    double un_H2_big    = un_solver_H2_big.  solveSilently(1e-14,1e4);
    double re_H2_big    = re_solver_H2_big.  solveSilently(1e-14,1e4);
    double HF_lim_H2    = -1.134;
    printf("%10.6g & %10.6g & %10.6g & %10.6g & %10.6g  \\\\ \n",
           re_H2_small,
           un_H2_small,
           re_H2_big,
           un_H2_big,
           HF_lim_H2); fflush(stdout);

    // CO
    // ========================================================================
    vec CO_nucleus1  {0, 0, 0};
    vec CO_nucleus2  {0, 0, 2.132};
    System* system_CO_small = new System(2);
    System* system_CO_big   = new System(2);
    system_CO_small->addAtom(new Carbon("3-21G",      CO_nucleus1));
    system_CO_small->addAtom(new Oxygen("3-21G",      CO_nucleus2));
    system_CO_big->  addAtom(new Carbon("6-311++G**", CO_nucleus1));
    system_CO_big->  addAtom(new Oxygen("6-311++G**", CO_nucleus2));
    UnrestrictedHartreeFock un_solver_CO_small(system_CO_small);
    RestrictedHartreeFock   re_solver_CO_small(system_CO_small);
    UnrestrictedHartreeFock un_solver_CO_big  (system_CO_big);
    RestrictedHartreeFock   re_solver_CO_big  (system_CO_big);

    double un_CO_small  = un_solver_CO_small.solveSilently(1e-10,1e4);
    double re_CO_small  = re_solver_CO_small.solveSilently(1e-10,1e4);
    double un_CO_big    = un_solver_CO_big.  solveSilently(1e-10,1e4);
    double re_CO_big    = re_solver_CO_big.  solveSilently(1e-10,1e4);
    double HF_lim_CO    = -112.791;
    printf("%10.6g & %10.6g & %10.6g & %10.6g & %10.6g  \\\\ \n",
           re_CO_small,
           un_CO_small,
           re_CO_big,
           un_CO_big,
           HF_lim_CO); fflush(stdout);



    // N2
    // ========================================================================
    vec N2_nucleus1  {0, 0, 0};
    vec N2_nucleus2  {0, 0, 2.074};
    System* system_N2_small = new System(2);
    System* system_N2_big   = new System(2);
    system_N2_small->addAtom(new Nitrogen("3-21G",      N2_nucleus1));
    system_N2_small->addAtom(new Nitrogen("3-21G",      N2_nucleus2));
    system_N2_big->  addAtom(new Nitrogen("6-311++G**", N2_nucleus1));
    system_N2_big->  addAtom(new Nitrogen("6-311++G**", N2_nucleus2));
    UnrestrictedHartreeFock un_solver_N2_small(system_N2_small);
    RestrictedHartreeFock   re_solver_N2_small(system_N2_small);
    UnrestrictedHartreeFock un_solver_N2_big  (system_N2_big);
    RestrictedHartreeFock   re_solver_N2_big  (system_N2_big);

    double un_N2_small  = un_solver_N2_small.solveSilently(1e-10,1e4);
    double re_N2_small  = re_solver_N2_small.solveSilently(1e-10,1e4);
    double un_N2_big    = un_solver_N2_big.  solveSilently(1e-10,1e4);
    double re_N2_big    = re_solver_N2_big.  solveSilently(1e-10,1e4);
    double HF_lim_N2    = -108.997;
    printf("%10.6g & %10.6g & %10.6g & %10.6g & %10.6g  \\\\ \n",
           re_N2_small,
           un_N2_small,
           re_N2_big,
           un_N2_big,
           HF_lim_N2); fflush(stdout);


    // CH4
    // ========================================================================
    vec CH4_nucleus1  {0,            0,             0};
    vec CH4_nucleus2  { 1.67381799,  0.        , -1.18356805};
    vec CH4_nucleus3  {-1.67381799,  0.        , -1.18356805};
    vec CH4_nucleus4  { 0.        ,  1.67381799,  1.18356805};
    vec CH4_nucleus5  { 0.        , -1.67381799,  1.18356805};
    System* system_CH4_small = new System(5);
    System* system_CH4_big   = new System(5);
    system_CH4_small->addAtom(new Carbon  ("3-21G",      CH4_nucleus1));
    system_CH4_small->addAtom(new Hydrogen("3-21G",      CH4_nucleus2));
    system_CH4_small->addAtom(new Hydrogen("3-21G",      CH4_nucleus3));
    system_CH4_small->addAtom(new Hydrogen("3-21G",      CH4_nucleus4));
    system_CH4_small->addAtom(new Hydrogen("3-21G",      CH4_nucleus5));
    system_CH4_big->  addAtom(new Carbon  ("6-311++G**", CH4_nucleus1));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus2));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus3));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus4));
    system_CH4_big->  addAtom(new Hydrogen("6-311++G**", CH4_nucleus5));
    UnrestrictedHartreeFock un_solver_CH4_small(system_CH4_small);
    RestrictedHartreeFock   re_solver_CH4_small(system_CH4_small);
    UnrestrictedHartreeFock un_solver_CH4_big  (system_CH4_big);
    RestrictedHartreeFock   re_solver_CH4_big  (system_CH4_big);

    double un_CH4_small  = un_solver_CH4_small.solveSilently(1e-10,1e4);
    double re_CH4_small  = re_solver_CH4_small.solveSilently(1e-10,1e4);
    double un_CH4_big    = un_solver_CH4_big.  solveSilently(1e-10,1e4);
    double re_CH4_big    = re_solver_CH4_big.  solveSilently(1e-10,1e4);
    double HF_lim_CH4    = -40.225;
    printf("%10.6g & %10.6g & %10.6g & %10.6g & %10.6g  \\\\ \n",
           re_CH4_small,
           un_CH4_small,
           re_CH4_big,
           un_CH4_big,
           HF_lim_CH4); fflush(stdout);

    // NH3
    // ========================================================================
    vec NH3_nucleus1 { 0.0000,    0.12661165,    0.0000};
    vec NH3_nucleus2 { 0.82770004, -0.58959455,  1.55902405};
    vec NH3_nucleus3 { 0.93541443, -0.58959455, -1.49666309};
    vec NH3_nucleus4 {-1.7650042 , -0.58959455, -0.06236096};
    System* system_NH3_small = new System(4);
    System* system_NH3_big   = new System(4);
    system_NH3_small->addAtom(new Nitrogen("3-21G",      NH3_nucleus1));
    system_NH3_small->addAtom(new Hydrogen("3-21G",      NH3_nucleus2));
    system_NH3_small->addAtom(new Hydrogen("3-21G",      NH3_nucleus3));
    system_NH3_small->addAtom(new Hydrogen("3-21G",      NH3_nucleus4));
    system_NH3_big->  addAtom(new Nitrogen("6-311++G**", NH3_nucleus1));
    system_NH3_big->  addAtom(new Hydrogen("6-311++G**", NH3_nucleus2));
    system_NH3_big->  addAtom(new Hydrogen("6-311++G**", NH3_nucleus3));
    system_NH3_big->  addAtom(new Hydrogen("6-311++G**", NH3_nucleus4));
    UnrestrictedHartreeFock un_solver_NH3_small(system_NH3_small);
    RestrictedHartreeFock   re_solver_NH3_small(system_NH3_small);
    UnrestrictedHartreeFock un_solver_NH3_big  (system_NH3_big);
    RestrictedHartreeFock   re_solver_NH3_big  (system_NH3_big);

    double un_NH3_small  = un_solver_NH3_small.solve(1e-10,1e4);
    double re_NH3_small  = re_solver_NH3_small.solve(1e-10,1e4);
    double un_NH3_big    = un_solver_NH3_big.  solve(1e-10,1e4);
    double re_NH3_big    = re_solver_NH3_big.  solve(1e-10,1e4);
    double HF_lim_NH3    = -56.225;
    printf("%10.6g & %10.6g & %10.6g & %10.6g & %10.6g  \\\\ \n",
           re_NH3_small,
           un_NH3_small,
           re_NH3_big,
           un_NH3_big,
           HF_lim_NH3); fflush(stdout);

    // H2O
    // ========================================================================
    vec H2O_nucleus1 { 0.000, 0.000, 0.000};
    vec H2O_nucleus2 {-1.430, 1.108, 0.000};
    vec H2O_nucleus3 { 1.430, 1.108, 0.000};
    System* system_H2O_small = new System(3);
    System* system_H2O_big   = new System(3);
    system_H2O_small->addAtom(new Oxygen  ("3-21G",      H2O_nucleus1));
    system_H2O_small->addAtom(new Hydrogen("3-21G",      H2O_nucleus2));
    system_H2O_small->addAtom(new Hydrogen("3-21G",      H2O_nucleus3));
    system_H2O_big->  addAtom(new Oxygen  ("6-311++G**", H2O_nucleus1));
    system_H2O_big->  addAtom(new Hydrogen("6-311++G**", H2O_nucleus2));
    system_H2O_big->  addAtom(new Hydrogen("6-311++G**", H2O_nucleus3));
    UnrestrictedHartreeFock un_solver_H2O_small(system_H2O_small);
    RestrictedHartreeFock   re_solver_H2O_small(system_H2O_small);
    UnrestrictedHartreeFock un_solver_H2O_big  (system_H2O_big);
    RestrictedHartreeFock   re_solver_H2O_big  (system_H2O_big);

    double un_H2O_small  = un_solver_H2O_small.solveSilently(1e-10,1e4);
    double re_H2O_small  = re_solver_H2O_small.solveSilently(1e-10,1e4);
    double un_H2O_big    = un_solver_H2O_big.  solveSilently(1e-10,1e4);
    double re_H2O_big    = re_solver_H2O_big.  solveSilently(1e-10,1e4);
    double HF_lim_H2O    = -76.065;
    printf("%10.6g & %10.6g & %10.6g & %10.6g & %10.6g  \\\\ \n",
           re_H2O_small,
           un_H2O_small,
           re_H2O_big,
           un_H2O_big,
           HF_lim_H2O); fflush(stdout);

    // HF
    // ========================================================================
    vec HF_nucleus1 { 0, 0, 0};
    vec HF_nucleus2 {1.733,0,0};
    System* system_HF_small = new System(2);
    System* system_HF_big   = new System(2);
    system_HF_small->addAtom(new Fluorine("3-21G",      HF_nucleus1));
    system_HF_small->addAtom(new Hydrogen("3-21G",      HF_nucleus2));
    system_HF_big->  addAtom(new Fluorine("6-311++G**", HF_nucleus1));
    system_HF_big->  addAtom(new Hydrogen("6-311++G**", HF_nucleus2));
    UnrestrictedHartreeFock un_solver_HF_small(system_HF_small);
    RestrictedHartreeFock   re_solver_HF_small(system_HF_small);
    UnrestrictedHartreeFock un_solver_HF_big  (system_HF_big);
    RestrictedHartreeFock   re_solver_HF_big  (system_HF_big);

    double un_HF_small  = un_solver_HF_small.solveSilently(1e-10,1e4);
    double re_HF_small  = re_solver_HF_small.solveSilently(1e-10,1e4);
    double un_HF_big    = un_solver_HF_big.  solveSilently(1e-10,1e4);
    double re_HF_big    = re_solver_HF_big.  solveSilently(1e-10,1e4);
    double HF_lim_HF    = -100.071;
    printf("%10.6g & %10.6g & %10.6g & %10.6g & %10.6g  \\\\ \n",
           re_HF_small,
           un_HF_small,
           re_HF_big,
           un_HF_big,
           HF_lim_HF); fflush(stdout);

    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}




void Examples::ValidationTableDissociation() {
    boost::timer t;

    // LiF
    // ========================================================================
    vec LiF_nucleus1 {0, 0, 0};
    vec LiF_nucleus2 {2.955,0,0};
    System* system   = new System(2);
    system->  addAtom(new Fluorine("6-311++G**", LiF_nucleus1));
    system->  addAtom(new Lithium ("6-311++G**", LiF_nucleus2));
    RestrictedHartreeFock* solver = new RestrictedHartreeFock(system);
    double un_LiF = solver->solve(1e-10,1e4);


    // Li
    vec Li_nucleus {0, 0, 0};
    system = new System(1);
    system->addAtom(new Lithium ("6-311++G**", Li_nucleus));
    UnrestrictedHartreeFock* unsolver = new UnrestrictedHartreeFock(system);
    double un_Li = unsolver->solve(1e-10,1e4);

    // F
    vec F_nucleus {0, 0, 0};
    system = new System(1);
    system->addAtom(new Fluorine("6-311++G**", F_nucleus));
    unsolver = new UnrestrictedHartreeFock(system);
    double un_F = unsolver->solve(1e-10,1e4);

    printf("Li: %10.6g \nF: %10.6g \nLi+F: %10.6g \nLiF: %10.6g \nLiF-(Li+F): %10.6g\n",
           un_Li, un_F, un_Li+un_F,un_LiF, un_LiF-(un_Li+un_F));
    fflush(stdout);

    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}

void Examples::ValidationH2plus() {
    boost::timer t;

    int     n  = 1;
    double L0  = 1.;
    double L1  = 7.;
    double dx  = (L1-L0)/n;

    for (int i=0; i<n; i++) {
        double x = L0 + dx*i;
        vec nucleus1  {0, 0, 0};
        vec nucleus2  {x, 0, 0};

        System* system = new System(2);
        //Hydrogen* Hm = new Hydrogen("cc-pVTZ", nucleus1);
        Hydrogen* Hm = new Hydrogen("6-311++G(2d,2p)", nucleus1);
        //Hydrogen* Hm = new Hydrogen("6-311++G**", nucleus1);
        //Hydrogen* Hm = new Hydrogen("3-21G", nucleus1);
        //Hydrogen* Hm = new Hydrogen("6-31G", nucleus1);
        //Hydrogen* Hm = new Hydrogen("3-21++G", nucleus1);
        //Hydrogen* Hm = new Hydrogen("6-31G**", nucleus1);
        //Hm->setNumberOfElectrons(0);
        system->addAtom(Hm);
        //system->addAtom(new Hydrogen("cc-pVTZ", nucleus2));
        //system->addAtom(new Hydrogen("6-311++G(2d,2p)", nucleus2));
        //system->addAtom(new Hydrogen("6-311++G**", nucleus2));
        //system->addAtom(new Hydrogen("3-21G", nucleus2));
        //system->addAtom(new Hydrogen("6-31G", nucleus2));
        //system->addAtom(new Hydrogen("6-31G**", nucleus2));
        UnrestrictedHartreeFock solver(system);
        double E = solver.solve(1e-14, 1e4);
        printf("%20.15g %20.15g; \n", x, E);
        if (i % 10 == 0) {
            fflush(stdout);
        }
    }

    double elapsedTime = t.elapsed();
    cout << "Elapsed time: " << elapsedTime << endl;
}



void Examples::SingleAtom(int Z,
                          std::string basis,
                          int electrons,
                          int maxIterations,
                          double tollerance,
                          std::string basisFileName) {
    boost::timer t;
    vec nucleus {0, 0, 0};
    System* system   = new System(1);
    Atom* atom;
    switch (Z) {
        case 1:
            atom = new Hydrogen(basis, nucleus);
            break;
        case 2:
            atom = new Helium(basis, nucleus);
            break;
        case 3:
            atom = new Lithium(basis, nucleus);
            break;
        case 4:
            atom = new Beryllium(basis, nucleus);
            break;
        case 5:
            atom = new Boron(basis, nucleus);
            break;
        case 6:
            atom = new Carbon(basis, nucleus);
            break;
        case 7:
            atom = new Nitrogen(basis, nucleus);
            break;
        case 8:
            atom = new Oxygen(basis, nucleus);
            break;
        case 9:
            atom = new Fluorine(basis, nucleus);
            break;
        case 10:
            atom = new Neon(basis, nucleus);
            break;
        default:
            std::cout << "Z = " << Z << " atom not yet implemented." << std::endl;
    }

    if (electrons==-1) electrons = Z;
    atom->setNumberOfElectrons(electrons);
    system->addAtom(atom);
    UnrestrictedHartreeFock* solver = new UnrestrictedHartreeFock(system);
    double E = solver->solve(tollerance,maxIterations);
    double elapsedTime = t.elapsed();
    if (! (basisFileName=="")) {
        solver->dumpBasisToFile(basisFileName);
    }
    cout << "Elapsed time: " << elapsedTime << endl;
}






















