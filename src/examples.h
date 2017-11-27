#pragma once
#include <string>

class Examples {
public:
    static void firstExample();
    static void secondExample();

    static void Hm();
    static void He();
    static void HeHp();
    static void H2();
    static void H20();
    static void ValidationTableEnergy();
    static void ValidationTableDissociation();
    static void ValidationH2plus();
    static void SingleAtom(int Z,
                           std::string basis="6-311++G**",
                           int electrons=-1,
                           int maxIterations=(int)1e4,
                           double tollerance=1e-8,
                           std::string basisFileName="");
};

