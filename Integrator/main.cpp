#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include "Orbitals/orbital.h"
#include "Orbitals/harmonicoscillator2d.h"
#include "Orbitals/hydrogen3d.h"
#include "montecarlointegrator.h"
#include "ran1.h"
#include "integraltable.h"

using std::cout;
using std::endl;

int main() {
    IntegralTable table;
  //table.createTwoBodyTable("../IntegralTables/test.dat", 0, 3); // HarmOsci_2D
    table.createTwoBodyTable("../IntegralTables/test.dat", 1, 3); // Hydrogen_3D

    return 0;
}
