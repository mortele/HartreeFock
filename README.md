[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
# Hartree-Fock

A complete Hartree-Fock framework employing Gaussian basis functions developed from scratch for [the master thesis of Morten Ledum](https://www.duo.uio.no/handle/10852/61196). The code base has been tested and validated on first and second row atoms and light molecules with s, p, d, and f-type Gaussian orbitals. 

A number of [Pople family basis sets](https://aip.scitation.org/doi/abs/10.1063/1.1677527) for atoms up to and including oxygen is included, along with an assortment of others. The task of adding more basis sets is semi-automated through the use of the python script `basisFileParser.py` which transforms a Tubomole basis file into the `C++` which can be directly input into the appropriate `Atom` subclass.

### Density Functional Theory

A simple DFT implementation is located in the DFT branch. This code has not been properly tested and should be considered unstable currently.

## Building and running

1. Download and install Qt Creator from http://www.qt.io/download-open-source/.
2. Clone the git repository: git@github.com:mortele/HartreeFock.git. 
3. Open Qt Creator and select `Open Project`, open `HF.pro` in the newly created local clone.
4. Select an appropriate compiler for your system when prompted.
5. Press Run to compile and run.

Please note that this has never been tested on Windows machines, but it *might* work. Anyway, if anything breaks you get to keep all the pieces.
