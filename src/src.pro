include(../defaults.pri)

CONFIG   += console c++11
CONFIG   -= app_bundle
CONFIG   -= qt

TEMPLATE = lib

TARGET = HF

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3


SOURCES +=  \
    Factorizations/hermitegaussian.cpp \
    Factorizations/hermitegaussianintegral.cpp \
    Integrators/electronnucleusintegrator.cpp \
    Integrators/contractedintegrator.cpp \
    Integrators/kineticintegrator.cpp \
    Integrators/overlapintegrator.cpp \
    Integrators/electronelectronintegrator.cpp \
    Math/boysfunction.cpp \
    Orbitals/contractedgaussian.cpp \
    Orbitals/gaussianprimitive.cpp \
    Atoms/atom.cpp \
    Atoms/hydrogen.cpp \
    Solvers/restrictedhartreefock.cpp \
    Parsers/basissetparser.cpp \
    Parsers/filenameparser.cpp \
    system.cpp \
    Solvers/unrestrictedhartreefock.cpp \
    Solvers/hartreefock.cpp \
    Atoms/oxygen.cpp \
    Atoms/helium.cpp \
    examples.cpp \
    Atoms/lithium.cpp \
    Atoms/beryllium.cpp \
    Atoms/boron.cpp \
    Atoms/carbon.cpp \
    Atoms/nitrogen.cpp \
    Atoms/fluorine.cpp \
    Atoms/neon.cpp \
    Solvers/restricteddft.cpp \
    Integrators/numericalintegrator.cpp \
    Integrators/grid.cpp

HEADERS += \
    Factorizations/hermitegaussian.h \
    Factorizations/hermitegaussianintegral.h \
    Integrators/electronnucleusintegrator.h \
    Integrators/contractedintegrator.h \
    Integrators/kineticintegrator.h \
    Integrators/overlapintegrator.h \
    Integrators/electronelectronintegrator.h \
    Math/boysfunction.h \
    Math/factorial.h \
    Orbitals/contractedgaussian.h \
    Orbitals/gaussianprimitive.h \
    Atoms/atom.h \
    Atoms/hydrogen.h \
    Solvers/restrictedhartreefock.h \
    Parsers/basissetparser.h \
    Parsers/filenameparser.h \
    system.h \
    Solvers/unrestrictedhartreefock.h \
    Solvers/hartreefock.h \
    Atoms/oxygen.h \
    Atoms/helium.h \
    examples.h \
    Atoms/lithium.h \
    Atoms/beryllium.h \
    Atoms/boron.h \
    Atoms/carbon.h \
    Atoms/nitrogen.h \
    Atoms/fluorine.h \
    Atoms/neon.h \
    Solvers/restricteddft.h \
    Integrators/numericalintegrator.h \
    Integrators/grid.h
