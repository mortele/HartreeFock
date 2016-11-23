TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

INCLUDEPATH += /usr/local/opt/gcc/lib/gcc/6
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -larmadillo -llapack -lblas -lboost_regex

SOURCES += main.cpp \
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
    Atoms/oxygen.cpp

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
    Atoms/oxygen.h
