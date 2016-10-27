TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -larmadillo -llapack -lblas

SOURCES += main.cpp \
    Factorizations/hermitegaussian.cpp \
    Factorizations/hermitegaussianintegral.cpp \
    Integrators/electronnucleusintegrator.cpp \
    Integrators/contractedintegrator.cpp \
    Integrators/kineticintegrator.cpp \
    Integrators/overlapintegrator.cpp \
    Math/boysfunction.cpp \
    Orbitals/contractedgaussian.cpp \
    Orbitals/gaussianprimitive.cpp \
    Integrators/electronelectronintegrator.cpp \
    Atoms/atom.cpp \
    Atoms/hydrogen.cpp \
    system.cpp

HEADERS += \
    Factorizations/hermitegaussian.h \
    Factorizations/hermitegaussianintegral.h \
    Integrators/electronnucleusintegrator.h \
    Integrators/contractedintegrator.h \
    Integrators/kineticintegrator.h \
    Integrators/overlapintegrator.h \
    Math/boysfunction.h \
    Orbitals/contractedgaussian.h \
    Orbitals/gaussianprimitive.h \
    Integrators/electronelectronintegrator.h \
    Atoms/atom.h \
    Atoms/hydrogen.h \
    Math/factorial.h \
    system.h
