TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -llapack -lblas -larmadillo

SOURCES += main.cpp \
    restrictedhartreefock.cpp \
    examples.cpp \
    ../Integrator/integraltable.cpp \
    ../Integrator/Orbitals/harmonicoscillator2d.cpp \
    ../Integrator/Orbitals/hydrogen3d.cpp \
    ../Integrator/Orbitals/orbital.cpp \
    ../Integrator/montecarlointegrator.cpp \
    ../Integrator/ran1.cpp

HEADERS += \
    restrictedhartreefock.h \
    examples.h \
    ../Integrator/integraltable.h \
    ../Integrator/Orbitals/harmonicoscillator2d.h \
    ../Integrator/Orbitals/hydrogen3d.h \
    ../Integrator/Orbitals/orbital.h \
    ../Integrator/montecarlointegrator.h \
    ../Integrator/ran1.h

