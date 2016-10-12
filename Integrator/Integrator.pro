TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3


SOURCES += main.cpp \
    ran1.cpp \
    Orbitals/orbital.cpp \
    Orbitals/harmonicoscillator2d.cpp \
    montecarlointegrator.cpp

HEADERS += \
    ran1.h \
    Orbitals/orbital.h \
    Orbitals/harmonicoscillator2d.h \
    montecarlointegrator.h
