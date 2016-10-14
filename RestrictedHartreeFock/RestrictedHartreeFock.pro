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
    ../Integrator/integraltable.cpp

HEADERS += \
    restrictedhartreefock.h \
    examples.h \
    ../Integrator/integraltable.h

