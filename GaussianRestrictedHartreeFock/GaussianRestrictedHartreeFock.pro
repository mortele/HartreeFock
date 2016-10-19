TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -larmadillo -llapack -lblas

SOURCES += main.cpp \
    gaussianprimitive.cpp \
    contractedgaussian.cpp

HEADERS += \
    gaussianprimitive.h \
    contractedgaussian.h
