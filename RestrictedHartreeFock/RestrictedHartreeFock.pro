TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    restrictedhartreefock.cpp \
    examples.cpp

HEADERS += \
    restrictedhartreefock.h \
    examples.h

