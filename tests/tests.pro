include(../defaults.pri)
TEMPLATE = app

CONFIG   += console c++11
CONFIG   -= app_bundle
CONFIG   -= qt

SOURCES += main.cpp

LIBS += -L../src -lHF

SOURCES += integraltester.cpp
HEADERS += integraltester.h
