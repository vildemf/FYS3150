TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp lib.cpp

include(deployment.pri)
qtcAddDeployment()

LIBS += -llapack -lblas -larmadillo

HEADERS +=

