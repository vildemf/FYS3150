TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp lib.cpp

include(deployment.pri)
qtcAddDeployment()

LIBS += -llapack -lblas -larmadillo

#INCLUDEPATH += /usr/local/include # header files LIBS += -L/usr/local/lib # library files
#LIBS += -L/usr/local/lib # library files
