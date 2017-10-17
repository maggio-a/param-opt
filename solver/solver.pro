VCGPATH = $(HOME)/include/vcglib
exists( $(HOME)/devel/vcglib ) {
     VCGPATH = $(HOME)/devel/vcglib
}

TARGET = solver

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle

QT += core svg

TEMPLATE = app

INCLUDEPATH += ../src $$VCGPATH $$VCGPATH/eigenlib

SOURCES += solver.cpp

SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp

HEADERS += \
    ../src/dcpsolver.h \
    ../src/mesh.h

