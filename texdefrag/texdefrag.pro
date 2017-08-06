VCGPATH = $(HOME)/include/vcglib

TARGET = texdefrag

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle
#CONFIG -= qt

QT = core gui

TEMPLATE = app

INCLUDEPATH += $$VCGPATH $$VCGPATH/eigenlib

unix {
    CONFIG += link_pkgconfig
    PKGCONFIG += glfw3
}

LIBS += -lGL -lGLEW

SOURCES += \
    texdefrag.cpp
SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp

HEADERS += \
    texture_rendering.h \
    mesh.h

