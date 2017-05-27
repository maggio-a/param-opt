VCGPATH = $(HOME)/include/vcglib

TARGET = diag

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app

INCLUDEPATH += $$VCGPATH

SOURCES += main.cpp
SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp
