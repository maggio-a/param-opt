VCGPATH = $(HOME)/include/vcglib

TARGET = pinfo

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle
#CONFIG -= qt

QT = core gui

TEMPLATE = app

INCLUDEPATH += $$VCGPATH

SOURCES += main.cpp
SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp

