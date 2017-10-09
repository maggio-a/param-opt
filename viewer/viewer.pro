VCGPATH = $(HOME)/include/vcglib

TARGET = viewer

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle

QT = core gui svg

TEMPLATE = app

INCLUDEPATH += ../src $$VCGPATH $$VCGPATH/eigenlib

#unix {
    CONFIG += link_pkgconfig
    PKGCONFIG += glfw3
#}

LIBS += -lGL -lGLEW

SOURCES += viewer.cpp

SOURCES += \
    ../src/mesh_viewer.cpp \
    ../src/mesh.cpp

SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp $${VCGPATH}/wrap/qt/Outline2ToQImage.cpp $${VCGPATH}/wrap/qt/outline2_rasterizer.cpp

HEADERS += \
    ../src/mesh.h \
    ../src/linmath.h \
    ../src/gl_util.h \
    ../src/timer.h \
    ../src/uv.h \
    ../src/mesh_graph.h \
    ../src/mesh_viewer.h

