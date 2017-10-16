VCGPATH = $(HOME)/include/vcglib

TARGET = texdefrag

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle

QT = core gui svg

TEMPLATE = app

INCLUDEPATH += ../src $$VCGPATH $$VCGPATH/eigenlib

unix {
    CONFIG += link_pkgconfig
    PKGCONFIG += glfw3
}

LIBS += -lGL -lGLEW

SOURCES += texdefrag.cpp

SOURCES += \
    ../src/mesh.cpp

SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp $${VCGPATH}/wrap/qt/Outline2ToQImage.cpp $${VCGPATH}/wrap/qt/outline2_rasterizer.cpp

HEADERS += \
    ../src/texture_rendering.h \
    ../src/mesh.h \
    ../src/optimizer.h \
    ../src/pushpull.h \
    ../src/linmath.h \
    ../src/gl_util.h \
    ../src/timer.h \
    ../src/uv.h \
    ../src/mesh_graph.h \

