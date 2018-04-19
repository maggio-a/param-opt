VCGPATH = $(HOME)/include/vcglib

TARGET = texdefrag

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle

QT = core gui svg

TEMPLATE = app

QMAKE_CXXFLAGS += -g

INCLUDEPATH += ../src $$VCGPATH $$VCGPATH/eigenlib

unix {
    CONFIG += link_pkgconfig
    PKGCONFIG += glfw3 openblas
}

LIBS += -lGL -lGLEW

DEFINES += "EIGEN_USE_BLAS=1"

SOURCES += texdefrag.cpp \
    ../src/energy.cpp \
    ../src/texture_optimization.cpp \
    ../src/mesh_graph.cpp \
    ../src/uv.cpp \
    ../src/gl_utils.cpp \
    ../src/mesh.cpp \
    ../src/iterative.cpp

SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp $${VCGPATH}/wrap/qt/Outline2ToQImage.cpp $${VCGPATH}/wrap/qt/outline2_rasterizer.cpp

HEADERS += \
    ../src/mesh.h \
    ../src/linmath.h \
    ../src/timer.h \
    ../src/uv.h \
    ../src/mesh_graph.h \
    ../src/vertex_position.h \
    ../src/dcp_solver.h \
    ../src/uniform_solver.h \
    ../src/ext/texcoord_optimization.h \
    ../src/texture_rendering.h \
    ../src/pushpull.h \
    ../src/iterative.h \
    ../src/energy.h \
    ../src/math_utils.h \
    ../src/texture_optimization.h \
    ../src/metric.h \
    ../src/parameterization_checker.h \
    ../src/gl_utils.h \
    ../src/texture.h

