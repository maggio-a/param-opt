VCGPATH = $(HOME)/include/vcglib
exists( $(HOME)/devel/vcglib ) {
     VCGPATH = $(HOME)/devel/vcglib
}

TARGET = viewer

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle

QT = core gui svg

TEMPLATE = app

INCLUDEPATH += ../imgui ../src $$VCGPATH $$VCGPATH/eigenlib

#unix {
    CONFIG += link_pkgconfig
     PKGCONFIG += glfw3 glew
#}

#LIBS += -lGL -lGLEW

SOURCES += viewer.cpp \
    ../imgui/imgui_demo.cpp \
    ../imgui/imgui_draw.cpp \
    ../imgui/imgui.cpp \
    ../imgui/imgui_glfw_gl3/imgui_impl_glfw_gl3.cpp \
    ../../include/vcglib/wrap/gui/trackmode.cpp \
    ../../include/vcglib/wrap/gui/trackball.cpp

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
    ../src/mesh_viewer.h \
    ../src/distortion_pos.h \
    ../src/optimizer.h \
    ../src/vertex_position.h \
    ../src/dcpsolver.h \
    ../src/fixed_border_bijective.h \
    ../src/ext/texcoord_optimization.h \
    ../src/texture_rendering.h \
    ../imgui/imconfig.h \
    ../imgui/imgui_internal.h \
    ../imgui/imgui.h \
    ../imgui/stb_rect_pack.h \
    ../imgui/stb_textedit.h \
    ../imgui/stb_truetype.h \
    ../imgui/imgui_glfw_gl3/imgui_impl_glfw_gl3.h \

DISTFILES += \
    readme.txt

