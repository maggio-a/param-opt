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

INCLUDEPATH += ../imgui ../src $$VCGPATH $$VCGPATH/eigenlib $(HOME)/include/earcut/include

QMAKE_CXXFLAGS += -g

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
    ../../include/vcglib/wrap/gui/trackball.cpp \
    ../src/energy.cpp \
    ../src/texture_optimization.cpp \
    ../src/mesh_graph.cpp \
    ../src/uv.cpp \
    ../src/gl_utils.cpp \
    ../src/mesh_utils.cpp \
    ../src/parameterization.cpp

SOURCES += \
    ../src/mesh_viewer.cpp \
    ../src/mesh.cpp \
    ../src/iterative.cpp

SOURCES += $${VCGPATH}/wrap/ply/plylib.cpp $${VCGPATH}/wrap/qt/Outline2ToQImage.cpp $${VCGPATH}/wrap/qt/outline2_rasterizer.cpp

HEADERS += \
    ../src/mesh.h \
    ../src/linmath.h \
    ../src/timer.h \
    ../src/uv.h \
    ../src/mesh_graph.h \
    ../src/mesh_viewer.h \
    ../src/vertex_position.h \
    ../src/texture_rendering.h \
    ../imgui/imconfig.h \
    ../imgui/imgui_internal.h \
    ../imgui/imgui.h \
    ../imgui/stb_rect_pack.h \
    ../imgui/stb_textedit.h \
    ../imgui/stb_truetype.h \
    ../imgui/imgui_glfw_gl3/imgui_impl_glfw_gl3.h \
    ../src/iterative.h \
    ../src/energy.h \
    ../src/math_utils.h \
    ../src/texture_optimization.h \
    ../src/metric.h \
    ../src/texture.h \
    ../src/pushpull.h \
    ../src/parameterization_checker.h \
    ../src/gl_utils.h \
    ../src/uniform_solver.h \
    ../src/dcp_solver.h \
    ../src/mesh_utils.h \
    ../src/mesh_attribute.h \
    ../src/polygon2_triangulator.h \
    ../src/parameterization.h \
    ../src/mean_value_param.h

DISTFILES += \
    readme.txt

