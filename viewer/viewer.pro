
include(../common.pri)

SOURCES += viewer.cpp \
    ../src/energy.cpp \
    ../src/texture_optimization.cpp \
    ../src/mesh_graph.cpp \
    ../src/uv.cpp \
    ../src/gl_utils.cpp \
    ../src/mesh_utils.cpp \
    ../src/parameterization.cpp \
    ../src/texture.cpp \
    ../src/mesh_viewer.cpp \
    ../src/mesh.cpp \
    ../src/texture_rendering.cpp \
    ../src/packing_utils.cpp \
    ../src/logging.cpp \
    ../src/iterative_solvers.cpp \
    ../src/parallel.cpp

SOURCES += \
    ../imgui/imgui_demo.cpp \
    ../imgui/imgui_draw.cpp \
    ../imgui/imgui.cpp \
    ../imgui/imgui_glfw_gl3/imgui_impl_glfw_gl3.cpp

SOURCES += \
    $${VCGPATH}/wrap/ply/plylib.cpp \
    $${VCGPATH}/wrap/gui/trackball.cpp \
    $${VCGPATH}/wrap/gui/trackmode.cpp \
    $${VCGPATH}/wrap/qt/outline2_rasterizer.cpp


HEADERS += \
    ../src/mesh.h \
    ../src/linmath.h \
    ../src/timer.h \
    ../src/uv.h \
    ../src/mesh_graph.h \
    ../src/mesh_viewer.h \
    ../src/vertex_position.h \
    ../src/texture_rendering.h \
    ../src/energy.h \
    ../src/math_utils.h \
    ../src/texture_optimization.h \
    ../src/metric.h \
    ../src/texture.h \
    ../src/pushpull.h \
    ../src/gl_utils.h \
    ../src/mesh_utils.h \
    ../src/mesh_attribute.h \
    ../src/polygon2_triangulator.h \
    ../src/parameterization.h \
    ../src/packing_utils.h \
    ../src/logging.h \
    ../src/parallel.h \
    ../src/linear_solvers.h \
    ../src/iterative_solvers.h

HEADERS += \
    ../imgui/imconfig.h \
    ../imgui/imgui_internal.h \
    ../imgui/imgui.h \
    ../imgui/stb_rect_pack.h \
    ../imgui/stb_textedit.h \
    ../imgui/stb_truetype.h \
    ../imgui/imgui_glfw_gl3/imgui_impl_glfw_gl3.h


DISTFILES += \
    readme.txt

