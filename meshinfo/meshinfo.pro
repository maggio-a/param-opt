
include(../common.pri)

SOURCES += meshinfo.cpp \
    ../src/mesh_graph.cpp \
    ../src/uv.cpp \
    ../src/gl_utils.cpp \
    ../src/mesh_utils.cpp \
    ../src/texture.cpp \
    ../src/texture_rendering.cpp \
    ../src/mesh.cpp

SOURCES += \
    $${VCGPATH}/wrap/ply/plylib.cpp \
    $${VCGPATH}/wrap/qt/Outline2ToQImage.cpp \
    $${VCGPATH}/wrap/qt/outline2_rasterizer.cpp

HEADERS += \
    ../src/mesh.h \
    ../src/linmath.h \
    ../src/timer.h \
    ../src/uv.h \
    ../src/mesh_graph.h \
    ../src/vertex_position.h \
    ../src/texture_rendering.h \
    ../src/iterative.h \
    ../src/energy.h \
    ../src/math_utils.h \
    ../src/texture_optimization.h \
    ../src/metric.h \
    ../src/texture.h \
    ../src/pushpull.h \
    ../src/gl_utils.h \
    ../src/uniform_solver.h \
    ../src/dcp_solver.h \
    ../src/mesh_utils.h \
    ../src/mesh_attribute.h \
    ../src/polygon2_triangulator.h \
    ../src/parameterization.h \
    ../src/mean_value_param.h \
    ../src/packing_utils.h


