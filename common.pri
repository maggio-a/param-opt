VCGPATH = $$PWD/vcglib

CONFIG += console
CONFIG += c++11
CONFIG -= app_bundle



#### QT STUFF ##################################################################

TEMPLATE = app
QT = core gui svg



##### INCLUDE PATH #############################################################

INCLUDEPATH += $$PWD/imgui $$PWD/src $$PWD/triangle $$VCGPATH $$VCGPATH/eigenlib



#### LIBS ######################################################################

unix {
  CONFIG += link_pkgconfig
  PKGCONFIG += glfw3 glew

  LIBS += -lGL -lGLEW

  !exists($$PWD/triangle/triangle.o) {
    error("Object file for triangle not found")
  }
  LIBS += $$PWD/triangle/triangle.o
}

win32 {
  WIN_GLFW_PATH = "." # set to glfw dir
  WIN_GLEW_PATH = "." # set to glew dir

  DEFINES += GLFW_DLL
  INCLUDEPATH += $$WIN_GLFW_PATH/include
  LIBS += -L$$WIN_GLFW_PATH/lib -lglfw3

  INCLUDEPATH += $$WIN_GLEW_PATH/include
  LIBS += -L$$WIN_GLEW_PATH/lib/Release/x64 -lglew32

  !exists($$PWD/triangle/triangle.obj) {
    error("Object file for triangle not found")
  }
  LIBS += $$PWD/triangle/triangle.obj
}



#### PLATFORM SPECIFIC FLAGS ###################################################

unix {
  CMAKE_CXXFLAGS += -g
}

win32 {
  DEFINES += NOMINMAX
}

