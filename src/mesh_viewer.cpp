#include <iostream>
#include <cstdio>

#include <QImage>

#include <vcg/space/ray3.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/intersection2.h>



#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <wrap/io_trimesh/export.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <wrap/gui/trackball.h>

#include <imgui.h>
#include <imgui_internal.h> // to disable elements
#include <imgui_glfw_gl3/imgui_impl_glfw_gl3.h>

#include "mesh.h"
#include "mesh_viewer.h"
#include "gl_utils.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "parameterization_checker.h"
#include "timer.h"
#include "mesh_attribute.h"
#include "mesh_utils.h"
#include "parameterization.h"

#include "linmath.h"

static vcg::Trackball trackball;

const char *vs_text_3D[] = {
    "#version 410 core                                                             \n"
    "                                                                              \n"
    "uniform mat4 modelViewMatrix;                                                 \n"
    "uniform mat4 projectionMatrix;                                                \n"
    "                                                                              \n"
    "in vec3 position;                                                             \n"
    "in vec2 texcoord;                                                             \n"
    "in vec4 distcolor;                                                            \n"
    "out vec2 uv;                                                                  \n"
    "out vec4 dcol;                                                                \n"
    "                                                                              \n"
    "void main(void)                                                               \n"
    "{                                                                             \n"
    "    uv = texcoord;                                                            \n"
    "    dcol = distcolor;                                                         \n"
    "    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0f);  \n"
    "}                                                                             \n"
};

const char *vs_text_texture[] = {
    "#version 410 core                                                             \n"
    "                                                                              \n"
    "uniform mat4 projectionMatrix;                                                \n"
    "uniform vec4 primitiveColor;                                                  \n"
    "                                                                              \n"
    "in vec2 position;                                                             \n"
    "in vec2 texcoord;                                                             \n"
    "in vec4 distcolor;                                                            \n"
    "                                                                              \n"
    "out vec2 uv;                                                                  \n"
    "out vec4 dcol;                                                                \n"
    "                                                                              \n"
    "void main(void)                                                               \n"
    "{                                                                             \n"
    "    uv = texcoord;                                                            \n"
    "    dcol = distcolor;                                                         \n"
    "    gl_Position = projectionMatrix * vec4(position, 0.5f, 1.0f);              \n"
    "}                                                                             \n"
};

const char *fs_text_texture[] = {
    "#version 410 core                                              \n"
    "                                                               \n"
    "const int COLOR_SRC_TEXTURE    = 1;                            \n"
    "const int COLOR_SRC_PRIMITIVE  = 2;                            \n"
    "const int COLOR_SRC_CHECKBOARD = 4;                            \n"
    "                                                               \n"
    "uniform sampler2D tex0;                                        \n"
    "uniform float resolution = 64.0f;                              \n"
    "uniform int colorMask = COLOR_SRC_TEXTURE;                     \n"
    "uniform vec4 weight = vec4(1.0f, 1.0f, 1.0f, 1.0f);            \n"
    "                                                               \n"
    "in vec2 uv;                                                    \n"
    "in vec4 dcol;                                                  \n"
    "out vec4 color;                                                \n"
    "                                                               \n"
    "void main(void)                                                \n"
    "{                                                              \n"
    "    color = weight;                                            \n"
    "    if ((colorMask & COLOR_SRC_TEXTURE) != 0) {                \n"
    "        vec2 tc = uv;                                          \n"
    "        if (tc.s < 0.0) tc = vec2(0.0, 0.0);                   \n"
    "        if (tc.s == 0.0 && tc.t == 0.0) color *= vec4(0.72, 0.2, 0.2, 1); \n"
    "        else color *= texture(tex0, tc);                                 \n"
    "    }                                                          \n"
    "    if ((colorMask & COLOR_SRC_PRIMITIVE) != 0) {              \n"
    "        color *= dcol;                                         \n"
    "    }                                                          \n"
    "    if ((colorMask & COLOR_SRC_CHECKBOARD) != 0) {             \n"
    "        float u_mod = mod(floor(uv.s * resolution), 2.0f);     \n"
    "        float v_mod = mod(floor(uv.t * resolution), 2.0f);     \n"
    "        float c = abs(u_mod - v_mod);                          \n"
    "        c = clamp(c+0.4f, 0.0f, 1.0f);                         \n"
    "        color *= vec4(c, c, c, 1.0f);                          \n"
    "    }                                                          \n"
    "}                                                              \n"
};

void TrackballGlfwMouseButtonCapture(GLFWwindow *window, int x, int y, int button, int action)
{
    int trackButton = vcg::Trackball::BUTTON_NONE;
    switch (button) {
    case GLFW_MOUSE_BUTTON_LEFT: trackButton = vcg::Trackball::BUTTON_LEFT; break;
    case GLFW_MOUSE_BUTTON_MIDDLE: trackButton = vcg::Trackball::BUTTON_MIDDLE; break;
    case GLFW_MOUSE_BUTTON_RIGHT:/* trackButton = vcg::Trackball::BUTTON_RIGHT;*/ break;
    default: assert(0 && "Mouse button input");
    }

    if (glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS) {
        trackButton |= vcg::Trackball::KEY_ALT;
    }

    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS) {
        trackButton |= vcg::Trackball::KEY_CTRL;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS) {
        trackButton |= vcg::Trackball::KEY_SHIFT;
    }

    if (trackButton == vcg::Trackball::BUTTON_NONE) return;

    if (action == GLFW_PRESS) trackball.MouseDown(x, y, trackButton);
    else if (action == GLFW_RELEASE) trackball.MouseUp(x, y, trackButton);
}

void TrackballGlfwKeyCapture(int key, int action)
{
    vcg::Trackball::Button trackKey = vcg::Trackball::BUTTON_NONE;
    switch (key) {
    case GLFW_KEY_LEFT_ALT:
    case GLFW_KEY_RIGHT_ALT:
        trackKey = vcg::Trackball::KEY_ALT; break;
    case GLFW_KEY_LEFT_CONTROL:
    case GLFW_KEY_RIGHT_CONTROL:
        trackKey = vcg::Trackball::KEY_CTRL; break;
    case GLFW_KEY_LEFT_SHIFT:
    case GLFW_KEY_RIGHT_SHIFT:
        trackKey = vcg::Trackball::KEY_SHIFT; break;
    default: return;
    }

    if (action == GLFW_PRESS) trackball.ButtonDown(trackKey);
    else if (action == GLFW_RELEASE) trackball.ButtonUp(trackKey);
}

// GLFW callbacks

void MeshViewer::MouseButtonCallback(GLFWwindow *window, int button, int action, int /*mods*/)
{
    static Timer clickTimer;
    static bool waitingDoubleClick = false;

    bool doubleClick = false;

    if (action == GLFW_PRESS && button >= 0 && button < 3) {
        g_MouseJustPressed[button] = true;

        float timeSinceLastClick = clickTimer.TimeElapsed();
        if (waitingDoubleClick && timeSinceLastClick < 0.6f) { // if after more than 0.6 sec, listen for a new double click
            if (clickTimer.TimeElapsed() < 0.25f) {
                doubleClick = true;
            }
            waitingDoubleClick = false;
        } else {
            clickTimer.Reset();
            waitingDoubleClick = true;
        }
    }


    ImGuiIO& io = ImGui::GetIO();

    // forward input to app logic only if gui is not involved or it is a release event
    if (io.WantCaptureMouse == false || action == GLFW_RELEASE) {
        MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);

        if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
            viewer->PickRegion();
        }

        if (viewer->InPerspectiveView()) {
            if (doubleClick) {
                viewer->CenterPerspectiveViewFromMouse();
            } else {
                TrackballGlfwMouseButtonCapture(window, viewer->_xpos, viewer->info.height - viewer->_ypos, button, action);
                if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                    viewer->_dragMode = DragMode::PERSPECTIVE;
                }
            }
        }
        else if (viewer->InTextureView()) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                viewer->_dragMode = DragMode::TEXTURE;
            }
        }
        else if (viewer->InDetailView()) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                viewer->_dragMode = DragMode::DETAIL;
            }
        }

        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
            TrackballGlfwMouseButtonCapture(window, viewer->_xpos, viewer->info.height - viewer->_ypos, button, action);
            viewer->_dragMode = DragMode::DISABLED;
        }
    }
}

void MeshViewer::CursorPositionCallback(GLFWwindow *window, double xpos, double ypos)
{
    MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);

    int windowWidth, windowHeight, fbWidth, fbHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);

    double scale = fbWidth / (double) windowWidth;

    xpos *= scale, ypos *= scale;

    trackball.MouseMove((int) xpos, viewer->info.height - ((int) ypos));

    if (viewer->_dragMode != DragMode::DISABLED) {
        viewer->_dragX += (int) (xpos - viewer->_xpos);
        viewer->_dragY += (int) (ypos - viewer->_ypos);
    }

    viewer->_xpos = xpos;
    viewer->_ypos = ypos;
}

void MeshViewer::ScrollCallback(GLFWwindow* window, double /*xoffset*/, double yoffset)
{
    g_MouseWheel += (float)yoffset;

    ImGuiIO& io = ImGui::GetIO();

    if (io.WantCaptureMouse == false) {
        MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);
        if (viewer->InPerspectiveView()) {
            trackball.MouseWheel(yoffset);
        }
        else if (viewer->InTextureView()) {
            float mx = ((viewer->_xpos-viewer->info.xSplit) / (viewer->info.height/2.0f)) - 0.5f;
            float my = 0.5f - (viewer->_ypos / (viewer->info.height/2.0f));
            yoffset > 0.0f ? viewer->_textureCamera.ZoomIn(mx, my) : viewer->_textureCamera.ZoomOut(mx, my);
        }
        else if (viewer->InDetailView()) {
            float mx = ((viewer->_xpos-viewer->info.xSplit) / (viewer->info.height/2.0f)) - 0.5f;
            float my = 0.5f - ((viewer->_ypos - viewer->info.height/2.0f) / (viewer->info.height/2.0f));
            yoffset > 0.0f ? viewer->_detailCamera.ZoomIn(mx, my) : viewer->_detailCamera.ZoomOut(mx, my);
        }
    }
}

void MeshViewer::KeyCallback(GLFWwindow* window, int key, int, int action, int /*mods*/)
{
    ImGuiIO& io = ImGui::GetIO();
    if (action == GLFW_PRESS)
        io.KeysDown[key] = true;
    if (action == GLFW_RELEASE)
        io.KeysDown[key] = false;

    io.KeyCtrl = io.KeysDown[GLFW_KEY_LEFT_CONTROL] || io.KeysDown[GLFW_KEY_RIGHT_CONTROL];
    io.KeyShift = io.KeysDown[GLFW_KEY_LEFT_SHIFT] || io.KeysDown[GLFW_KEY_RIGHT_SHIFT];
    io.KeyAlt = io.KeysDown[GLFW_KEY_LEFT_ALT] || io.KeysDown[GLFW_KEY_RIGHT_ALT];
    io.KeySuper = io.KeysDown[GLFW_KEY_LEFT_SUPER] || io.KeysDown[GLFW_KEY_RIGHT_SUPER];

    if (io.WantCaptureKeyboard == false) {
        MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);
        if (viewer->InPerspectiveView()) {
            TrackballGlfwKeyCapture(key, action);
        }
    }
}

void MeshViewer::FramebufferSizeCallback(GLFWwindow *window, int width, int height)
{
    MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);
    if (width <= 0) width = 1;
    if (height <= 0) height = 1;
    int borderToSplitWidth = width > (height/2) ? (width - (height/2)) : 1;
    viewer->info.width = width;
    viewer->info.height = height;
    viewer->info.xSplit = borderToSplitWidth;

    // for the perspective viewport, discard 10% of the borderToSplitWidth to make room for the controls

    viewer->info.perspectiveViewport[0] = std::min(160, (int) (0.2f * borderToSplitWidth));
    viewer->info.perspectiveViewport[1] = 0;
    viewer->info.perspectiveViewport[2] = borderToSplitWidth - viewer->info.perspectiveViewport[0];
    viewer->info.perspectiveViewport[3] = viewer->info.height;

    viewer->info.perspectiveViewAspect = viewer->info.perspectiveViewport[2] / (float) viewer->info.perspectiveViewport[3];

    viewer->info.detailViewport[0] = viewer->info.xSplit;
    viewer->info.detailViewport[1] = 0;
    viewer->info.detailViewport[2] = viewer->info.height/2;
    viewer->info.detailViewport[3] = viewer->info.height/2;

    viewer->info.textureViewport[0] = viewer->info.xSplit;
    viewer->info.textureViewport[1] = viewer->info.height/2;
    viewer->info.textureViewport[2] = viewer->info.height/2;
    viewer->info.textureViewport[3] = viewer->info.height/2;
}


// Member functions

MeshViewer::MeshViewer(GraphHandle gh, const std::string& fileName_)
    : graph{gh},
      _currentTexture{gh->textureObject},
      gm{std::make_shared<GraphManager>(gh, std::unique_ptr<EdgeWeightFunction>(new W3D(gh->mesh)))},
      fileName{fileName_},
      minRegionSize{0},
      _textureCamera{},
      _detailCamera{},
      strategy(DefaultStrategy())
{
    assert(_currentTexture->ArraySize() == 1 && "Only single texture meshes are supported by the viewer");
    Mesh& m = graph->mesh;
    if (m.FN() < 100000)
        minRegionSize = 1000;
    else if (m.FN() < 300000)
        minRegionSize = 5000;
    else
        minRegionSize = 10000;

    //std::cout << "fixme" << std::endl;
    //minRegionSize = 100000;

    std::size_t numRegions = graph->Count();
    regionColors.reserve(numRegions);
    for (const auto& c : graph->charts) {
        auto color = vcg::Color4f::Scatter(20, c.first % 20, 0.75f);
        regionColors.insert(std::make_pair(c.first, color/255.0f));
    }
}

bool MeshViewer::InPerspectiveView()
{
    return _xpos < info.xSplit;
}

bool MeshViewer::InTextureView()
{
    return _xpos > info.xSplit && _ypos < (info.height/2);
}

bool MeshViewer::InDetailView()
{
    return _xpos > info.xSplit && _ypos > (info.height/2);
}

void MeshViewer::TexturePick()
{
    // Perform intersection tests in uv space

    float dx = ((_xpos-info.xSplit) / (info.height/2.0f)) - 0.5f;
    dx *= _textureCamera.viewSize;
    float dy = 0.5f - (_ypos / (info.height/2.0f));
    dy *= _textureCamera.viewSize;

    Point2d p{_textureCamera.x + dx, _textureCamera.y + dy};

    const Mesh& m = graph->mesh;
    const MeshFace *fp = nullptr;
    for (auto &f : m.face) {
        vcg::Triangle2<double> uvFace{f.cWT(0).P(), f.cWT(1).P(), f.cWT(2).P()};
        if (vcg::IsInsideTrianglePoint(uvFace, p)) {
            fp = &f;
            break;
        }
    }

    if (fp != nullptr) {
        auto CCIDh = tri::Allocator<Mesh>::GetPerFaceAttribute<RegionID>(graph->mesh, "ConnectedComponentID");
        RegionID selectionID = CCIDh[fp];
        Select(selectionID);
    } else {
        ClearSelection();
    }
}

void MeshViewer::ClearSelection()
{
    if (selectedRegions.size() > 0) {
        for (auto vao : selectionVao) glDeleteVertexArrays(1, &vao);
        for (auto vbo : _vertexBuffers.selection) glDeleteBuffers(1, &vbo);

        for (auto vao : highlightVao) glDeleteVertexArrays(1, &vao);
        for (auto vbo : _vertexBuffers.highlight) glDeleteBuffers(1, &vbo);

        glDeleteVertexArrays(1, &_detailView.vao);
        glDeleteVertexArrays(1, &_detailView.borderVao);
        glDeleteBuffers(1, &_vertexBuffers.detail);
        _detailView.count = 0;
        glDeleteBuffers(1, &_vertexBuffers.detailBorder);
        _detailView.borderCount = 0;
        _detailView.seamCount = 0;

        for (auto& entry : selectedRegions) glDeleteTextures(1, &entry.second.texicon);

        selectionVao.clear();
        _vertexBuffers.selection.clear();

        highlightVao.clear();
        _vertexBuffers.highlight.clear();

        _detailView.vao = 0;
        _detailView.borderVao = 0;
        _vertexBuffers.detail = 0;
        _vertexBuffers.detailBorder = 0;

        selectedRegions.clear();
        primaryCharts.clear();

        shellGroup = nullptr;
        parameterizer = nullptr;
    }
}

void MeshViewer::Select(const RegionID id)
{
    assert(graph->GetChart(id) != nullptr);
    UpdateSelection(id);
}

void MeshViewer::UpdateDetailBuffers()
{
    glUseProgram(_detailView.program);

    Mesh& shell = parameterizer->Shell();
    Mesh& m = graph->mesh;
    auto ia = GetFaceIndexAttribute(shell);
    auto wtcs = GetWedgeTexCoordStorageAttribute(m);

    // Shell buffer

    glBindVertexArray(_detailView.vao);
    glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.detail);

    GLint buffersz;
    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &buffersz);
    GLint requiredsz = 15 * shell.FN() * sizeof(float);
    if (buffersz < requiredsz) { // resize buffer
        glBufferData(GL_ARRAY_BUFFER, requiredsz, NULL, GL_DYNAMIC_DRAW);
    }

    std::vector<float> borderVertexData;
    std::vector<float> seamVertexData;

    float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto& sf : shell.face) {
        for (int i = 0; i < 3; ++i) {
            if (face::IsBorder(sf, i)) {
                borderVertexData.push_back(sf.P(i).X());
                borderVertexData.push_back(sf.P(i).Y());
                borderVertexData.push_back(sf.P((i+1)%3).X());
                borderVertexData.push_back(sf.P((i+1)%3).Y());
            } else if (sf.IsF(i)) {
                seamVertexData.push_back(sf.P(i).X());
                seamVertexData.push_back(sf.P(i).Y());
                seamVertexData.push_back(sf.P((i+1)%3).X());
                seamVertexData.push_back(sf.P((i+1)%3).Y());
            }
            *buffptr++ = sf.P(i).X();
            *buffptr++ = sf.P(i).Y();
            vcg::Color4b color = vcg::Color4b::Black;
            if (sf.IsMesh()) {
                if (ia[sf] == -1) {
                    std::cout << tri::Index<Mesh>(shell, sf) << std::endl;
                    assert(ia[sf] != -1);
                }
                auto& f = m.face[ia[sf]];
                *buffptr++ = wtcs[f].tc[i].U() / (double) _currentTexture->TextureWidth(0);
                *buffptr++ = wtcs[f].tc[i].V() / (double) _currentTexture->TextureHeight(0);
                /*
                if (uvRatio > 1) {
                    *buffptr++ = wtcs[f].tc[i].U() / uvRatio;
                    *buffptr++ = wtcs[f].tc[i].V();
                } else if (uvRatio <= 1) {
                    *buffptr++ = wtcs[f].tc[i].U();
                    *buffptr++ = wtcs[f].tc[i].V() * uvRatio;
                }
                */
            } else {
                *buffptr++ = 0.0;
                *buffptr++ = 0.0;
            }

            if (shellColorMode == NONE)
                color = vcg::Color4b::White;
            else if (shellColorMode == FACE)
                color = sf.cC();
            else if (shellColorMode == VERTEX) {
                color = sf.V(i)->C();
            }

            unsigned char *colorptr = (unsigned char *) buffptr;
            *colorptr++ = color[0];
            *colorptr++ = color[1];
            *colorptr++ = color[2];
            *colorptr++ = color[3];
            buffptr++;
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    _detailView.count = (GLsizei) (shell.FN() * 3);

    _detailView.attributes.loc_position = glGetAttribLocation(_detailView.program, "position");
    glVertexAttribPointer(_detailView.attributes.loc_position, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), 0);
    glEnableVertexAttribArray(_detailView.attributes.loc_position);

    _detailView.attributes.loc_texcoord = glGetAttribLocation(_detailView.program, "texcoord");
    glVertexAttribPointer(_detailView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (const GLvoid *) (2*sizeof(float)));
    glEnableVertexAttribArray(_detailView.attributes.loc_texcoord);

    _detailView.attributes.loc_color = glGetAttribLocation(_detailView.program, "distcolor");
    glVertexAttribPointer(_detailView.attributes.loc_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, 5 * sizeof(float), (const GLvoid *) (4*sizeof(float)));
    glEnableVertexAttribArray(_detailView.attributes.loc_color);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Border buffer

    glBindVertexArray(_detailView.borderVao);
    glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.detailBorder);

    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &buffersz);
    requiredsz = (borderVertexData.size() + seamVertexData.size()) * sizeof(float);
    if (buffersz < requiredsz) { // resize buffer
        glBufferData(GL_ARRAY_BUFFER, requiredsz, NULL, GL_DYNAMIC_DRAW);
    }

    CheckGLError();
    glBufferSubData(GL_ARRAY_BUFFER, 0, borderVertexData.size() * sizeof(float), &borderVertexData[0]);
    CheckGLError();
    glBufferSubData(GL_ARRAY_BUFFER, borderVertexData.size() * sizeof(float), seamVertexData.size() * sizeof(float), &seamVertexData[0]);
    _detailView.borderCount = (GLsizei) (borderVertexData.size()/2);
    _detailView.seamCount = (GLsizei) (seamVertexData.size()/2);

    CheckGLError();
    glVertexAttribPointer(_detailView.attributes.loc_position, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(_detailView.attributes.loc_position);

    CheckGLError();
    glVertexAttribPointer(_detailView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(_detailView.attributes.loc_texcoord);

    CheckGLError();
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    CheckGLError();
}

void MeshViewer::UpdateSelection(const RegionID id)
{
    // Allocate new buffers if necessary
    std::set<ChartHandle> newCharts;
    std::size_t newElements = 0;

    if (primaryCharts.count(id) == 1) {
        // de-select
        primaryCharts.erase(id);
        selectedRegions[id].referenceCount--;
        for (auto& chart : graph->GetChart(id)->adj) {
            selectedRegions[chart->id].referenceCount--;
        }
        if (primaryCharts.size() == 0) {
            ClearSelection();
            return;
        }
    }
    else {
        // make the selected chart primary
        primaryCharts[id] = 1;

        if (selectedRegions.count(id) == 0) {
            // new buffers will be allocated later
            newCharts.insert(graph->GetChart(id));
            newElements += graph->GetChart(id)->FN();
        }
        else {
            selectedRegions[id].referenceCount++;
        }
        for (auto chart : graph->GetChart(id)->adj) {
            if (selectedRegions.count(chart->id) == 0) {
                newCharts.insert(chart);
                newElements += chart->FN();
            } else {
                selectedRegions[chart->id].referenceCount++;
            }
        }
    }
    if (newElements > 0) {
        std::size_t bufferIndex = _vertexBuffers.selection.size();
        GLint first = 0;
        GLuint selection_vbo;
        glGenBuffers(1, &selection_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, selection_vbo);
        glBufferData(GL_ARRAY_BUFFER, newElements * 15 * sizeof(float), NULL, GL_STATIC_DRAW);
        float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        for (auto chart : newCharts) {
            selectedRegions.insert(std::make_pair(chart->id,
                                                 SelectionBufferInfo{chart, bufferIndex, first, (GLsizei)(chart->FN() * 3), -1, 0, 1}));
            for (auto fptr : chart->fpVec) {
                for (int i = 0; i < 3; ++i) {
                    *buffptr++ = fptr->cV(i)->P().X();
                    *buffptr++ = fptr->cV(i)->P().Y();
                    *buffptr++ = fptr->cV(i)->P().Z();
                    *buffptr++ = fptr->cWT(i).U() / (double) _currentTexture->TextureWidth(0);
                    *buffptr++ = fptr->cWT(i).V() / (double) _currentTexture->TextureHeight(0);
                }
            }
            first += chart->FN() * 3;
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);

        GLuint selection_vao;
        glGenVertexArrays(1, &selection_vao);
        glBindVertexArray(selection_vao);

        _perspectiveView.selection.attributes.loc_position = glGetAttribLocation(_perspectiveView.selection.program, "position");
        glVertexAttribPointer(_perspectiveView.selection.attributes.loc_position, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
        glEnableVertexAttribArray(_perspectiveView.selection.attributes.loc_position);

        _perspectiveView.selection.attributes.loc_texcoord = glGetAttribLocation(_perspectiveView.selection.program, "texcoord");
        glVertexAttribPointer(_perspectiveView.selection.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float),
                              (const GLvoid *) (3*sizeof(float)));
        glEnableVertexAttribArray(_perspectiveView.selection.attributes.loc_texcoord);

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        CheckGLError();

        GLuint highlight_vbo, highlight_vao;
        glGenBuffers(1, &highlight_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, highlight_vbo);
        glBufferData(GL_ARRAY_BUFFER, newCharts.size() * 16 * sizeof(float), NULL, GL_STATIC_DRAW);
        float *bufferBase = buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        for (auto& chart : newCharts) {
            selectedRegions[chart->id].first_highlight = (buffptr - bufferBase) / 2;

            Point2d center = chart->UVBox().Center();
            *buffptr++ = center.X(); *buffptr++ = center.Y() - 0.1f;
            *buffptr++ = center.X(); *buffptr++ = center.Y() + 0.1f;
            *buffptr++ = center.X() - 0.1f; *buffptr++ = center.Y();
            *buffptr++ = center.X() + 0.1f; *buffptr++ = center.Y();

            *buffptr++ = center.X(); *buffptr++ = center.Y() - 0.1f;
            *buffptr++ = center.X(); *buffptr++ = center.Y() + 0.1f;
            *buffptr++ = center.X() - 0.1f; *buffptr++ = center.Y();
            *buffptr++ = center.X() + 0.1f; *buffptr++ = center.Y();
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);

        glGenVertexArrays(1, &highlight_vao);
        glBindVertexArray(highlight_vao);

        _textureView.highlight.attributes.loc_texcoord = glGetAttribLocation(_textureView.program, "texcoord");
        glVertexAttribPointer(_textureView.highlight.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(_textureView.highlight.attributes.loc_texcoord);

        _textureView.highlight.attributes.loc_position = glGetAttribLocation(_textureView.program, "position");
        glVertexAttribPointer(_textureView.highlight.attributes.loc_position, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(_textureView.highlight.attributes.loc_position);

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        CheckGLError();

        selectionVao.push_back(selection_vao);
        highlightVao.push_back(highlight_vao);
        _vertexBuffers.selection.push_back(selection_vbo);
        _vertexBuffers.highlight.push_back(highlight_vbo);

        // generate icons - reuse detail view program
        /*
        glUseProgram(_detailView.program);

        glBindBuffer(GL_ARRAY_BUFFER, selection_vbo);
        GLuint icon_vao;
        glGenVertexArrays(1, &icon_vao);
        glBindVertexArray(icon_vao);

        GLint loc_texcoord = glGetAttribLocation(_detailView.program, "texcoord");
        glVertexAttribPointer(loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (const GLvoid *) (3*sizeof(float)));
        glEnableVertexAttribArray(loc_texcoord);

        GLint loc_position = glGetAttribLocation(_detailView.program, "position");
        glVertexAttribPointer(loc_position, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (const GLvoid *) (3*sizeof(float)));
        glEnableVertexAttribArray(loc_position);

        GLint loc_projection = glGetUniformLocation(_detailView.program, "projectionMatrix");

        GLint loc_colorMask = glGetUniformLocation(_detailView.program, "colorMask");
        GLint loc_weight = glGetUniformLocation(_detailView.program, "weight");
        glUniform1i(loc_colorMask, ColorMask_TEXTURE);
        glUniform4f(loc_weight, 1.0f, 1.0f, 1.0f, 1.0f);

        GLuint fbo;
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        for (auto& chart : newCharts) {
            SelectionBufferInfo& sbi = selectedRegions[chart->id];
            vcg::Box2d box = chart->UVBox();

            mat4x4 projection;
            float halfsz = std::max(box.Dim().X(), box.Dim().Y()) / 2.0f;
            mat4x4_ortho(projection, box.Center().X() - halfsz, box.Center().X() + halfsz, box.Center().Y() - halfsz, box.Center().Y() + halfsz, 1.0f, -1.0f);

            glUniformMatrix4fv(loc_projection, 1, GL_FALSE, (const float *) projection);

            // the icon will be a 64x64 texture
            GLuint texicon = 0;
            glGenTextures(1, &texicon);
            glBindTexture(GL_TEXTURE_2D, texicon);
            glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, 64, 64);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texicon, 0);
            glBindTexture(GL_TEXTURE_2D, 0);

            glActiveTexture(GL_TEXTURE0);
            _currentTexture->Bind();

            glViewport(0, 0, 64, 64);

            glDisable(GL_DEPTH_TEST);
            glDisable(GL_STENCIL_TEST);

            glDrawBuffer(GL_COLOR_ATTACHMENT0);

            if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) assert(0 && "Framebuffer state not supported");

            glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            glDrawArrays(GL_TRIANGLES, sbi.first, sbi.count);

            sbi.texicon = texicon;
        }
        glUseProgram(0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDeleteVertexArrays(1, &icon_vao);
        glDeleteFramebuffers(1, &fbo);
        */

        CheckGLError();
    }

    // parameterize selection and generate rendering data

    /*
    if (_vertexBuffers.detail != 0) {
        glDeleteBuffers(1, &_vertexBuffers.detail);
        _vertexBuffers.detail = 0;
        _detailView.count = 0;
        glDeleteBuffers(1, &_vertexBuffers.detailBorder);
        _detailView.borderCount = 0;
        _detailView.seamCount = 0;
    }
    */

    std::vector<ChartHandle> charts;
    for (auto& entry : primaryCharts) charts.push_back(graph->GetChart(entry.first));
    auto status = gm->CollapseAllowed(charts.begin(), charts.end());
    if (status.first == gm->Collapse_OK) { // the regions can be parameterized together
        // Build the face group that needs to be parameterized
        shellGroup = std::make_shared<FaceGroup>(graph->mesh, INVALID_ID);
        for (auto& c : charts)
            for (auto fptr : c->fpVec)
                shellGroup->AddFace(fptr);
        // Parameterize the aggregate chart, build the vertex buffer and restore the original state
        parameterizer = std::make_shared<ParameterizerObject>(shellGroup, strategy);

        // if necessary, initialized detail vaos and buffers
        if (_detailView.vao == 0) {
            glGenVertexArrays(1, &_detailView.vao);
            glGenBuffers(1, &_vertexBuffers.detail);
            glGenVertexArrays(1, &_detailView.borderVao);
            glGenBuffers(1, &_vertexBuffers.detailBorder);
        }

        UpdateDetailBuffers();

        SetupDetailView(parameterizer->Shell());
        CheckGLError();
    } else {
        /// todo check error code
        std::cout << "Warning: current selection cannot be parameterized: ";
        if (status.first == gm->Collapse_ERR_DISCONNECTED) std::cout << "disconnected selection";
        else if (status.first == gm->Collapse_ERR_UNFEASIBLE) std:: cout << "unfeasible selection";
        else std::cout << "unhandled error code";
        std::cout << std::endl;
    }
}

bool MeshViewer::IntersectionMouseRayModel(Mesh::ConstFacePointer *fp, float &u, float &v)
{
    using vcg::Point3f;
    using vcg::Ray3f;

    // Compute ray from mouse position in perspective view

    mat4x4 model, invModel;
    mat4x4_dup(model, _meshTransform.trackballMatrix);
    //mat4x4_dup(model, _meshTransform.positionMatrix);
    //mat4x4_mul(model, _meshTransform.trackballMatrix, model);
    //mat4x4_mul(model, _meshTransform.scaleMatrix, model);
    mat4x4_invert(invModel, model);

    mat4x4 invProj, invView, ndcToWorld;
    mat4x4_invert(invProj, _meshTransform.projectionMatrix);
    mat4x4_invert(invView, _meshTransform.viewMatrix);
    mat4x4_mul(ndcToWorld, invView, invProj);

    float x = ((_xpos - info.perspectiveViewport[0]) / (float) info.perspectiveViewport[2]) * 2.0f - 1.0f;
    float y = ((info.height - _ypos) / (float) info.height) * 2.0f - 1.0f;

    vec4 ndc_ray = {x, y, -1.0f, 1.0f}, worldRay, modelRay;
    mat4x4_mul_vec4(worldRay, ndcToWorld, ndc_ray);
    worldRay[2] = -1.0f; worldRay[3] = 0.0f;
    mat4x4_mul_vec4(modelRay, invModel, worldRay);
    vec4_norm(worldRay, modelRay);

    vec4 eye = {_perspectiveCamera.eye[0], _perspectiveCamera.eye[1], _perspectiveCamera.eye[2], 1.0f};
    vec4 modelEye;
    mat4x4_mul_vec4(modelEye, invModel, eye);

    Ray3d selectionRay {
        Point3d(modelEye[0], modelEye[1], modelEye[2]),
        Point3d(modelRay[0], modelRay[1], modelRay[2])
    };

    // Perform intersection tests
    const Mesh& m = graph->mesh;
    *fp = nullptr;
    double tMin = std::numeric_limits<double>::max();
    for (auto &f : m.face) {
        double t, uface, vface;
        if (vcg::IntersectionRayTriangle(selectionRay, f.cP(0), f.cP(1), f.cP(2), t, uface, vface) && t < tMin) {
            *fp = &f;
            tMin = t;
            u = uface;
            v = vface;
        }
    }
    return tMin != std::numeric_limits<double>::max();
}

void MeshViewer::PerspectivePick()
{
    Mesh::ConstFacePointer fp;
    float u, v;

    //if (fp != nullptr) {
    if (IntersectionMouseRayModel(&fp, u, v)) {
        auto CCIDh = tri::Allocator<Mesh>::GetPerFaceAttribute<RegionID>(graph->mesh, "ConnectedComponentID");
        RegionID selectionID = CCIDh[fp];
        Select(selectionID);
    } else {
        ClearSelection();
    }
}

void MeshViewer::PickRegion()
{
    if (InPerspectiveView()) PerspectivePick();
    else if (InTextureView()) TexturePick();
}

void MeshViewer::CenterPerspectiveViewFromMouse()
{
    if (InPerspectiveView()) {
        Mesh::ConstFacePointer fp;
        float u;
        float v;
        if (IntersectionMouseRayModel(&fp, u, v)) {
            Point3f centerPoint = Point3f::Construct(fp->cP(0) * (1.0f - u - v) + fp->cP(1) * u + fp->cP(2) * v);
            Matrix44f transform = trackball.Matrix();
            trackball.Translate(-(transform*centerPoint));
        }
    }
}

void MeshViewer::InitBuffers()
{
    const Mesh& m = graph->mesh;

    // Load data, for each vertex position and texture coords
    glGenBuffers(1, &_vertexBuffers.mesh);
    glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.mesh);
    glBufferData(GL_ARRAY_BUFFER, m.FN()*18*sizeof(float), NULL, GL_DYNAMIC_DRAW);
    float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for(auto &f : m.face) {
        for (int i = 0; i < 3; ++i) {
            *buffptr++ = f.cV(i)->P().X();
            *buffptr++ = f.cV(i)->P().Y();
            *buffptr++ = f.cV(i)->P().Z();
            /*if (f.cWT(i).U() < 0.0) {
                *buffptr++ = 0.0;
                *buffptr++ = 0.0;
            } else { */
                *buffptr++ = f.cWT(i).U() / (double) _currentTexture->TextureWidth(0);
                *buffptr++ = f.cWT(i).V() / (double) _currentTexture->TextureHeight(0);
            //}
            unsigned char *colorptr = (unsigned char *) buffptr;
            *colorptr++ = f.cC()[0];
            *colorptr++ = f.cC()[1];
            *colorptr++ = f.cC()[2];
            *colorptr++ = f.cC()[3];
            buffptr++;
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    CheckGLError();

    // Load texture data
    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind(0);

    CheckGLError();

    // Initialize vertex array objects

    // Perspective view
    glGenVertexArrays(1, &_perspectiveView.vao);
    glBindVertexArray(_perspectiveView.vao);

    glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.mesh);

    _perspectiveView.attributes.loc_position = glGetAttribLocation(_perspectiveView.program, "position");
    glVertexAttribPointer(_perspectiveView.attributes.loc_position, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), 0);
    glEnableVertexAttribArray(_perspectiveView.attributes.loc_position);

    _perspectiveView.attributes.loc_texcoord = glGetAttribLocation(_perspectiveView.program, "texcoord");
    glVertexAttribPointer(_perspectiveView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void *) (3*sizeof(float)));
    glEnableVertexAttribArray(_perspectiveView.attributes.loc_texcoord);

    _perspectiveView.attributes.loc_distcolor = glGetAttribLocation(_perspectiveView.program, "distcolor");
    glVertexAttribPointer(_perspectiveView.attributes.loc_distcolor, 4, GL_UNSIGNED_BYTE, GL_TRUE, 6*sizeof(float), (void *) (5*sizeof(float)));
    glEnableVertexAttribArray(_perspectiveView.attributes.loc_distcolor);

    // Texture view
    glGenVertexArrays(1, &_textureView.vao);
    glBindVertexArray(_textureView.vao);

    _textureView.attributes.loc_position = glGetAttribLocation(_textureView.program, "position");
    glVertexAttribPointer(_textureView.attributes.loc_position, 2, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void *) (3*sizeof(float)));
    glEnableVertexAttribArray(_textureView.attributes.loc_position);

    _textureView.attributes.loc_texcoord = glGetAttribLocation(_textureView.program, "texcoord");
    glVertexAttribPointer(_textureView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void *) (3*sizeof(float)));
    glEnableVertexAttribArray(_textureView.attributes.loc_texcoord);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    CheckGLError();
}

void MeshViewer::SetupViews()
{
    Mesh& m = graph->mesh;

    // Setup perspective view

    tri::UpdateBounding<Mesh>::Box(m);

    // Define model transform: the mesh will be translated to the origin of the world space,
    // rotated according to accumulated rotations and then scaled

    /// todo remove all unused matrices
    mat4x4_identity(_meshTransform.positionMatrix);
    //mat4x4_translate(_meshTransform.positionMatrix, -m.bbox.Center().X(), -m.bbox.Center().Y(), -m.bbox.Center().Z());
    trackball.SetIdentity();
    trackball.center = {0.0f, 0.0f, 0.0f};
    trackball.radius = 1.0f;
    trackball.track.sca = 3.0f / m.bbox.Diag();
    trackball.track.tra.Import(-m.bbox.Center());
    mat4x4_identity(_meshTransform.trackballMatrix);
    /*float scale = 1.0f / m.bbox.Diag();
    mat4x4_identity(_meshTransform.scaleMatrix);
    mat4x4_scale_aniso(_meshTransform.scaleMatrix, _meshTransform.scaleMatrix, scale, scale, scale);*/

    // view matrix is constant
    mat4x4_identity(_meshTransform.viewMatrix);
    mat4x4_look_at(_meshTransform.viewMatrix, _perspectiveCamera.eye, _perspectiveCamera.target, _perspectiveCamera.up);

    // view and projection matrices will be updated at each draw
    mat4x4_identity(_meshTransform.projectionMatrix);

    _perspectiveView.uniforms.loc_modelView = glGetUniformLocation(_perspectiveView.program, "modelViewMatrix");
    _perspectiveView.uniforms.loc_projection = glGetUniformLocation(_perspectiveView.program, "projectionMatrix");
    _perspectiveView.uniforms.loc_colorMask = glGetUniformLocation(_perspectiveView.program, "colorMask");
    _perspectiveView.uniforms.loc_weight = glGetUniformLocation(_perspectiveView.program, "weight");

    _perspectiveView.selection.uniforms.loc_modelView = glGetUniformLocation(_perspectiveView.selection.program, "modelViewMatrix");
    _perspectiveView.selection.uniforms.loc_projection = glGetUniformLocation(_perspectiveView.selection.program, "projectionMatrix");
    _perspectiveView.selection.uniforms.loc_colorMask = glGetUniformLocation(_perspectiveView.selection.program, "colorMask");
    _perspectiveView.selection.uniforms.loc_weight = glGetUniformLocation(_perspectiveView.selection.program, "weight");

    // Setup texture view

    _textureCamera.Reset();
    _textureView.uniforms.loc_projection = glGetUniformLocation(_textureView.program, "projectionMatrix");
    _textureView.highlight.uniforms.loc_projection = glGetUniformLocation(_textureView.highlight.program, "projectionMatrix");
    _textureView.highlight.uniforms.loc_primitiveColor = glGetUniformLocation(_textureView.highlight.program, "primitiveColor");
}

//void MeshViewer::SetupDetailView(ChartHandle chart)
void MeshViewer::SetupDetailView(const Mesh& detailMesh)
{
    Box2d bbox;
    for (auto& sv : detailMesh.vert) {
        bbox.Add(sv.T().P());
    }

    //_detailCamera.Reset();
    _detailView.uniforms.loc_projection = glGetUniformLocation(_detailView.program, "projectionMatrix");

    //std::cout << "UVBox = " << bbox.DimX() << " by " << bbox.DimY() << std::endl;

    float scale = 1.0f / std::max(bbox.Dim().X(), bbox.Dim().Y());

    mat4x4 positionMatrix, scaleMatrix, pbm;
    mat4x4_translate(positionMatrix, -bbox.Center().X(), -bbox.Center().Y(), 0.0f);
    mat4x4_identity(scaleMatrix);
    mat4x4_scale_aniso(scaleMatrix, scaleMatrix, scale, scale, 1.0f);
    mat4x4_translate(pbm, 0.5f, 0.5f, 0.0f);

    mat4x4_mul(_detailTransform.modelMatrix, scaleMatrix, positionMatrix);
    mat4x4_mul(_detailTransform.modelMatrix, pbm, _detailTransform.modelMatrix);
}

void MeshViewer::UpdateTransforms()
{
    // Copy current state to the trackball
    Matrix44f proj;
    memcpy(&proj, &_meshTransform.projectionMatrix, 16*sizeof(float));
    proj.transposeInPlace();
    Matrix44f mv;
    memcpy(&mv, &_meshTransform.viewMatrix, 16*sizeof(float));
    mv.transposeInPlace();
    int windowView[] = {0, 0, info.width, info.height};
    trackball.camera.SetView(proj.V(), mv.V(), windowView);

    // set the tracball transform
    Matrix44f trackballMatrix = trackball.Matrix().transpose();
    memcpy(&_meshTransform.trackballMatrix, trackballMatrix.V(), 16*sizeof(float));

    switch (_dragMode) {
    case DragMode::PERSPECTIVE:
        break;
    case DragMode::TEXTURE:
        _textureCamera.MoveX(-_dragX / (info.height/2.0f));
        _textureCamera.MoveY(_dragY / (info.height/2.0f));
        break;
    case DragMode::DETAIL:
        _detailCamera.MoveX(-_dragX / (info.height/2.0f));
        _detailCamera.MoveY(_dragY / (info.height/2.0f));
        break;
    default:
        break;
    }
    _dragX = 0.0f;
    _dragY = 0.0f;
}

void MeshViewer::Draw3DView()
{
    mat4x4 model;
    mat4x4 modelView;

    // Build model matrix: translate, orient and scale
    //mat4x4_dup(model, _meshTransform.positionMatrix);
    //mat4x4_mul(model, _meshTransform.trackballMatrix, model);
    //mat4x4_mul(model, _meshTransform.scaleMatrix, model);



   // mat4x4_mul(modelView, _meshTransform.viewMatrix, model);
    mat4x4_mul(modelView, _meshTransform.viewMatrix, _meshTransform.trackballMatrix);

    mat4x4_perspective(_meshTransform.projectionMatrix, 60.0f * M_PI / 180.0f, info.perspectiveViewAspect,
                       _perspectiveCamera.near, _perspectiveCamera.far);

    int *vp = info.perspectiveViewport;
    glViewport(vp[0], vp[1], vp[2], vp[3]);
    glScissor(vp[0], vp[1], vp[2], vp[3]);

    glDrawBuffer(GL_BACK);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind(0);

    glEnable(GL_DEPTH_TEST);

    if (selectedRegions.size() > 0) {
        //glBindVertexArray(_perspectiveView.selection.vao);
        glUseProgram(_perspectiveView.selection.program);

        glUniformMatrix4fv(_perspectiveView.selection.uniforms.loc_modelView, 1, GL_FALSE, (const GLfloat *) modelView);
        glUniformMatrix4fv(_perspectiveView.selection.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *)_meshTransform.projectionMatrix);

        glEnable(GL_STENCIL_TEST);
        glStencilFuncSeparate(GL_FRONT, GL_ALWAYS, 1, 0xFF);
        glStencilOpSeparate(GL_FRONT, GL_KEEP, GL_KEEP, GL_REPLACE);

        glStencilFuncSeparate(GL_BACK, GL_ALWAYS, 0, 0xFF);
        glStencilOpSeparate(GL_BACK, GL_KEEP, GL_KEEP, GL_REPLACE);

        for (const auto& sel : selectedRegions) {
            RegionID id = sel.first;
            const SelectionBufferInfo& sbi = sel.second;
            if (sbi.referenceCount > 0) {
                glBindVertexArray(selectionVao[sbi.bufferIndex]);
                if (primaryCharts.count(id) == 1) {
                    float one[] = {1.0f, 1.0f, 1.0f, 1.0f};
                    glUniform4fv(_perspectiveView.selection.uniforms.loc_weight, 1, one);
                } else {
                    glUniform4fv(_perspectiveView.selection.uniforms.loc_weight, 1, regionColors[id].V());
                }
                glDrawArrays(GL_TRIANGLES, sbi.first, sbi.count);
            }
        }

        // Now draw transparent layer
        glBindVertexArray(_perspectiveView.vao);
        glUseProgram(_perspectiveView.program);

        glUniformMatrix4fv(_perspectiveView.uniforms.loc_modelView, 1, GL_FALSE, (const GLfloat *) modelView);
        glUniformMatrix4fv(_perspectiveView.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *)_meshTransform.projectionMatrix);

        glUniform1i(_perspectiveView.uniforms.loc_colorMask, _perspectiveView.colorMask);
        glUniform4f(_perspectiveView.uniforms.loc_weight, 0.4f, 0.4f, 0.4f, 0.2f);

        // if stencil buffer == 1 the fragment must be discarded since the selection layer has already been drawn
        glStencilFuncSeparate(GL_FRONT_AND_BACK, GL_EQUAL, 0, 0xFF);
        glStencilOpSeparate(GL_FRONT_AND_BACK, GL_KEEP, GL_KEEP, GL_KEEP); // Never change the stencil buffer

        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glDrawArrays(GL_TRIANGLES, 0, graph->mesh.FN() * 3);
        glBindVertexArray(0);

        glDisable(GL_BLEND);
        glDisable(GL_STENCIL_TEST);
        glEnable(GL_DEPTH_TEST);
    }
    else {
        glBindVertexArray(_perspectiveView.vao);
        glUseProgram(_perspectiveView.program);

        glUniformMatrix4fv(_perspectiveView.uniforms.loc_modelView, 1, GL_FALSE, (const GLfloat *) modelView);
        glUniformMatrix4fv(_perspectiveView.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *)_meshTransform.projectionMatrix);

        glUniform1i(_perspectiveView.uniforms.loc_colorMask, _perspectiveView.colorMask);
        glUniform4f(_perspectiveView.uniforms.loc_weight, 1.0f, 1.0f, 1.0f, 1.0f);

        glDrawArrays(GL_TRIANGLES, 0, graph->mesh.FN() * 3);
        glBindVertexArray(0);
    }

    CheckGLError();
}

void MeshViewer::DrawTextureView()
{
    glBindVertexArray(_textureView.vao);
    glUseProgram(_textureView.program);

    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind(0);

    mat4x4 projection;
    float l = _textureCamera.x - (_textureCamera.viewSize / 2.0f);
    float r = _textureCamera.x + (_textureCamera.viewSize / 2.0f);
    float b = _textureCamera.y - (_textureCamera.viewSize / 2.0f);
    float t = _textureCamera.y + (_textureCamera.viewSize / 2.0f);
    mat4x4_ortho(projection, l, r, b, t, 1.0f, -1.0f);

    glUniformMatrix4fv(_textureView.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *) projection);

    int *vp = info.textureViewport;
    glViewport(vp[0], vp[1], vp[2], vp[3]);
    glScissor(vp[0], vp[1], vp[2], vp[3]);

    glDrawBuffer(GL_BACK);
    glDisable(GL_DEPTH_TEST);

    //glClearColor(0.2f, 0.7f, 0.2f, 1.0f);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 0, graph->mesh.FN() * 3);

    if (selectedRegions.size() > 0) {
        glUseProgram(_textureView.highlight.program);

        glUniformMatrix4fv(_textureView.highlight.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *) projection);

        GLint loc_colorMask = glGetUniformLocation(_textureView.highlight.program, "colorMask");
        glUniform1i(loc_colorMask, ColorMask_PRIMITIVE); // set primitive color as color source

        for (const auto& sel : selectedRegions) {
            RegionID id = sel.first;
            const SelectionBufferInfo& sbi = sel.second;
            if (sbi.referenceCount > 0) {
                glBindVertexArray(highlightVao[sbi.bufferIndex]);
                if (primaryCharts.count(id) == 1) {
                    float white[] = {1.0f, 1.0f, 1.0f, 1.0f};
                    glUniform4fv(_textureView.highlight.uniforms.loc_primitiveColor, 1, white);
                } else {
                    glUniform4fv(_textureView.highlight.uniforms.loc_primitiveColor, 1, regionColors[id].V());
                }
                glDrawArrays(GL_LINES, sbi.first_highlight, 4);
            }
        }
    }

    glBindVertexArray(0);

    CheckGLError();
}

void MeshViewer::DrawDetailView()
{
    SetupDetailView(parameterizer->Shell());
    glBindVertexArray(_detailView.vao);
    glUseProgram(_detailView.program);

    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind(0);

    // mixing model and projection together...
    mat4x4 transform, projection;
    float l = _detailCamera.x - (_detailCamera.viewSize / 2.0f);
    float r = _detailCamera.x + (_detailCamera.viewSize / 2.0f);
    float b = _detailCamera.y - (_detailCamera.viewSize / 2.0f);
    float t = _detailCamera.y + (_detailCamera.viewSize / 2.0f);
    mat4x4_ortho(projection, l, r, b, t, 1.0f, -1.0f);

    mat4x4_mul(transform, projection, _detailTransform.modelMatrix);

    glUniformMatrix4fv(_detailView.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *) transform);

    GLint loc_colorMask = glGetUniformLocation(_detailView.program, "colorMask");
    GLint loc_weight = glGetUniformLocation(_detailView.program, "weight");
    glUniform1i(loc_colorMask, ColorMask_TEXTURE | ColorMask_PRIMITIVE);
    glUniform4f(loc_weight, 1.0f, 1.0f, 1.0f, 1.0f);

    int *vp = info.detailViewport;
    glViewport(vp[0], vp[1], vp[2], vp[3]);
    glScissor(vp[0], vp[1], vp[2], vp[3]);

    glDrawBuffer(GL_BACK);
    glDisable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    if (_detailView.wireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    glDrawArrays(GL_TRIANGLES, 0, _detailView.count);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glBindVertexArray(_detailView.borderVao);
    glUniform1i(loc_colorMask, ColorMask_EMPTY);

    // draw boundary
    glUniform4f(loc_weight, 0.2f, 0.9f, 0.2f, 1.0f);
    glDrawArrays(GL_LINES, 0, _detailView.borderCount);

    // draw seams
    glUniform4f(loc_weight, 0.9f, 0.9f, 0.9f, 1.0f);
    glDrawArrays(GL_LINES, _detailView.borderCount, _detailView.seamCount);

    glBindVertexArray(0);
    CheckGLError();
}

void MeshViewer::DrawViews()
{
    Draw3DView();
    DrawTextureView();
    if (selectedRegions.size() > 0) {
        DrawDetailView();
    }
}


void MeshViewer::Run()
{
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_VISIBLE, GL_TRUE);

    _window = glfwCreateWindow(1280, 720, "Mesh viewer", NULL, NULL);
    if (!_window) {
        cout << "Failed to create window or context" << endl;
        std::exit(-1);
    }

    // Set user pointer
    glfwSetWindowUserPointer(_window, static_cast<void*>(this));

    int width, height;
    glfwGetFramebufferSize(_window, &width, &height);

    // Set callbacks
    glfwSetMouseButtonCallback(_window, MeshViewer::MouseButtonCallback);
    glfwSetCursorPosCallback(_window, MeshViewer::CursorPositionCallback);
    glfwSetScrollCallback(_window, MeshViewer::ScrollCallback);
    glfwSetKeyCallback(_window, ImGui_ImplGlfwGL3_KeyCallback);
    glfwSetCharCallback(_window, ImGui_ImplGlfwGL3_CharCallback);
    glfwSetFramebufferSizeCallback(_window, MeshViewer::FramebufferSizeCallback);
    FramebufferSizeCallback(_window, width, height);

    glfwMakeContextCurrent(_window);
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err) {
        cout << "Failed to initialize glew: " << glewGetErrorString(err) << endl;
        std::exit(-1);
    }

    glGetError(); // suppress possible error on glew init

    glfwSwapInterval(1);

    ImGui_ImplGlfwGL3_Init(_window);

    GLint loc_tex0;
    _perspectiveView.program = CompileShaders(vs_text_3D, fs_text_texture);
    _perspectiveView.selection.program = CompileShaders(vs_text_3D, fs_text_texture);
    _textureView.program = CompileShaders(vs_text_texture, fs_text_texture);
    _textureView.highlight.program = CompileShaders(vs_text_texture, fs_text_texture);
    _detailView.program = CompileShaders(vs_text_texture, fs_text_texture);

    loc_tex0 = glGetUniformLocation(_perspectiveView.program, "tex0");
    glUseProgram(_perspectiveView.program);
    glUniform1i(loc_tex0, 0);

    loc_tex0 = glGetUniformLocation(_perspectiveView.selection.program, "tex0");
    glUseProgram(_perspectiveView.selection.program);
    glUniform1i(loc_tex0, 0);

    loc_tex0 = glGetUniformLocation(_textureView.program, "tex0");
    glUseProgram(_textureView.program);
    glUniform1i(loc_tex0, 0);

    loc_tex0 = glGetUniformLocation(_textureView.highlight.program, "tex0");
    glUseProgram(_textureView.highlight.program);
    glUniform1i(loc_tex0, 0);

    loc_tex0 = glGetUniformLocation(_detailView.program, "tex0");
    glUseProgram(_detailView.program);
    glUniform1i(loc_tex0, 0);

    glUseProgram(0);

    CheckGLError();

    InitBuffers();
    SetupViews();

    glEnable(GL_SCISSOR_TEST);

    bool show_test_window = false;

    while (!glfwWindowShouldClose(_window)) {
        glDrawBuffer(GL_BACK);
        CheckGLError();
        glScissor(0, 0, info.width, info.height);
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        CheckGLError();

        glfwPollEvents();
        CheckGLError();
        ImGui_ImplGlfwGL3_NewFrame();
        CheckGLError();
        ManageImGuiState();
        CheckGLError();
        UpdateTransforms();
        CheckGLError();
        DrawViews();

        if (show_test_window) {
            ImGui::SetNextWindowPosCenter(ImGuiCond_FirstUseEver);
            ImGui::ShowTestWindow(&show_test_window);
        }
        ImGui::Render();

        glfwSwapBuffers(_window);
    }

    ImGui_ImplGlfwGL3_Shutdown();

    glfwDestroyWindow(_window);
}

void MeshViewer::ManageImGuiState()
{
    static bool showInfoArea = true;

    // distortion display
    static bool distortionFromTexture = false;
    static DistortionMetric::Type distortion[] = {
        DistortionMetric::Area,
        DistortionMetric::Angle,
    };
    static int distortionIndex = 0;
    static int activeDistIndex = -1;

    // parameterization strategy
    static DirectParameterizer dirparameterizer[] = { DCP, FixedBorderBijective };
    static EnergyType optimizer[] = { SymmetricDirichlet};
    static DescentType descent[] = { Gradient, LimitedMemoryBFGS, ScalableLocallyInjectiveMappings, CompositeMajorization };
    static ParameterizationGeometry geometry[] = { Model, Texture };

    static int parameterizerInUse = 0;
    static int optimizerInUse = 0;
    static int descentTypeInUse = 2;
    static int geometryInUse = 0;

    // data to update
    bool updateTexcoord = false;
    bool updateColor = false;

    // draw main controls
    {

        ImGui::SetNextWindowPos(ImVec2{0, 0}, ImGuiCond_Once);

        ImGui::Begin("Controls", nullptr, 0);

        ImGui::Checkbox("Info area##checkbox", &showInfoArea);
        ImGui::SameLine();
        ImGui::Checkbox("Wireframe", &_detailView.wireframe);

        ImGui::Separator();

        if (ImGui::Button("Select next greedy merge")) {
            ClearSelection();
            if (gmHasNextEdge()) {
                const auto& next = gmPeekNextEdge();
                std::cout << "Next edge has weight " << next.second << std::endl;
                Select(next.first.a->id);
                Select(next.first.b->id);
            }
            else {
                std::cout << "No next edge available" << std::endl;
            }
        }

        bool disabled = false;
        if (primaryCharts.size() < 2) {
            ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
            ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
            disabled = true;
        }
        if (ImGui::Button("Merge current selection")) {
            std::cout << "Merging selected charts" << std::endl;
            std::vector<ChartHandle> cm;
            for (auto& entry : primaryCharts) {
                cm.push_back(graph->GetChart(entry.first));
            }
            std::pair<int,ChartHandle> result = gmCollapse(cm.begin(), cm.end());
            if (result.first == gm->Collapse_OK) {
                ClearSelection();
                Select(result.second->id);
                std::cout << "Collapse succeeded" << std::endl;
            } else {
                if (result.first == gm->Collapse_ERR_DISCONNECTED) {
                    std::cout << "Cannot merge, selected charts are disconnected" << std::endl;
                } else if (result.first == gm->Collapse_ERR_UNFEASIBLE) {
                    std::cout << "Cannot merge, result is unfeasible" << std::endl;
                } else {
                    assert(0 && "GraphManager::Collapse() exit status");
                }
            }
        }
        if (disabled) {
            ImGui::PopItemFlag();
            ImGui::PopStyleVar();
        }

        if (ImGui::Button("Close islands")) {
            ClearSelection();
            gmClose();
        }

        if (ImGui::Button("Invoke greedy algorithm")) {
            ClearSelection();
            RecomputeSegmentation(*gm, 20, minRegionSize);
        }

        static bool retry = true;
        static float tau = 0.0f;
        //bool clickPack = ImGui::Button("Parameterize graph");
        bool clickPack = ImGui::Button("Pack the atlas");
        ImGui::SameLine();
        ImGui::Checkbox("Retry", &retry);
        ImGui::InputFloat("tau", &tau, 0.002f, 1.0f);
        if (tau < 0) tau = 0;
        if (tau > 1) tau = 1;
        if (clickPack) {
            ClearSelection();
            if (graph->MergeCount() > 0) {
                int c = ParameterizeGraph(*gm, strategy, retry ? tau : -1);
                if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;
                PackingOptions opts = { RasterizationBasedPacker::Parameters::CostFuncEnum::MinWastedSpace, true, true, true };
                Pack(gm->Graph(), opts);
                _currentTexture->Release(0);
                _currentTexture = RenderTexture(graph->mesh, graph->textureObject, true, InterpolationMode::Linear, _window);
                for (auto& f : graph->mesh.face) {
                    for (int i = 0; i < 3; ++i) {
                        f.WT(i).P()[0] *= _currentTexture->TextureWidth(0);
                        f.WT(i).P()[1] *= _currentTexture->TextureHeight(0);
                    }
                }
                ScaleTextureCoordinatesToImage(graph->mesh, _currentTexture);
                updateTexcoord = true;
                if (activeDistIndex != -1) {
                    ParameterizationGeometry distGeom = distortionFromTexture ? Texture : Model;
                    graph->MapDistortion(distortion[activeDistIndex], distGeom);
                    updateColor = true;
                }

           } else {
               std::cout << "No merges, nothing to do" << std::endl;
           }
        }

        if (primaryCharts.size() == 1) {
            if (ImGui::Button("Save current chart")) {
                Mesh pm;
                auto chart = graph->GetChart(primaryCharts.begin()->first);
                Box2d b = chart->UVBox();

                auto f = [&pm, &b](typename Mesh::FacePointer fptr) {
                    auto f = tri::Allocator<Mesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
                    for (int i = 0; i < 3; ++i) {
                        f->WT(i) = fptr->WT(i);
                        f->WT(i).P() = (f->WT(i).P() - b.min) / std::max(b.DimX(), b.DimY());
                    }
                };

                std::for_each(chart->fpVec.begin(), chart->fpVec.end(), f);

                tri::Clean<Mesh>::RemoveDuplicateVertex(pm);
                tri::Allocator<Mesh>::CompactEveryVector(pm);
                tri::UpdateTopology<Mesh>::FaceFace(pm);

                tri::io::ExporterOBJ<Mesh>::Save(pm, "chart.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
            }
        }

        ImGui::Separator();

        static int selId = 0;
        ImGui::InputInt("##SelectId", &selId, 1, 1);
        ImGui::SameLine();
        if (ImGui::Button("S") && graph->GetChart(selId) != nullptr) {
            Select(selId);
        }

        ImGui::Separator();

        static char exportFileName[256] = "";
        if (*exportFileName == 0) std::snprintf(exportFileName, 256, "%s", fileName.c_str());
        if (ImGui::Button("Export mesh")) {
            ImGui::OpenPopup("Export mesh...");
        }
        if (ImGui::BeginPopupModal("Export mesh...", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::InputText("file name", exportFileName, 256);
            if (ImGui::Button("Export", ImVec2(120,0))) {
                Mesh& m = graph->mesh;
                if (SaveMesh(m, exportFileName, _currentTexture, true) == false) {
                    std::cout << "Model not saved correctly" << std::endl;
                }
                ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120,0))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }

        ImGui::Separator();

        ImGui::Text("Geometry");
        ImGui::RadioButton("Vertex position", &geometryInUse, 0);
        ImGui::RadioButton("Texture coords", &geometryInUse, 1);

        static bool padBoundaries = true;
        static bool applyCut = true;
        static bool scaffold = false;
        ImGui::Text("Parameterizer");
        ImGui::RadioButton("Discrete Conformal", &parameterizerInUse, 0);
        ImGui::RadioButton("Circular border bijective", &parameterizerInUse, 1);
        if (parameterizerInUse == 1) {
            ImGui::Checkbox("Pad inner boundaries", &padBoundaries);
            ImGui::Checkbox("Apply cut", &applyCut);
            ImGui::Checkbox("Scaffold", &scaffold);
        }

        //ImGui::Text("Descent method energy");
        //ImGui::RadioButton("Area preserving", &optimizerInUse, 0);
        //ImGui::RadioButton("Symmetric Dirichlet", &optimizerInUse, 1);
        //ImGui::RadioButton("MIPS", &optimizerInUse, 2);

        ImGui::Text("Descent method");
        ImGui::RadioButton("Gradient descent", &descentTypeInUse, 0);
        ImGui::RadioButton("LBFGS", &descentTypeInUse, 1);
        ImGui::RadioButton("SLIM", &descentTypeInUse, 2);
        ImGui::RadioButton("CM", &descentTypeInUse, 3);

        strategy.directParameterizer = dirparameterizer[parameterizerInUse];
        strategy.energy = optimizer[optimizerInUse];
        strategy.geometry = geometry[geometryInUse];
        strategy.descent = descent[descentTypeInUse];
        strategy.padBoundaries = padBoundaries;
        strategy.applyCut = applyCut;
        strategy.scaffold = scaffold;

        ImGui::Text("Max descent iterations");
        ImGui::InputInt("##Optimizer iterations", &strategy.optimizerIterations, 1, 100);
        if (strategy.optimizerIterations < 0) strategy.optimizerIterations = 0;
        if (ImGui::Button("Reselect")) {
            std::unordered_map<RegionID,int> savedPrimary = primaryCharts;
            int numIter = strategy.optimizerIterations; strategy.optimizerIterations = 0;
            ClearSelection();
            std::size_t numPrimary = savedPrimary.size();
            for (auto& entry : savedPrimary) {
                if (numPrimary == 1) strategy.optimizerIterations = numIter;
                Select(entry.first);
                numPrimary--;
            }
        }

        ImGui::Separator();

        ImGui::Text("Distortion metric");
        bool clicked = ImGui::Checkbox("Target original texture", &distortionFromTexture);
        ImGui::RadioButton("Area Distortion", &distortionIndex, 0);
        ImGui::RadioButton("Angle Distortion", &distortionIndex, 1);
        if (distortionIndex != activeDistIndex || clicked) {
            ParameterizationGeometry distGeom = distortionFromTexture ? Texture : Model;
            graph->MapDistortion(distortion[distortionIndex], distGeom);
            activeDistIndex = distortionIndex;
            updateColor = true;
        }

        ImGui::Separator();

        ImGui::Text("Color source");
        static bool mixTexture = true;
        static bool mixDistortion = false;
        static bool mixCheckboard = false;
        ImGui::Checkbox("Texture", &mixTexture);
        ImGui::Checkbox("Distortion", &mixDistortion);
        ImGui::Checkbox("Checkboard", &mixCheckboard);
        int colorMask = ColorMask_EMPTY;
        if (mixTexture) colorMask |= ColorMask_TEXTURE;
        if (mixDistortion) colorMask |= ColorMask_PRIMITIVE;
        if (mixCheckboard) colorMask |= ColorMask_CHECKBOARD;
        _perspectiveView.colorMask = colorMask;

        ImGui::End();
    }

    // info area
    {
        if (showInfoArea) {
            ImGui::Begin("Info area", &showInfoArea);
            if (primaryCharts.size() == 1) {
                // display info about the chart
                auto chart = graph->GetChart(primaryCharts.begin()->first);
                ImGui::Text("Chart %lu (%lu faces, %lu adjacencies)", chart->id, chart->FN(), chart->NumAdj());
                ImGui::Text("Aggregate count: %d", chart->numMerges + 1);
                ImGui::Text("Area 3D: %.4f", chart->Area3D());
                ImGui::Text("Area UV: %.4f | Border UV: %.6f", chart->AreaUV(), chart->BorderUV());
                ImGui::Text("Distortion range: %.4f , %.4f", chart->minMappedFaceValue, chart->maxMappedFaceValue);
            } else {
                // display info about the whole graph
                ImGui::Text("Filename: %s", graph->mesh.name.c_str());
                ImGui::Text("Parameterization charts: %lu", graph->Count());
                ImGui::Text("Area 3D: %.4f", graph->Area3D());
                ImGui::Text("Area UV: %.4f | Border UV: %.6f", graph->AreaUV(), graph->BorderUV());
                auto range = graph->DistortionRange();
                ImGui::Text("Distortion range: %.4f , %.4f", range.first, range.second);
            }
            ImGui::End();
        }
    }

    // selection widget
    /*
    {
        if (selectedRegions.size() > 0) {
            ImGui::Begin("Selection widget", nullptr, 0);
            int i = 0;
            std::vector<RegionID> selected;
            for (const auto& entry : selectedRegions) {
                RegionID id = entry.first;
                const SelectionBufferInfo& sbi = entry.second;
                if (sbi.referenceCount > 0) {
                    ImVec4 bg = (primaryCharts.count(id) == 1) ? ImVec4(0.4f, 0.4f, 0.4f, 1.0f) : ImVec4(0.1f, 0.1f, 0.1f, 1.0f);
                    if (ImGui::ImageButton((void *)(long unsigned)sbi.texicon, ImVec2(64, 64), ImVec2(0, 1), ImVec2(1, 0), -1, bg)) {
                        selected.push_back(id);
                    }
                    if ((++i % 4) != 0) ImGui::SameLine();
                }
            }
            for (auto id : selected) Select(id);
            ImGui::End();
        }
    }
    */

    // shell parameterization controls
    bool shellChanged = false;
    static int shellColor = 0;
    if (selectedRegions.size() > 0) {
        ImGui::Begin("Shell parameterization", nullptr, 0);
        static const char *colorOptions[] = {
            "None",
            "Energy value",
            "Gradient",
            "Descent direction",
            "Var(Local Grad)",
            "Conformal scaling"
        };
        if (ImGui::Combo("Shell color", &shellColor, colorOptions, IM_ARRAYSIZE(colorOptions))) {
            shellChanged = true;
        }
        if (ImGui::Button("Parameterize##ParamObject")) {
            parameterizer->Parameterize();
            shellChanged = true;
        }
        if (ImGui::Button("Reset")) {
            parameterizer->Reset();
            shellChanged = true;
        }
        static int iternum = 0;
        ImGui::InputInt("##SelectId", &iternum, 1, 1);
        if (iternum < 1) iternum = 1;
        ImGui::SameLine();
        ImGui::Text("Iteration count");
        if (ImGui::Button("Iterate")) {
            IterationInfo info;
            for (int i = 0; i < iternum; ++i) {
                std::cout << "Iteration " << i << std::endl;
                info = parameterizer->Iterate();
            }
            std::cout << "Energy = " << info.energyVal << ", DeltaE = " << info.energyDiff << ", Gradient norm = " << info.gradientNorm << std::endl;
            shellChanged = true;
        }
        ImGui::SameLine();
        ImGui::Text("Iterations: %d", parameterizer->IterationCount());

        if (ImGui::Button("Place cut")) {
            parameterizer->PlaceCut();
            shellChanged = true;
        }


        ImGui::Separator();
        static int ncones = 0;
        ImGui::InputInt("##ncones", &ncones, 1, 1);
        if (ncones < 1) ncones = 1;


        if (ImGui::Button("Place cut with cone singularities")) {
            if (parameterizer->PlaceCutWithConesUntilThreshold(0.01)) {
                //parameterizer->InitializeSolution();
                shellChanged = true;
            } else {
                std::cout << "Unable to place cut" << std::endl;
            }
        }
        ImGui::Separator();

        if (ImGui::Button("Fix uvs")) {
            if (primaryCharts.size() != 1) {
                std::cout << "Unable to fix uv for non-unique primary chart" << std::endl;
            } else {
                parameterizer->SyncChart();
                RegionID selId = primaryCharts.begin()->first;
                ClearSelection();
                Select(selId);
                parameterizer->ForceWarmStart();
                shellChanged = true;
            }
        }

        if (ImGui::Button("Split chart")) {
            if (primaryCharts.size() != 1) {
                std::cout << "Unable to fix uv for non-unique primary chart" << std::endl;
            } else {
                auto id = primaryCharts.begin()->first;
                ClearSelection();
                std::vector<ChartHandle> splitCharts;
                gm->Split(id, splitCharts);
                std::vector<ChartHandle> newCharts;
                RecoverFromSplit(splitCharts, *gm, newCharts, true);
                for (auto& c : newCharts)
                    std::cout << "Chart " << c->id << " is new" << std::endl;
            }


        }

        ImGui::End();
    }

    if (shellChanged) {
        switch (shellColor) {
        case 0:
            if (shellColorMode != NONE) {
                parameterizer->ClearShellFaceColor();
                shellColorMode = NONE;
            }
            break;
        case 1:
            parameterizer->MapEnergyToShellFaceColor();
            shellColorMode = FACE;
            break;
        case 2:
            parameterizer->MapEnergyGradientToShellVertexColor();
            shellColorMode = VERTEX;
            break;
        case 3:
            parameterizer->MapDescentDirectionToShellVertexColor();
            shellColorMode = VERTEX;
            break;
        case 4:
            parameterizer->MapLocalGradientVarianceToShellVertexColor();
            shellColorMode = VERTEX;
            break;
        case 5:
            parameterizer->MapConformalScalingFactorsToShellVertexColor();
            shellColorMode = VERTEX;
            break;
        default: break;
        }

        UpdateDetailBuffers();
    }

    if (updateTexcoord || updateColor) {
        glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.mesh);
        float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        for (const auto &f : graph->mesh.face) {
            for (int i = 0; i < 3; ++i) {
                buffptr += 3;
                if (updateTexcoord) {
                    *buffptr++ = f.cWT(i).U() / (double) _currentTexture->TextureWidth(0);
                    *buffptr++ = f.cWT(i).V() / (double) _currentTexture->TextureHeight(0);
                }
                else buffptr += 2;
                if (updateColor) {
                    unsigned char *colorptr = (unsigned char *) buffptr;
                    *colorptr++ = f.cC()[0];
                    *colorptr++ = f.cC()[1];
                    *colorptr++ = f.cC()[2];
                    *colorptr++ = f.cC()[3];
                }
                buffptr++;
            }
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}


// GraphManager proxy interface

void MeshViewer::gmClose()
{
    int count = gm->CloseMacroRegions(minRegionSize);
    if (count > 0) std::cout << "Close operation merged " << count << " regions" << std::endl;
}

bool MeshViewer::gmHasNextEdge()
{
    return gm->HasNextEdge();
}

const std::pair<GraphManager::Edge,double>& MeshViewer::gmPeekNextEdge()
{
    return gm->PeekNextEdge();
}

void MeshViewer::gmRemoveNextEdge()
{
    gm->RemoveNextEdge();
}

GraphManager::ChartHandle MeshViewer::gmCollapse(const GraphManager::Edge& e)
{
    return gm->Collapse(e);
}

template <typename ChartInputIterator>
std::pair<int,GraphManager::ChartHandle> MeshViewer::gmCollapse(ChartInputIterator first, ChartInputIterator last)
{
    return gm->Collapse(first, last);
}
