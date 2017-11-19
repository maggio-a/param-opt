#include <iostream>
#include <cstdio>

#include <QImage>

#include <vcg/space/ray3.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/intersection2.h>



#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/export.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <wrap/gui/trackball.h>

#include <imgui.h>
#include <imgui_internal.h> // to disable elements
#include <imgui_glfw_gl3/imgui_impl_glfw_gl3.h>

#include "mesh.h"
#include "mesh_viewer.h"
#include "gl_util.h"
#include "optimizer.h"
#include "texture_rendering.h"

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
    "                                                                              \n"
    "out vec2 uv;                                                                  \n"
    "out vec4 dcol;                                                                \n"
    "                                                                              \n"
    "void main(void)                                                               \n"
    "{                                                                             \n"
    "    uv = texcoord;                                                            \n"
    "    dcol = primitiveColor;                                                    \n"
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
    "        color *= texture(tex0, uv);                            \n"
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
    case GLFW_MOUSE_BUTTON_RIGHT: trackButton = vcg::Trackball::BUTTON_RIGHT; break;
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

    viewer->info.perspectiveViewport[0] = (int) (0.2f * borderToSplitWidth);
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

MeshViewer::MeshViewer(std::shared_ptr<MeshGraph> meshParamData_, std::size_t minRegionSize_, const std::string& fileName_)
    : meshParamData{meshParamData_}, _currentTexture{meshParamData_->textureObject}, gm{std::make_shared<GraphManager>(meshParamData_)},
      fileName{fileName_}, minRegionSize{minRegionSize_}, _textureCamera{}, _detailCamera{}
{
    std::size_t numRegions = meshParamData->Count();
    regionColors.reserve(numRegions);
    for (const auto& c : meshParamData->charts) {
        auto color = vcg::Color4f::Scatter(20, c.first % 20, 0.95f);
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

    Point2f p{_textureCamera.x + dx, _textureCamera.y + dy};

    const Mesh& m = meshParamData->mesh;
    const MeshFace *fp = nullptr;
    for (auto &f : m.face) {
        vcg::Triangle2<float> uvFace{f.cWT(0).P(), f.cWT(1).P(), f.cWT(2).P()};
        if (vcg::IsInsideTrianglePoint(uvFace, p)) {
            fp = &f;
            break;
        }
    }

    if (fp != nullptr) {
        auto CCIDh = tri::Allocator<Mesh>::GetPerFaceAttribute<RegionID>(meshParamData->mesh, "ConnectedComponentID");
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
        glDeleteBuffers(1, &_vertexBuffers.detailBorder);

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
    }
}

void MeshViewer::Select(const RegionID id)
{
    if (primaryCharts.count(id) == 1) return;
    else UpdateSelection(id);
}

void MeshViewer::UpdateSelection(const RegionID id)
{
    // Allocate new buffers if necessary
    std::set<ChartHandle> newCharts;
    std::size_t  newElements = 0;
    if (selectedRegions.count(id) == 0) {
        newCharts.insert(meshParamData->GetChart(id));
        newElements += meshParamData->GetChart(id)->FN();
    }
    else {
        selectedRegions[id].referenceCount++;
    }
    for (auto chart : meshParamData->GetChart(id)->adj) {
        if (selectedRegions.count(chart->id) == 0) {
            newCharts.insert(chart);
            newElements += chart->FN();
        } else {
            selectedRegions[chart->id].referenceCount++;
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
                    *buffptr++ = fptr->cWT(i).U();
                    *buffptr++ = fptr->cWT(i).V();
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

            Point2f center = chart->UVBox().Center();
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
            vcg::Box2f box = chart->UVBox();

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

        CheckGLError();
    }

    // make the selected chart primary
    primaryCharts[id] = 1;

    if (_vertexBuffers.detail != 0) {
        glDeleteBuffers(1, &_vertexBuffers.detail);
        _vertexBuffers.detail = 0;
    }

    std::vector<ChartHandle> charts;
    for (auto& entry : primaryCharts) charts.push_back(meshParamData->GetChart(entry.first));

    auto status = gm->CollapseAllowed(charts.begin(), charts.end());

    auto CCIDh = tri::Allocator<Mesh>::FindPerFaceAttribute<RegionID>(meshParamData->mesh, "ConnectedComponentID");
    assert(tri::Allocator<Mesh>::IsValidHandle<RegionID>(meshParamData->mesh, CCIDh));

    if (status.first == gm->Collapse_OK) { // the regions can be parameterized together
        // Note that this passes even if there is only a single primary chart

        struct FaceData { RegionID id; TexCoordStorage wt; };
        std::unordered_map<Mesh::FacePointer,FaceData> fd;
        constexpr RegionID tempID = 0xffffffff;

        auto aggregate = std::make_shared<FaceGroup>(meshParamData->mesh, tempID);

        for (auto& c : charts) {
            for (auto fptr : c->fpVec) {
                TexCoordStorage wt;
                for (int i = 0; i < 3; ++i) {
                    wt.tc[i] = fptr->WT(i);
                }
                fd.insert(std::make_pair(fptr, FaceData{CCIDh[fptr], wt}));
                CCIDh[fptr] = tempID;
                aggregate->AddFace(fptr);
            }
        }
        // Parameterize the aggregate chart, build the vertex buffer and restore the original state

        ParameterizeChartFromInitialTexCoord(meshParamData->mesh, aggregate);

        glUseProgram(_detailView.program);

        std::vector<float> borderVertexData;

        glGenVertexArrays(1, &_detailView.vao);
        glBindVertexArray(_detailView.vao);
        glGenBuffers(1, &_vertexBuffers.detail);
        glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.detail);
        glBufferData(GL_ARRAY_BUFFER, 12 * aggregate->FN() * sizeof(float), NULL, GL_STATIC_DRAW);
        float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        for (auto fptr : aggregate->fpVec) {
            for (int i = 0; i < 3; ++i) {
                if (face::IsBorder(*fptr, i) || CCIDh[fptr->FFp(i)] != tempID) {
                    borderVertexData.push_back(fptr->WT(i).P().X());
                    borderVertexData.push_back(fptr->WT(i).P().Y());
                    borderVertexData.push_back(fptr->WT((i+1)%3).P().X());
                    borderVertexData.push_back(fptr->WT((i+1)%3).P().Y());
                }
                *buffptr++ = fptr->WT(i).P().X();
                *buffptr++ = fptr->WT(i).P().Y();
                *buffptr++ = fd[fptr].wt.tc[i].P().X();
                *buffptr++ = fd[fptr].wt.tc[i].P().Y();
            }
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);

        _detailView.count = (GLsizei) (aggregate->FN() * 3);

        _detailView.attributes.loc_position = glGetAttribLocation(_detailView.program, "position");
        glVertexAttribPointer(_detailView.attributes.loc_position, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
        glEnableVertexAttribArray(_detailView.attributes.loc_position);

        _detailView.attributes.loc_texcoord = glGetAttribLocation(_detailView.program, "texcoord");
        glVertexAttribPointer(_detailView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (const GLvoid *) (2*sizeof(float)));
        glEnableVertexAttribArray(_detailView.attributes.loc_texcoord);

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glGenVertexArrays(1, &_detailView.borderVao);
        glBindVertexArray(_detailView.borderVao);
        glGenBuffers(1, &_vertexBuffers.detailBorder);
        glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.detailBorder);
        glBufferData(GL_ARRAY_BUFFER, borderVertexData.size() * sizeof(float), &borderVertexData[0], GL_STATIC_DRAW);
        _detailView.borderCount = (GLsizei) (borderVertexData.size()/2);

        glVertexAttribPointer(_detailView.attributes.loc_position, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(_detailView.attributes.loc_position);

        glVertexAttribPointer(_detailView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(_detailView.attributes.loc_texcoord);

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        SetupDetailView(aggregate);

        // restore state
        for (auto fptr : aggregate->fpVec) {
            CCIDh[fptr] = fd[fptr].id;
            for (int i = 0; i < 3; ++i) {
                fptr->WT(i) = fd[fptr].wt.tc[i];
            }
        }
    } else {
        /// todo check error code
        std::cout << "Warning: current selection cannot be parameterized" << std::endl;
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

    Ray3f selectionRay {
        Point3f{modelEye[0], modelEye[1], modelEye[2]},
        Point3f{modelRay[0], modelRay[1], modelRay[2]}
    };

    // Perform intersection tests
    const Mesh& m = meshParamData->mesh;
    *fp = nullptr;
    float tMin = std::numeric_limits<float>::max();
    for (auto &f : m.face) {
        float t, uface, vface;
        if (vcg::IntersectionRayTriangle(selectionRay, f.cP(0), f.cP(1), f.cP(2), t, uface, vface) && t < tMin) {
            *fp = &f;
            tMin = t;
            u = uface;
            v = vface;
        }
    }
    return tMin != std::numeric_limits<float>::max();
}

void MeshViewer::PerspectivePick()
{
    Mesh::ConstFacePointer fp;
    float u, v;

    //if (fp != nullptr) {
    if (IntersectionMouseRayModel(&fp, u, v)) {
        auto CCIDh = tri::Allocator<Mesh>::GetPerFaceAttribute<RegionID>(meshParamData->mesh, "ConnectedComponentID");
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
            Point3f centerPoint = fp->cP(0) * (1.0f - u - v) + fp->cP(1) * u + fp->cP(2) * v;
            Matrix44f transform = trackball.Matrix();
            trackball.Translate(-(transform*centerPoint));
        }
    }
}

GLuint MeshViewer::CompileShaders(const GLchar **vs_text, const GLchar **fs_text)
{
    GLint status;
    char infoLog[1024] = {0};

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, vs_text, NULL);
    glCompileShader(vs);
    glGetShaderInfoLog(vs, 1024, NULL, infoLog);
    if (*infoLog) {
        std::cout << infoLog << std::endl;
        memset(infoLog, 0, 1024);
    }
    glGetShaderiv(vs, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        std::cout << "Vertex shader compilation failed" << std::endl;
    }

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, fs_text, NULL);
    glCompileShader(fs);
    glGetShaderInfoLog(fs, 1024, NULL, infoLog);
    if (*infoLog) {
        std::cout << infoLog << std::endl;
        memset(infoLog, 0, 1024);
    }
    glGetShaderiv(fs, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        std::cout << "Fragment shader compilation failed" << std::endl;
    }

    GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);
    glValidateProgram(program);
    glGetProgramInfoLog(program, 1024, NULL, infoLog);
    if (*infoLog) {
        std::cout << infoLog << std::endl;
    }
    glGetProgramiv(program, GL_LINK_STATUS, &status);
    if (status == GL_FALSE) {
        std::cout << "Shader program link failed" << std::endl;
    }

    glDeleteShader(vs);
    glDeleteShader(fs);

    CheckGLError();

    return program;
}

void MeshViewer::InitBuffers()
{
    const Mesh& m = meshParamData->mesh;

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
            *buffptr++ = f.cWT(i).U();
            *buffptr++ = f.cWT(i).V();
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

    // Load texture data
    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind();

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
    Mesh& m = meshParamData->mesh;

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
    _perspectiveView.selection.uniforms.loc_weight = glGetUniformLocation(_perspectiveView.selection.program, "weight");

    // Setup texture view

    _textureCamera.Reset();
    _textureView.uniforms.loc_projection = glGetUniformLocation(_textureView.program, "projectionMatrix");
    _textureView.highlight.uniforms.loc_projection = glGetUniformLocation(_textureView.highlight.program, "projectionMatrix");
    _textureView.highlight.uniforms.loc_primitiveColor = glGetUniformLocation(_textureView.highlight.program, "primitiveColor");
}

void MeshViewer::SetupDetailView(ChartHandle chart)
{
    Box2f bbox = chart->UVBox();

    _detailCamera.Reset();
    _detailView.uniforms.loc_projection = glGetUniformLocation(_detailView.program, "projectionMatrix");

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
    _currentTexture->Bind();

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
            glBindVertexArray(selectionVao[sbi.bufferIndex]);
            if (primaryCharts.count(id) == 1) {
                float one[] = {1.0f, 1.0f, 1.0f, 1.0f};
                glUniform4fv(_perspectiveView.selection.uniforms.loc_weight, 1, one);
            } else {
                glUniform4fv(_perspectiveView.selection.uniforms.loc_weight, 1, regionColors[id].V());
            }
            glDrawArrays(GL_TRIANGLES, sbi.first, sbi.count);
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

        glDrawArrays(GL_TRIANGLES, 0, meshParamData->mesh.FN() * 3);
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

        glDrawArrays(GL_TRIANGLES, 0, meshParamData->mesh.FN() * 3);
        glBindVertexArray(0);
    }

    CheckGLError();
}

void MeshViewer::DrawTextureView()
{
    glBindVertexArray(_textureView.vao);
    glUseProgram(_textureView.program);

    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind();

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

    glDrawArrays(GL_TRIANGLES, 0, meshParamData->mesh.FN() * 3);

    if (selectedRegions.size() > 0) {
        glUseProgram(_textureView.highlight.program);

        glUniformMatrix4fv(_textureView.highlight.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *) projection);

        GLint loc_colorMask = glGetUniformLocation(_textureView.highlight.program, "colorMask");
        glUniform1i(loc_colorMask, ColorMask_PRIMITIVE); // set primitive color as color source

        for (const auto& sel : selectedRegions) {
            RegionID id = sel.first;
            const SelectionBufferInfo& sbi = sel.second;
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

    glBindVertexArray(0);

    CheckGLError();
}

// as of now it shows the atlas of primary charts
void MeshViewer::DrawDetailView()
{
    glBindVertexArray(_detailView.vao);
    glUseProgram(_detailView.program);

    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind();

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
    glUniform1i(loc_colorMask, ColorMask_TEXTURE);
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
    glUniform4f(loc_weight, 0.3f, 1.0f, 0.3f, 1.0f);
    glDrawArrays(GL_LINES, 0, _detailView.borderCount);

    glBindVertexArray(0);
    CheckGLError();
}

void MeshViewer::DrawViews()
{
    Draw3DView();
    DrawTextureView();
    if (selectedRegions.size() > 0) DrawDetailView();
}


void MeshViewer::Run()
{
    if (!glfwInit()) {
        std::cout << "Failed to initialize glfw" << std::endl;
        std::exit(-1);
    }

    glfwSetErrorCallback(ErrorCallback);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    _window = glfwCreateWindow(980, 720, "Mesh viewer", NULL, NULL);
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
    glUseProgram(_perspectiveView.program);
    glUniform1i(loc_tex0, 0);

    loc_tex0 = glGetUniformLocation(_detailView.program, "tex0");
    glUseProgram(_detailView.program);
    glUniform1i(loc_tex0, 0);

    glUseProgram(0);

    InitBuffers();
    SetupViews();

    glEnable(GL_SCISSOR_TEST);

    bool show_test_window = false;

    while (!glfwWindowShouldClose(_window)) {

        glDrawBuffer(GL_BACK);
        glScissor(0, 0, info.width, info.height);
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glfwPollEvents();

        ImGui_ImplGlfwGL3_NewFrame();

        ManageImGuiState();
        UpdateTransforms();
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
    glfwTerminate();

}

void MeshViewer::ManageImGuiState()
{
    static bool showInfoArea = true;

    static DistortionWedge::DistType distortion[] = {
        DistortionWedge::AreaDist,
        DistortionWedge::EdgeDist,
        DistortionWedge::AngleDist,
    };
    static int distortionIndex = 0;
    static int activeDistIndex = -1;

    bool updateTexcoord = false;
    bool updateColor = false;
    // draw main controls
    {

        ImGui::SetNextWindowPos(ImVec2{0, 0}, ImGuiCond_Once);

        ImGui::Begin("Controls", nullptr, 0);

        ImGui::Checkbox("Info area##checkbox", &showInfoArea);

        ImGui::Checkbox("Toggle wireframe", &_detailView.wireframe);

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
                cm.push_back(meshParamData->GetChart(entry.first));
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
            ReduceTextureFragmentation_NoPacking(*gm, minRegionSize);
        }

        if (ImGui::Button("Pack the atlas")) {
            ClearSelection();
            if (meshParamData->MergeCount() > 0) {
                int c = ParameterizeGraph(meshParamData);
                if (c > 0) std::cout << "WARNING: " << c << " regions were not parameterized correctly" << std::endl;
                updateTexcoord = true;
                if (activeDistIndex != -1) {
                    meshParamData->MapDistortion(distortion[activeDistIndex]);
                    updateColor = true;
                }
                _currentTexture->Release();
                _currentTexture = RenderTexture(meshParamData->mesh, meshParamData->textureObject, true, _window);
           } else {
               std::cout << "No merges, nothing to do" << std::endl;
           }
        }

        if (primaryCharts.size() == 1) {
            if (ImGui::Button("Save current chart")) {
                PMesh pm;
                auto chart = meshParamData->GetChart(primaryCharts.begin()->first);
                Box2f b = chart->UVBox();

                auto f = [&pm, &b](typename Mesh::FacePointer fptr) {
                    auto f = tri::Allocator<PMesh>::AddFace(pm, fptr->P(0), fptr->P(1), fptr->P(2));
                    for (int i = 0; i < 3; ++i) {
                        f->WT(i) = fptr->WT(i);
                        f->WT(i).P() = (f->WT(i).P() - b.min) / std::max(b.DimX(), b.DimY());
                    }
                };

                std::for_each(chart->fpVec.begin(), chart->fpVec.end(), f);

                tri::Clean<PMesh>::RemoveDuplicateVertex(pm);
                tri::Allocator<PMesh>::CompactEveryVector(pm);

                tri::io::ExporterOBJ<PMesh>::Save(pm, "chart.obj", tri::io::Mask::IOM_WEDGTEXCOORD);
            }
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
                Mesh& m = meshParamData->mesh;
                if(SaveMesh(m, exportFileName, _currentTexture) == false) {
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

        ImGui::Text("Distortion metric");
        ImGui::RadioButton("Area Distortion", &distortionIndex, 0);
        ImGui::RadioButton("Edge Distortion", &distortionIndex, 1);
        ImGui::RadioButton("Angle Distortion", &distortionIndex, 2);
        if (distortionIndex != activeDistIndex) {
            meshParamData->MapDistortion(distortion[distortionIndex]);
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
                auto chart = meshParamData->GetChart(primaryCharts.begin()->first);
                ImGui::Text("Chart %lu (%lu faces, %lu adjacencies)", chart->id, chart->FN(), chart->NumAdj());
                ImGui::Text("Aggregate count: %d", chart->numMerges + 1);
                ImGui::Text("Area 3D: %.4f", chart->Area3D());
                ImGui::Text("Area UV: %.4f Border UV: %.6f", chart->AreaUV(), chart->BorderUV());
                ImGui::Text("Distortion range: %.4f , %.4f", chart->minMappedFaceValue, chart->maxMappedFaceValue);
            } else {
                // display info about the whole graph
                ImGui::Text("Filename: %s", fileName.c_str());
                ImGui::Text("Parameterization charts: %lu", meshParamData->Count());
                ImGui::Text("Area 3D: %.4f", meshParamData->Area3D());
                ImGui::Text("Area UV: %.4f Border UV: %.6f", meshParamData->AreaUV(), meshParamData->BorderUV());
                auto range = meshParamData->DistortionRange();
                ImGui::Text("Distortion range: %.4f , %.4f", range.first, range.second);
            }
            ImGui::End();
        }
    }

    // selection widget
    {
        if (selectedRegions.size() > 0) {
            ImGui::Begin("Selection widget", nullptr, 0);
            int i = 0;
            for (const auto& entry : selectedRegions) {
                RegionID id = entry.first;
                const SelectionBufferInfo& sbi = entry.second;
                ImVec4 bg = (primaryCharts.count(id) == 1) ? ImVec4(0.4f, 0.4f, 0.4f, 1.0f) : ImVec4(0.1f, 0.1f, 0.1f, 1.0f);
                if (ImGui::ImageButton((void *)(long unsigned)sbi.texicon, ImVec2(64, 64), ImVec2(0, 1), ImVec2(1, 0), -1, bg)) {
                    Select(id);
                }
                if ((++i % 4) != 0) ImGui::SameLine();
            }
            ImGui::End();
        }
    }


    if (updateTexcoord || updateColor) {
        glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.mesh);
        float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        for (const auto &f : meshParamData->mesh.face) {
            for (int i = 0; i < 3; ++i) {
                buffptr += 3;
                if (updateTexcoord) {
                    *buffptr++ = f.cWT(i).U();
                    *buffptr++ = f.cWT(i).V();
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
