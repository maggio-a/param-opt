#include <iostream>
#include <cstdio>

#include <QImage>

#include <vcg/space/ray3.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/intersection2.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <wrap/gui/trackball.h>

#include <imgui.h>
#include <imgui_glfw_gl3/imgui_impl_glfw_gl3.h>

#include "mesh.h"
#include "mesh_viewer.h"
#include "gl_util.h"
#include "optimizer.h"
#include "texture_rendering.h"

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
    "                                                                              \n"
    "in vec2 texcoord;                                                             \n"
    "out vec2 uv;                                                                  \n"
    "out vec4 dcol;                                                                \n"
    "                                                                              \n"
    "void main(void)                                                               \n"
    "{                                                                             \n"
    "    uv = texcoord;                                                            \n"
    "    dcol = vec4(1.0f, 1.0f, 1.0f, 1.0f); // unused but required by the fs     \n"
    "    gl_Position = projectionMatrix * vec4(texcoord, 0.5f, 1.0f);              \n"
    "}                                                                             \n"
};

const char *fs_text_texture[] = {
    "#version 410 core                                              \n"
    "                                                               \n"
    "const int COLOR_SRC_TEXTURE    = 0;                            \n"
    "const int COLOR_SRC_DISTORTION = 1;                            \n"
    "const int COLOR_SRC_CHECKBOARD = 2;                            \n"
    "                                                               \n"
    "uniform sampler2D tex0;                                        \n"
    "uniform float resolution = 64.0f;                              \n"
    "uniform int colorSource = COLOR_SRC_TEXTURE;                   \n"
    "uniform vec4 weight = vec4(1.0f, 1.0f, 1.0f, 1.0f);            \n"
    "                                                               \n"
    "in vec2 uv;                                                    \n"
    "in vec4 dcol;                                                  \n"
    "out vec4 color;                                                \n"
    "                                                               \n"
    "void main(void)                                                \n"
    "{                                                              \n"
    "    if (colorSource == COLOR_SRC_TEXTURE) {                    \n"
    "        color = texture(tex0, uv)*weight;                      \n"
    "    }                                                          \n"
    "    else if (colorSource == COLOR_SRC_DISTORTION) {            \n"
    "        color = dcol*weight;                                   \n"
    "    }                                                          \n"
    "    else if (colorSource == COLOR_SRC_CHECKBOARD) {            \n"
    "        float u_mod = mod(floor(uv.s * resolution), 2.0f);     \n"
    "        float v_mod = mod(floor(uv.t * resolution), 2.0f);     \n"
    "        float c = abs(u_mod - v_mod);                          \n"
    "        color = vec4(c, c, c, 1.0f);                           \n"
    "    } else {                                                   \n"
    "        color = vec4(0.0f, 1.0f, 0.0f, 1.0f);                  \n"
    "    }                                                          \n"
    "}                                                              \n"
};

// GLFW callbacks

void MeshViewer::MouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
    if (action == GLFW_PRESS && button >= 0 && button < 3)
        g_MouseJustPressed[button] = true;

    ImGuiIO& io = ImGui::GetIO();

    // forward input to app logic only if gui is not involved or it is a release event
    if (io.WantCaptureMouse == false || action == GLFW_RELEASE) {
        MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);
        if (mods & GLFW_MOD_SHIFT) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) { // Center view on the point
                viewer->CenterPerspectiveViewFromMouse();
            }
        } else {
            // Drag mode // TODO: viewer->enable/disable dragmode
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                if (viewer->InPerspectiveView()) {
                    viewer->_dragMode = DragMode::PERSPECTIVE;
                    trackball.MouseDown(viewer->_xpos, viewer->info.height - viewer->_ypos, vcg::Trackball::BUTTON_LEFT);
                }
                else if (viewer->InTextureView()) viewer->_dragMode = DragMode::TEXTURE;
                if (viewer->InDetailView()) viewer->_dragMode = DragMode::DETAIL;
            } else {
                if (viewer->_dragMode == DragMode::PERSPECTIVE)
                    trackball.MouseUp(viewer->_xpos, viewer->info.height - viewer->_ypos, vcg::Trackball::BUTTON_LEFT);
                viewer->_dragMode = DragMode::DISABLED;
            }
            // Pick mode, select region
            if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
                viewer->PickRegion();
            }
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

    if (viewer->_dragMode != DragMode::DISABLED) {
        viewer->_dragX += (int) (xpos - viewer->_xpos);
        viewer->_dragY += (int) (ypos - viewer->_ypos);

        if (viewer->_dragMode == DragMode::PERSPECTIVE) {
            trackball.center = {0, 0, 0};
            trackball.radius = 1;
            Matrix44f proj;
            memcpy(&proj, &viewer->_meshTransform.projectionMatrix, 16*sizeof(float));
            proj.transposeInPlace();
            Matrix44f mv;
            mv.SetIdentity();
            int viewport[4] = {0, 0, viewer->info.xSplit, viewer->info.height};

            trackball.camera.SetView(proj.V(), mv.V(), viewport);
            trackball.MouseMove((int) xpos, viewer->info.height - ((int) ypos));
        }
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
            float factor = (yoffset > 0.0f) ? 0.9f : 1.1f;
            viewer->_perspectiveCamera.eye[2] *= factor;
            viewer->_perspectiveCamera.near *= factor;
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

void MeshViewer::FramebufferSizeCallback(GLFWwindow *window, int width, int height)
{
    MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);
    if (width <= 0) width = 1;
    if (height <= 0) height = 1;
    viewer->info.width = width;
    viewer->info.height = height;
    viewer->info.xSplit = width > (height/2) ? (width - (height/2)) : 1;
    viewer->info.perspectiveViewAspect = viewer->info.xSplit / (float) height;
}


// Member functions


MeshViewer::MeshViewer(std::shared_ptr<MeshGraph> meshParamData_, std::size_t minRegionSize_, const std::string& fileName_)
    : meshParamData{meshParamData_}, gm{std::make_shared<GraphManager>(meshParamData_)},
      minRegionSize{minRegionSize_}, _currentTexture{meshParamData_->textureObject}, fileName{fileName_}, _textureCamera{}, _detailCamera{}
{
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
    if (selectionVector.size() > 0) {
        glDeleteVertexArrays(1, &_perspectiveView.selection.vao);
        glDeleteVertexArrays(1, &_detailView.vao);
        glDeleteBuffers(1, &_vertexBuffers.selection);

        CheckGLError();

        _perspectiveView.selection.vao = 0;
        _detailView.vao = 0;
        _vertexBuffers.selection = 0;

        selectionType = SelectionType::None;
        selectionVector.clear();
    }
}

void MeshViewer::Select(const RegionID id)
{
    auto region = meshParamData->GetChart(id);
    std::vector<std::pair<RegionID,vcg::Color4f>> vsel;
    vsel.reserve(1 + region->NumAdj());

    vsel.emplace_back(std::make_pair(region->id, vcg::Color4f{1.0f, 1.0f, 1.0f, 1.0f}));

    for (auto adj : region->adj) {
        float c = adj->id / (float) meshParamData->charts.size();
        vsel.emplace_back(std::make_pair(adj->id, vcg::Color4f{c*0.6f, 1.0f*0.6f, (1.0f-c)*0.6f, 1.0f}));
    }

    InitializeSelection(vsel);
    selectionType = SelectionType::Chart;

    std::cout << "Selected region " << id << " (uv area = " << region->AreaUV() << ")" << std::endl;
}

void MeshViewer::Select(const GraphManager::Edge& e)
{
    InitializeSelection({
        {e.a->id, vcg::Color4f{1.0f, 0.1f, 0.1f, 1.0f}},
        {e.b->id, vcg::Color4f{0.1f, 0.1f, 1.0f, 1.0f}}
    });
    selectionType = SelectionType::Edge;

    std::cout << "Selected edge (" << e.a->id << " , " << e.b->id << ")" << std::endl;
}

void MeshViewer::InitializeSelection(const std::vector<std::pair<RegionID,vcg::Color4f>>& vsel)
{
    ClearSelection();
    std::size_t selectionCount = 0;
    for (const auto& p : vsel) {
        selectionCount += meshParamData->GetChart(p.first)->FN();
    }

    selectionVector.reserve(vsel.size());
    GLint first = 0;

    glGenBuffers(1, &_vertexBuffers.selection);
    glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.selection);
    glBufferData(GL_ARRAY_BUFFER, selectionCount*15*sizeof(float), NULL, GL_STATIC_DRAW);
    float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (const auto& p : vsel) {
        const auto& region = meshParamData->GetChart(p.first);
        selectionVector.emplace_back(SelectedRegionInfo{region, first, static_cast<GLsizei>(region->FN() * 3), p.second});
        for (auto &fptr : region->fpVec) {
            for (int i = 0; i < 3; ++i) {
                *buffptr++ = fptr->cV(i)->P().X();
                *buffptr++ = fptr->cV(i)->P().Y();
                *buffptr++ = fptr->cV(i)->P().Z();
                *buffptr++ = fptr->cWT(i).U();
                *buffptr++ = fptr->cWT(i).V();
            }
        }
        first += region->FN() * 3;
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    glGenVertexArrays(1, &_perspectiveView.selection.vao);
    glBindVertexArray(_perspectiveView.selection.vao);

    _perspectiveView.selection.attributes.loc_position = glGetAttribLocation(_perspectiveView.selection.program, "position");
    glVertexAttribPointer(_perspectiveView.selection.attributes.loc_position, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
    glEnableVertexAttribArray(_perspectiveView.selection.attributes.loc_position);
    _perspectiveView.selection.attributes.loc_texcoord = glGetAttribLocation(_perspectiveView.selection.program, "texcoord");
    glVertexAttribPointer(_perspectiveView.selection.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float),
                          (const GLvoid *) (3*sizeof(float)));
    glEnableVertexAttribArray(_perspectiveView.selection.attributes.loc_texcoord);

    glGenVertexArrays(1, &_detailView.vao);
    glBindVertexArray(_detailView.vao);

    _detailView.attributes.loc_texcoord = glGetAttribLocation(_detailView.program, "texcoord");
    glVertexAttribPointer(_detailView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (const GLvoid *) (3*sizeof(float)));
    glEnableVertexAttribArray(_detailView.attributes.loc_texcoord);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    CheckGLError();

    SetupDetailView(vsel[0].first);
}

bool MeshViewer::IntersectionMouseRayModel(Mesh::ConstFacePointer *fp, float &u, float &v)
{
    using vcg::Point3f;
    using vcg::Ray3f;

    // Compute ray from mouse position

    mat4x4 model, invModel;
    mat4x4_dup(model, _meshTransform.positionMatrix);
    mat4x4_mul(model, _meshTransform.orientationMatrix, model);
    mat4x4_mul(model, _meshTransform.scaleMatrix, model);
    mat4x4_invert(invModel, model);

    mat4x4 invProj, invView, ndcToWorld;
    mat4x4_invert(invProj, _meshTransform.projectionMatrix);
    mat4x4_invert(invView, _meshTransform.viewMatrix);
    mat4x4_mul(ndcToWorld, invView, invProj);

    float x = (_xpos / (float) info.xSplit) * 2.0f - 1.0f;
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
            mat4x4_translate(_meshTransform.positionMatrix, -centerPoint.X(), -centerPoint.Y(), -centerPoint.Z());
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

    _textureView.attributes.loc_texcoord = glGetAttribLocation(_textureView.program, "texcoord");
    glVertexAttribPointer(_textureView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void *) (3*sizeof(float)));
    glEnableVertexAttribArray(_textureView.attributes.loc_texcoord);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    CheckGLError();
}

void MeshViewer::SetupViews()
{
    const Mesh& m = meshParamData->mesh;

    // Setup perspective view

    vcg::Box3f bbox = m.bbox;
    if (bbox.IsNull()) {
        for (const auto& v : m.vert) {
            bbox.Add(v.cP());
        }
    }

    // Define model transform: the mesh will be translated to the origin of the world space,
    // rotated according to accumulated rotations and then scaled

    mat4x4_translate(_meshTransform.positionMatrix, -bbox.Center().X(), -bbox.Center().Y(), -bbox.Center().Z());
    mat4x4_identity(_meshTransform.orientationMatrix); // Initially the mesh is not rotated
    float scale = _perspectiveCamera.eye[2] / bbox.Diag();
    mat4x4_identity(_meshTransform.scaleMatrix);
    mat4x4_scale_aniso(_meshTransform.scaleMatrix, _meshTransform.scaleMatrix, scale, scale, scale);

    // view and projection matrices will be updated at each draw
    mat4x4_identity(_meshTransform.viewMatrix);
    mat4x4_identity(_meshTransform.projectionMatrix);

    _perspectiveView.uniforms.loc_modelView = glGetUniformLocation(_perspectiveView.program, "modelViewMatrix");
    _perspectiveView.uniforms.loc_projection = glGetUniformLocation(_perspectiveView.program, "projectionMatrix");
    _perspectiveView.uniforms.loc_colorSource = glGetUniformLocation(_perspectiveView.program, "colorSource");
    _perspectiveView.uniforms.loc_weight = glGetUniformLocation(_perspectiveView.program, "weight");

    _perspectiveView.selection.uniforms.loc_modelView = glGetUniformLocation(_perspectiveView.selection.program, "modelViewMatrix");
    _perspectiveView.selection.uniforms.loc_projection = glGetUniformLocation(_perspectiveView.selection.program, "projectionMatrix");
    _perspectiveView.selection.uniforms.loc_weight = glGetUniformLocation(_perspectiveView.selection.program, "weight");

    // Setup texture view

    _textureCamera.Reset();
    _textureView.uniforms.loc_projection = glGetUniformLocation(_textureView.program, "projectionMatrix");
}

void MeshViewer::SetupDetailView(const RegionID id)
{
    auto region = meshParamData->GetChart(id);

    _detailCamera.Reset();
    _detailView.uniforms.loc_projection = glGetUniformLocation(_detailView.program, "projectionMatrix");

    Box2f bbox;
    for (auto fptr : region->fpVec) {
        for (int i = 0; i < fptr->VN(); ++i) {
            bbox.Add(fptr->cWT(i).P());
        }
    }

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
    mat4x4 invO, rotation;
    vec4 x = {1.0f, 0.0f, 0.0f, 0.0f}, y = {0.0f, 1.0f, 0.0f, 0.0f}, orientedX, orientedY;


    switch (_dragMode) {
    case DragMode::PERSPECTIVE:
        {/*
            trackball.center = {0, 0, 0};
            trackball.radius = 1;
            Matrix44f proj;
            memcpy(&proj, &_meshTransform.projectionMatrix, 16*sizeof(float));
            proj.transposeInPlace();
            Matrix44f mv;
            mv.Identity();
            int viewport[4] = {0, 0, info.xSplit, info.height};


            trackball.camera.SetView(proj.V(), mv.V(), viewport);
            trackball.MouseMove((int) _xpos, (int) _ypos);*/

            Matrix44f rot = trackball.Matrix().transpose();
            //rot.print();
            memcpy(&_meshTransform.orientationMatrix, rot.V(), 16*sizeof(float));
        }
        /*mat4x4_identity(rotation);
        mat4x4_invert(invO, _meshTransform.orientationMatrix);
        mat4x4_mul_vec4(orientedX, invO, x);
        mat4x4_mul_vec4(orientedY, invO, y);
        vec4_norm(orientedX, orientedX);
        mat4x4_rotate(rotation, rotation, orientedX[0], orientedX[1], orientedX[2], _dragY * M_PI / info.height);
        mat4x4_mul(_meshTransform.orientationMatrix, _meshTransform.orientationMatrix, rotation);
        mat4x4_rotate(_meshTransform.orientationMatrix, _meshTransform.orientationMatrix,
                      orientedY[0], orientedY[1], orientedY[2],
                      _dragX * M_PI / info.width);*/
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
    mat4x4_dup(model, _meshTransform.positionMatrix);
    mat4x4_mul(model, _meshTransform.orientationMatrix, model);
    mat4x4_mul(model, _meshTransform.scaleMatrix, model);

    mat4x4_look_at(_meshTransform.viewMatrix, _perspectiveCamera.eye, _perspectiveCamera.target, _perspectiveCamera.up);

    mat4x4_mul(modelView, _meshTransform.viewMatrix, model);

    mat4x4_perspective(_meshTransform.projectionMatrix, 60.0f * M_PI / 180.0f, info.perspectiveViewAspect,
                       _perspectiveCamera.near, _perspectiveCamera.far);

    glViewport(0, 0, info.xSplit, info.height);
    glScissor(0, 0, info.xSplit, info.height);

    glDrawBuffer(GL_BACK);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0);
    _currentTexture->Bind();

    glEnable(GL_DEPTH_TEST);

    if (selectionVector.size() > 0) {
        glBindVertexArray(_perspectiveView.selection.vao);
        glUseProgram(_perspectiveView.selection.program);

        glUniformMatrix4fv(_perspectiveView.selection.uniforms.loc_modelView, 1, GL_FALSE, (const GLfloat *) modelView);
        glUniformMatrix4fv(_perspectiveView.selection.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *)_meshTransform.projectionMatrix);

        glEnable(GL_STENCIL_TEST);
        glStencilFuncSeparate(GL_FRONT, GL_ALWAYS, 1, 0xFF);
        glStencilOpSeparate(GL_FRONT, GL_KEEP, GL_KEEP, GL_REPLACE);

        glStencilFuncSeparate(GL_BACK, GL_ALWAYS, 0, 0xFF);
        glStencilOpSeparate(GL_BACK, GL_KEEP, GL_KEEP, GL_REPLACE);

        for (const SelectedRegionInfo& sri : selectionVector) {
            glUniform4fv(_perspectiveView.selection.uniforms.loc_weight, 1, sri.color.V());
            glDrawArrays(GL_TRIANGLES, sri.first, sri.count);

        }

        // Now draw transparent layer
        glBindVertexArray(_perspectiveView.vao);
        glUseProgram(_perspectiveView.program);

        glUniformMatrix4fv(_perspectiveView.uniforms.loc_modelView, 1, GL_FALSE, (const GLfloat *) modelView);
        glUniformMatrix4fv(_perspectiveView.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *)_meshTransform.projectionMatrix);

        glUniform1i(_perspectiveView.uniforms.loc_colorSource, _perspectiveView.colorSource);
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

        glUniform1i(_perspectiveView.uniforms.loc_colorSource, _perspectiveView.colorSource);
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

    glViewport(info.xSplit, info.height/2, info.height/2, info.height/2);
    glScissor(info.xSplit, info.height/2, info.height/2, info.height/2);

    glDrawBuffer(GL_BACK);
    glDisable(GL_DEPTH_TEST);

    glClearColor(0.2f, 0.7f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 0, meshParamData->mesh.FN() * 3);
    glBindVertexArray(0);

    CheckGLError();
}

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

    glViewport(info.xSplit, 0, info.height/2, info.height/2);
    glScissor(info.xSplit, 0, info.height/2, info.height/2);

    glDrawBuffer(GL_BACK);
    glDisable(GL_DEPTH_TEST);

    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    if (_detailView.wireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }

    glDrawArrays(GL_TRIANGLES, selectionVector[0].first, selectionVector[0].count);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    CheckGLError();
}

void MeshViewer::DrawViews()
{
    Draw3DView();
    DrawTextureView();
    if (selectionVector.size() > 0) DrawDetailView();
}


void MeshViewer::Run()
{
    trackball.SetIdentity();

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
    _detailView.program = CompileShaders(vs_text_texture, fs_text_texture);

    loc_tex0 = glGetUniformLocation(_perspectiveView.program, "tex0");
    glUseProgram(_perspectiveView.program);
    glUniform1i(loc_tex0, 0);
    loc_tex0 = glGetUniformLocation(_textureView.program, "tex0");
    glUseProgram(_textureView.program);
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

        glfwPollEvents();
        ImGui_ImplGlfwGL3_NewFrame();

        ManageImGuiState();

        glScissor(0, 0, info.width, info.height);
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

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
    static bool showInfoArea = false;

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

        if (ImGui::Button("Toggle wireframe")) {
            _detailView.wireframe ^= true;
        }

        if (ImGui::Button("Highlight next merge")) {
            if (gmHasNextEdge()) {
                const auto& next = gmPeekNextEdge();
                Select(next.first);
                std::cout << "Next edge has weight " << next.second << std::endl;
            }
            else {
                std::cout << "No next edge available" << std::endl;
            }
        }


        if (ImGui::Button("Close islands")) {
            ClearSelection();
            gmClose();
        }

        if (ImGui::Button("Perform next merge")) {
            if (gmHasNextEdge()) {
                auto next = gmPeekNextEdge();
                if (next.first.a->FN() > minRegionSize && next.first.b->FN() > minRegionSize) {
                    std::cout << "Next edge cannot be collapsed (charts are too large)" << std::endl;
                }
                else {
                    std::cout << "Merging " << next.first.a->id << " and " << next.first.b->id << " (weight=" << next.second << ")" << std::endl;
                    gmRemoveNextEdge();
                    auto chart = gmCollapse(next.first);
                    Select(chart->id);
                }
            }
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

        ImGui::RadioButton("Area Distortion", &distortionIndex, 0);
        ImGui::RadioButton("Edge Distortion", &distortionIndex, 1);
        ImGui::RadioButton("Angle Distortion", &distortionIndex, 2);
        if (ImGui::Button("Map distortion to faces")) {
            meshParamData->MapDistortion(distortion[distortionIndex]);
            activeDistIndex = distortionIndex;
            updateColor = true;
        }

        ImGui::Separator();

        ImGui::Combo("Color source", &_perspectiveView.colorSource, "Texture\0Distortion\0Checkboard\0\0");

        ImGui::End();
    }

    // info area
    {
        if (showInfoArea) {
            ImGui::Begin("Info area", &showInfoArea);
            if (selectionType == SelectionType::Chart) {
                // display info about the chart
                auto chart = selectionVector[0].chart;
                ImGui::Text("Chart %lu (%lu faces)", chart->id, chart->FN());
                ImGui::Text("Area 3D: %.4f", chart->Area3D());
                ImGui::Text("Area UV: %.4f", chart->AreaUV());
                ImGui::Text("Aggregate count: %d", chart->numMerges + 1);
                ImGui::Text("Distortion range: %.4f , %.4f", chart->minMappedFaceValue, chart->maxMappedFaceValue);
            } else {
                // display info about the whole graph
                ImGui::Text("Filename: %s", fileName.c_str());
                ImGui::Text("Parameterization charts: %lu", meshParamData->Count());
                auto range = meshParamData->DistortionRange();
                ImGui::Text("Distortion range: %.4f , %.4f", range.first, range.second);
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
