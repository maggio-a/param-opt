#ifndef MESHVIEWER_H
#define MESHVIEWER_H

#include <iostream>

#include <QImage>

#include <vcg/space/ray3.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/intersection2.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "mesh_viewer.h"
#include "gl_util.h"

const char *vs_text_3D[] = {
    "#version 410 core                                                             \n"
    "                                                                              \n"
    "uniform mat4 modelViewMatrix;                                                 \n"
    "uniform mat4 projectionMatrix;                                                \n"
    "                                                                              \n"
    "in vec3 position;                                                             \n"
    "in vec2 texcoord;                                                             \n"
    "out vec2 uv;                                                                  \n"
    "                                                                              \n"
    "void main(void)                                                               \n"
    "{                                                                             \n"
    "    uv = texcoord;                                                            \n"
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
    "                                                                              \n"
    "void main(void)                                                               \n"
    "{                                                                             \n"
    "    uv = texcoord;                                                            \n"
    "    gl_Position = projectionMatrix * vec4(texcoord, 0.5f, 1.0f);              \n"
    "}                                                                             \n"
};

const char *fs_text_texture[] = {
    "#version 410 core                            \n"
    "                                             \n"
    "uniform sampler2D tex0; \n"
    "                                             \n"
    "in vec2 uv;                                  \n"
    "out vec4 color;                              \n"
    "                                             \n"
    "void main(void)                              \n"
    "{                                            \n"
    "    color = texture(tex0, uv);               \n"
    "}                                            \n"
};

const char *fs_text_selection[] = {
    "#version 410 core                            \n"
    "                                             \n"
    "in vec2 uv;                                  \n"
    "out vec4 color;                              \n"
    "                                             \n"
    "void main(void)                              \n"
    "{                                            \n"
    "    color = vec4(1.0, 0.0, 0.0, 0.5);        \n"
    "}                                            \n"
};

MeshViewer::MeshViewer(std::shared_ptr<MeshGraph> meshParamData_)
    : meshParamData{meshParamData_}, _textureCamera{}, _detailCamera{}
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

    float x = (_xpos - info.xSplit) / (info.height/2.0f);
    x += (_textureCamera.x - 0.5f);
    x /= _textureCamera.viewSize;
    float y = _ypos / (info.height/2.0f);
    y += (_textureCamera.y - 0.5f);
    y /= _textureCamera.viewSize;

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
        InitializeSelection(selectionID);
    } else {
        ClearSelection();
    }
}

void MeshViewer::ClearSelection()
{
    if (_selectionCount > 0) {
        glDeleteVertexArrays(1, &_perspectiveView.selection.vao);
        glDeleteVertexArrays(1, &_detailView.vao);
        glDeleteBuffers(1, &_vertexBuffers.selectedRegion);

        CheckGLError();

        _selectionCount = 0;
        _perspectiveView.selection.vao = 0;
        _detailView.vao = 0;
        _vertexBuffers.selectedRegion = 0;
    }
}

void MeshViewer::InitializeSelection(const RegionID id)
{
    auto region = meshParamData->GetChart(id);

    ClearSelection();
    _selectionCount = region->FN();

    glGenBuffers(1, &_vertexBuffers.selectedRegion);
    glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.selectedRegion);
    glBufferData(GL_ARRAY_BUFFER, region->FN()*15*sizeof(float), NULL, GL_STATIC_DRAW);
    float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto fptr : region->fpVec) {
        for (int i = 0; i < 3; ++i) {
            *buffptr++ = fptr->cV(i)->P().X();
            *buffptr++ = fptr->cV(i)->P().Y();
            *buffptr++ = fptr->cV(i)->P().Z();
            *buffptr++ = fptr->cWT(i).U();
            *buffptr++ = fptr->cWT(i).V();
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    glGenVertexArrays(1, &_perspectiveView.selection.vao);
    glBindVertexArray(_perspectiveView.selection.vao);

    _perspectiveView.selection.attributes.loc_position = glGetAttribLocation(_perspectiveView.selection.program, "position");
    glVertexAttribPointer(_perspectiveView.selection.attributes.loc_position, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
    glEnableVertexAttribArray(_perspectiveView.selection.attributes.loc_position);

    glGenVertexArrays(1, &_detailView.vao);
    glBindVertexArray(_detailView.vao);

    _detailView.attributes.loc_texcoord = glGetAttribLocation(_detailView.program, "texcoord");
    glVertexAttribPointer(_detailView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void *) (3*sizeof(float)));
    glEnableVertexAttribArray(_detailView.attributes.loc_texcoord);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    CheckGLError();

    SetupDetailView(id);
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
        InitializeSelection(selectionID);
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

void MeshViewer::MouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
    MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);

    if (mods & GLFW_MOD_SHIFT) {

        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) { // Center view on the point
            viewer->CenterPerspectiveViewFromMouse();
        }

    } else {

        // Drag mode // TODO: viewer->enable/disable dragmode
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
            if (viewer->InPerspectiveView()) viewer->_dragMode = DragMode::PERSPECTIVE;
            else if (viewer->InTextureView()) viewer->_dragMode = DragMode::TEXTURE;
            if (viewer->InDetailView()) viewer->_dragMode = DragMode::DETAIL;
        } else {
            viewer->_dragMode = DragMode::DISABLED;
        }

        // Pick mode, select region
        if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
            viewer->PickRegion();
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
    }

    viewer->_xpos = xpos;
    viewer->_ypos = ypos;
}

void MeshViewer::ScrollCallback(GLFWwindow* window, double xoffset, double yoffset)
{
    MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);
    if (viewer->InTextureView()) yoffset > 0.0f ? viewer->_textureCamera.ZoomIn() : viewer->_textureCamera.ZoomOut();
    if (viewer->InDetailView()) yoffset > 0.0f ? viewer->_detailCamera.ZoomIn() : viewer->_detailCamera.ZoomOut();
}

void MeshViewer::KeyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    MeshViewer *viewer = (MeshViewer *) glfwGetWindowUserPointer(window);

    // w toggles wireframe mode in detail view
    if (key == GLFW_KEY_W && action == GLFW_PRESS) {
        viewer->_detailView.wireframe = !viewer->_detailView.wireframe;
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
    glBufferData(GL_ARRAY_BUFFER, m.FN()*15*sizeof(float), NULL, GL_STATIC_DRAW);
    float *buffptr = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for(auto &f : m.face) {
        for (int i = 0; i < 3; ++i) {
            *buffptr++ = f.cV(i)->P().X();
            *buffptr++ = f.cV(i)->P().Y();
            *buffptr++ = f.cV(i)->P().Z();
            *buffptr++ = f.cWT(i).U();
            *buffptr++ = f.cWT(i).V();
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Load texture data
    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &_texture);
    LoadTexture2DFromQImage(*(meshParamData->textures[0]), _texture); // leaves the texture bound

    CheckGLError();

    // Initialize vertex array objects

    // Perspective view
    glGenVertexArrays(1, &_perspectiveView.vao);
    glBindVertexArray(_perspectiveView.vao);

    glBindBuffer(GL_ARRAY_BUFFER, _vertexBuffers.mesh);

    _perspectiveView.attributes.loc_position = glGetAttribLocation(_perspectiveView.program, "position");
    glVertexAttribPointer(_perspectiveView.attributes.loc_position, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
    glEnableVertexAttribArray(_perspectiveView.attributes.loc_position);

    CheckGLError();

    _perspectiveView.attributes.loc_texcoord = glGetAttribLocation(_perspectiveView.program, "texcoord");
    glVertexAttribPointer(_perspectiveView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void *) (3*sizeof(float)));
    glEnableVertexAttribArray(_perspectiveView.attributes.loc_texcoord);

    // Texture view
    glGenVertexArrays(1, &_textureView.vao);
    glBindVertexArray(_textureView.vao);

    _textureView.attributes.loc_texcoord = glGetAttribLocation(_textureView.program, "texcoord");
    glVertexAttribPointer(_textureView.attributes.loc_texcoord, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void *) (3*sizeof(float)));
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

    // view matrix is static
    mat4x4_look_at(_meshTransform.viewMatrix, _perspectiveCamera.eye, _perspectiveCamera.target, _perspectiveCamera.up);

    // projection will be updated at each draw
    mat4x4_identity(_meshTransform.projectionMatrix);

    _perspectiveView.uniforms.loc_modelView = glGetUniformLocation(_perspectiveView.program, "modelViewMatrix");
    _perspectiveView.uniforms.loc_projection = glGetUniformLocation(_perspectiveView.program, "projectionMatrix");

    _perspectiveView.selection.uniforms.loc_modelView = glGetUniformLocation(_perspectiveView.selection.program, "modelViewMatrix");
    _perspectiveView.selection.uniforms.loc_projection = glGetUniformLocation(_perspectiveView.selection.program, "projectionMatrix");

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
        mat4x4_identity(rotation);
        mat4x4_invert(invO, _meshTransform.orientationMatrix);
        mat4x4_mul_vec4(orientedX, invO, x);
        mat4x4_mul_vec4(orientedY, invO, y);
        vec4_norm(orientedX, orientedX);
        mat4x4_rotate(rotation, rotation, orientedX[0], orientedX[1], orientedX[2], _dragY * M_PI / info.height);
        mat4x4_mul(_meshTransform.orientationMatrix, _meshTransform.orientationMatrix, rotation);
        mat4x4_rotate(_meshTransform.orientationMatrix, _meshTransform.orientationMatrix,
                      orientedY[0], orientedY[1], orientedY[2],
                      _dragX * M_PI / info.width);
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
    glBindVertexArray(_perspectiveView.vao);
    glUseProgram(_perspectiveView.program);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _texture);

    mat4x4 model;
    mat4x4 modelView;

    // Build model matrix: translate, orient and scale
    mat4x4_dup(model, _meshTransform.positionMatrix);
    mat4x4_mul(model, _meshTransform.orientationMatrix, model);
    mat4x4_mul(model, _meshTransform.scaleMatrix, model);

    mat4x4_mul(modelView, _meshTransform.viewMatrix, model);

    mat4x4_perspective(_meshTransform.projectionMatrix, 60.0f * M_PI / 180.0f, info.perspectiveViewAspect, 0.1f, 2000.0f);

    glUniformMatrix4fv(_perspectiveView.uniforms.loc_modelView, 1, GL_FALSE, (const GLfloat *) modelView);
    glUniformMatrix4fv(_perspectiveView.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *)_meshTransform.projectionMatrix);

    glViewport(0, 0, info.xSplit, info.height);
    glScissor(0, 0, info.xSplit, info.height);

    glDrawBuffer(GL_BACK);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 0, meshParamData->mesh.FN() * 3);
    glBindVertexArray(0);

    if (_selectionCount > 0) {
        glBindVertexArray(_perspectiveView.selection.vao);
        glUseProgram(_perspectiveView.selection.program);
        glUniformMatrix4fv(_perspectiveView.selection.uniforms.loc_modelView, 1, GL_FALSE, (const GLfloat *) modelView);
        glUniformMatrix4fv(_perspectiveView.selection.uniforms.loc_projection, 1, GL_FALSE, (const GLfloat *)_meshTransform.projectionMatrix);
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDrawArrays(GL_TRIANGLES, 0, _selectionCount * 3);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
    }

    CheckGLError();
}

void MeshViewer::DrawTextureView()
{
    glBindVertexArray(_textureView.vao);
    glUseProgram(_textureView.program);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _texture);

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
    glBindTexture(GL_TEXTURE_2D, _texture);

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

    glClearColor(0.2f, 0.7f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    if (_detailView.wireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }

    glDrawArrays(GL_TRIANGLES, 0, _selectionCount * 3);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    CheckGLError();
}

void MeshViewer::DrawViews()
{
    Draw3DView();
    DrawTextureView();
    if (_selectionCount > 0) DrawDetailView();
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
    glfwSetKeyCallback(_window, MeshViewer::KeyCallback);
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

    GLint loc_tex0;
    _perspectiveView.program = CompileShaders(vs_text_3D, fs_text_texture);
    _perspectiveView.selection.program = CompileShaders(vs_text_3D, fs_text_selection);
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

    while (!glfwWindowShouldClose(_window)) {

        glScissor(0, 0, info.width, info.height);
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        UpdateTransforms();
        DrawViews();
        glfwSwapBuffers(_window);

        glfwPollEvents();
    }

    glfwDestroyWindow(_window);
    glfwTerminate();

}

#endif // MESHVIEWER_H
