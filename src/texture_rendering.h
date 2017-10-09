#ifndef TEXTURE_RENDERING_H
#define TEXTURE_RENDERING_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>

#include <QImage>

#include "mesh.h"
#include "pushpull.h"

#include "gl_util.h"

using namespace std;

static const char *vs_text[] = {
    "#version 420 core                                           \n"
    "                                                            \n"
    "in vec2 position;                                           \n"
    "in vec2 texcoord;                                           \n"
    "out vec2 uv;                                                \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    uv = texcoord;                                          \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);               \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                        \n"
    "}                                                           \n"
};

static const char *fs_text[] = {
    "#version 420 core                                           \n"
    "                                                            \n"
    "layout (binding = 0) uniform sampler2D img0;                \n"
    "in vec2 uv;                                                 \n"
    "                                                            \n"
    "out vec4 color;                                             \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    color = texture(img0, uv);                              \n"
    "}                                                           \n"
};

static std::shared_ptr<QImage> RenderTexture(Mesh &m, std::vector<std::shared_ptr<QImage>> imgVec, bool filter)
{
    assert(imgVec.size() == 1); // FIXME multiple texture images support
    QImage &img = *imgVec[0];

    if (!glfwInit())
    {
        cout << "Failed to initialize glfw" << endl;
        exit(-1);
    }

    glfwSetErrorCallback(ErrorCallback);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    GLFWwindow *window = glfwCreateWindow(512, 512, "Window", NULL, NULL);
    if (!window)
    {
        cout << "Failed to create window or context" << endl;
    }

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err)
    {
        cout << "glew init error " << glewGetErrorString(err) << endl;
    }

    glGetError(); // suppress possible error on glew init

    glfwSwapInterval(1);

    // OpenGL setup

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLint status;
    char str[1024] = {0};

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, vs_text, NULL);
    glCompileShader(vs);

    glGetShaderInfoLog(vs, 1024, NULL, str);
    cout << str;
    memset(str, 0, 1024);

    glGetShaderiv(vs, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE)
    {
        cout << "vs compilation failed" << endl;
    }

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, fs_text, NULL);
    glCompileShader(fs);

    glGetShaderInfoLog(fs, 1024, NULL, str);
    cout << str << endl;
    memset(str, 0, 1024);

    glGetShaderiv(fs, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE)
    {
        cout << "fs compilation failed" << endl;
    }

    GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    glValidateProgram(program);

    glGetProgramInfoLog(program, 1024, NULL, &str[0]);
    cout << str << endl;

    glGetProgramiv(program, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        cout << "link failed" << endl;
    }

    glDeleteShader(vs);
    glDeleteShader(fs);

    CheckGLError();

    glUseProgram(program);

    CheckGLError();

    // Allocate buffers for vertex data, map it and copy tex coords

    GLuint buffer[2];
    glGenBuffers(2, buffer);

    // load vertex coordinates
    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER, m.FN()*6*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            *p++ = f.cWT(i).U();
            *p++ = f.cWT(i).V();
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint pos_location = glGetAttribLocation(program, "position");
    glVertexAttribPointer(pos_location, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(pos_location);

    // load texture coordinates
    glBindBuffer(GL_ARRAY_BUFFER, buffer[1]);
    glBufferData(GL_ARRAY_BUFFER, m.FN()*6*sizeof(float), NULL, GL_STATIC_DRAW);
    Mesh::PerFaceAttributeHandle<WedgeTexCoordStorage> WTCSh
            = tri::Allocator<Mesh>::GetPerFaceAttribute<WedgeTexCoordStorage>(m, "WedgeTexCoordStorage");
    p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            *p++ = WTCSh[&f].tc[i].U();
            *p++ = WTCSh[&f].tc[i].V();
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);
    GLint tc_location = glGetAttribLocation(program, "texcoord");
    glVertexAttribPointer(tc_location, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(tc_location);

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0); // done, unbind

    // Setup FBO

    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glViewport(0, 0, img.width(), img.height());

    GLuint renderTarget;
    glGenTextures(1, &renderTarget);
    glBindTexture(GL_TEXTURE_2D, renderTarget);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, img.width(), img.height());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, renderTarget, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    // Load texture image

    GLuint texture;
    glGenTextures(1, &texture);
    glActiveTexture(GL_TEXTURE0);
    LoadTexture2DFromQImage(img, texture); // mirror qimage to match opengl convention

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    QImage textureImage(img.width(), img.height(), QImage::Format_ARGB32);

    glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 0, m.FN()*3);

    CheckGLError();

    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, textureImage.width(), textureImage.height(), GL_BGRA, GL_UNSIGNED_BYTE, textureImage.bits());

    glfwPollEvents();


    // cleanup TODO

    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);
    glfwTerminate();

    if (filter) vcg::PullPush(textureImage, qRgba(0, 255, 0, 255));

    return std::make_shared<QImage>(textureImage.mirrored());
}


#endif // TEXTURE_RENDERING_H

