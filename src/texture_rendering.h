/*
 * References for the bicubic interpolated texture lookup:
 *  - GPU gems 2 ch 20 (Sigg and Hadwiger 2005)
 *  - Efficient GPU-Based Texture Interpolation using Uniform B-Splines  (Ruijters et al. 2009)
 * */

#ifndef TEXTURE_RENDERING_H
#define TEXTURE_RENDERING_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>

#include <QImage>

#include "mesh.h"
#include "pushpull.h"
#include "uv.h"

#include "gl_utils.h"

using namespace std;

static const char *vs_text[] = {
    "#version 410 core                                           \n"
    "                                                            \n"
    "in vec2 position;                                           \n"
    "in vec2 texcoord;                                           \n"
    "out vec2 uv;                                                \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    uv = texcoord;                                          \n"
    "    if (uv.s < 0) uv = vec2(0.0, 0.0);                      \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);               \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                        \n"
    "}                                                           \n"
};

static const char *fs_text[] = {
    "#version 410 core                                                     \n"
    "                                                                      \n"
    "uniform sampler2D img0;                                               \n"
    "                                                                      \n"
    "uniform vec2 texture_size;                                            \n"
    "uniform int use_cubic_interpolation;                                  \n"
    "                                                                      \n"
    "in vec2 uv;                                                           \n"
    "                                                                      \n"
    "out vec4 color;                                                       \n"
    "                                                                      \n"
    "void main(void)                                                       \n"
    "{                                                                     \n"
    "    if (use_cubic_interpolation == 0)                                 \n"
    "        color = texture2D(img0, uv);                                  \n"
    "    else {                                                            \n"
    "        vec2 coord = uv * texture_size - vec2(0.5, 0.5);              \n"
    "        vec2 idx = floor(coord);                                      \n"
    "        vec2 fraction = coord - idx;                                  \n"
    "        vec2 one_frac = vec2(1.0, 1.0) - fraction;                    \n"
    "        vec2 one_frac2 = one_frac * one_frac;                         \n"
    "        vec2 fraction2 = fraction * fraction;                         \n"
    "        vec2 w0 = (1.0/6.0) * one_frac2 * one_frac;                   \n"
    "        vec2 w1 = (2.0/3.0) - 0.5 * fraction2 * (2.0 - fraction);     \n"
    "        vec2 w2 = (2.0/3.0) - 0.5 * one_frac2 * (2.0 - one_frac);     \n"
    "        vec2 w3 = (1.0/6.0) * fraction2 * fraction;                   \n"
    "        vec2 g0 = w0 + w1;                                            \n"
    "        vec2 g1 = w2 + w3;                                            \n"
    "        vec2 h0 = (w1 / g0) - 0.5 + idx;                              \n"
    "        vec2 h1 = (w3 / g1) + 1.5 + idx;                              \n"
    "        vec4 tex00 = texture2D(img0, vec2(h0.x, h0.y) / texture_size);\n"
    "        vec4 tex10 = texture2D(img0, vec2(h1.x, h0.y) / texture_size);\n"
    "        vec4 tex01 = texture2D(img0, vec2(h0.x, h1.y) / texture_size);\n"
    "        vec4 tex11 = texture2D(img0, vec2(h1.x, h1.y) / texture_size);\n"
    "        tex00 = mix(tex00, tex01, g1.y);                              \n"
    "        tex10 = mix(tex10, tex11, g1.y);                              \n"
    "        color = mix(tex00, tex10, g1.x);                              \n"
    "    }                                                                 \n"
    "}                                                                     \n"
};


// if parent == nullptr this function creates an exclusive context, otherwise set
// up a shared context so textures arent duplicated
static TextureObjectHandle RenderTexture(Mesh::FaceIterator fbegin, Mesh::FaceIterator fend,
                                         Mesh &m, TextureObjectHandle textureObject,
                                         bool filter, bool useCubicInterpolation, GLFWwindow *parentWindow);

static TextureObjectHandle RenderTexture(Mesh &m, TextureObjectHandle textureObject, bool filter, bool useCubicInterpolation, GLFWwindow *parentWindow)
{
    return RenderTexture(m.face.begin(), m.face.end(), m, textureObject, filter, useCubicInterpolation, parentWindow);
}

static TextureObjectHandle RenderTexture(Mesh::FaceIterator fbegin, Mesh::FaceIterator fend,
                                         Mesh &m, TextureObjectHandle textureObject,
                                         bool filter, bool useCubicInterpolation, GLFWwindow *parentWindow)
{
    assert(textureObject->ArraySize() == 1); // TODO multiple texture images support

    bool sharedContext = (parentWindow != nullptr);

    QImage &img = *(textureObject->imgVec[0]);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    GLFWwindow *window = glfwCreateWindow(512, 512, "Window", NULL, parentWindow);
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

    // Allocate vertex data

    // Note that if the viewer is running texture coords are already in a gpu buffer so I could re-use those
    auto WTCSh = tri::Allocator<Mesh>::FindPerFaceAttribute<TexCoordStorage>(m, "WedgeTexCoordStorage");
    assert(tri::Allocator<Mesh>::IsValidHandle<TexCoordStorage>(m, WTCSh));

    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, m.FN()*12*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto it = fbegin; it != fend; ++it) {
        for (int i = 0; i < 3; ++i) {
            *p++ = it->cWT(i).U();
            *p++ = it->cWT(i).V();
            *p++ = WTCSh[&*it].tc[i].U();
            *p++ = WTCSh[&*it].tc[i].V();
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint pos_location = glGetAttribLocation(program, "position");
    glVertexAttribPointer(pos_location, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), 0);
    glEnableVertexAttribArray(pos_location);

    GLint tc_location = glGetAttribLocation(program, "texcoord");
    glVertexAttribPointer(tc_location, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void *) (2*sizeof(float)));
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
    glActiveTexture(GL_TEXTURE0);
    textureObject->Bind();

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    GLint loc_img0 = glGetUniformLocation(program, "img0");
    glUniform1i(loc_img0, 0);

    CheckGLError();

    GLint loc_texture_size = glGetUniformLocation(program, "texture_size");
    glUniform2f(loc_texture_size, float(img.width()), float(img.height()));

    GLint loc_cubic_flag = glGetUniformLocation(program, "use_cubic_interpolation");
    if (useCubicInterpolation) {
        std::cout << "[LOG] Using cubic interpolation to generate the new texture" << std::endl;
        glUniform1i(loc_cubic_flag, 1);
    } else
        glUniform1i(loc_cubic_flag, 0);

    QImage textureImage(img.width(), img.height(), QImage::Format_ARGB32);

    // disable depth and stencil test (if they were enabled) as the render target does not have the buffers attached
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 0, m.FN()*3);

    CheckGLError();

    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, textureImage.width(), textureImage.height(), GL_BGRA, GL_UNSIGNED_BYTE, textureImage.bits());

    glfwPollEvents();

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    if (sharedContext == false) textureObject->Release();
    glDeleteTextures(1, &renderTarget);
    glDeleteFramebuffers(1, &fbo);
    glDeleteBuffers(1, &vertexbuf);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);
    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }

    if (filter) vcg::PullPush(textureImage, qRgba(0, 255, 0, 255));

    TextureObjectHandle newTextureObject = std::make_shared<TextureObject>();
    newTextureObject->AddImage(std::make_shared<QImage>(textureImage.mirrored()));

    return newTextureObject;
}


#endif // TEXTURE_RENDERING_H

