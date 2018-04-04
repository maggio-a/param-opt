#include "gl_utils.h"

#include <cassert>
#include <iostream>
#include <memory>

#include <GLFW/glfw3.h>

#include <QImage>
#include <QFileInfo>

static void ErrorCallback(int error, const char *message)
{
    std::cout << "GLFW error: " << message << " (code " << error << ")" << std::endl;
}

void GLInit()
{
    if (!glfwInit()) {
        std::cout << "Failed to initialize glfw" << std::endl;
        exit(-1);
    }
    glfwSetErrorCallback(ErrorCallback);
}

void GLTerminate()
{
    glfwTerminate();
}

void CheckGLError()
{
    GLenum error = glGetError();
    if (error != GL_NO_ERROR)
    {
        std::cout << "OpenGL error " << error  << " ";
        if (error == GL_INVALID_VALUE) std::cout << "GL_INVALID_VALUE";
        if (error == GL_INVALID_OPERATION) std::cout << "GL_INVALID_OPERATION";
        std::cout << std::endl;
        //assert(error == GL_NO_ERROR);
    }
}

GLuint CompileShaders(const GLchar **vs_text, const GLchar **fs_text)
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



// TextureObject class implementation
// ==================================

TextureObject::TextureObject()
    : _texture{0}
{
}

TextureObject::~TextureObject()
{
    Release();
}

void TextureObject::AddImage(std::shared_ptr<QImage> img)
{
    imgVec.push_back(img);
}

/// !!! NOTE !!!
/// using GL_BGRA when reading from-to QImages as the QImage format is AA RR GG BB
/// and on little endian machines this is stored as [BB] [GG] [RR] [AA] when reinterpreted
/// as a pointer to the first byte of the four components
void TextureObject::Bind()
{
    assert(imgVec.size() > 0);
    // load texture from qimage on first use
    if (_texture == 0) {
        QImage& img = *imgVec[0];
        glGenTextures(1, &_texture);

        GLenum format, channels, type;

        switch(img.format()) {
        case QImage::Format_RGB32:
        case QImage::Format_ARGB32:
            format = GL_RGBA8; channels = GL_BGRA; type = GL_UNSIGNED_BYTE; break;
        default:
            std::cout << "Unsupported texture format" << std::endl; std::exit(-1);
        }
        QImage mirrored = img.mirrored(); // mirror to match opengl convention
        glBindTexture(GL_TEXTURE_2D, _texture);
        glTexStorage2D(GL_TEXTURE_2D, 1, format, img.width(), img.height());
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width(), img.height(), channels, type, mirrored.constBits());
        CheckGLError();
    }
    else {
        glBindTexture(GL_TEXTURE_2D, _texture);
        CheckGLError();
    }
}

void TextureObject::Release()
{
    if (_texture) {
        glDeleteTextures(1, &_texture);
        _texture = 0;
    }
}

int TextureObject::TextureWidth(std::size_t i)
{
    assert(i < imgVec.size());
    return imgVec[i]->width();
}

int TextureObject::TextureHeight(std::size_t i)
{
    assert(i < imgVec.size());
    return imgVec[i]->height();
}

std::size_t TextureObject::ArraySize() {
    return imgVec.size();
}

