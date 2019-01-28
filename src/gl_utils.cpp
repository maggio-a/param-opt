#include "gl_utils.h"
#include "logging.h"

#include <cassert>
#include <iostream>
#include <memory>

#include <GLFW/glfw3.h>

#include <QImage>
#include <QFileInfo>

static void ErrorCallback(int error, const char *message)
{
    LOG_ERR << "GLFW error: " << message << " (code " << error << ")";
}

void GLInit()
{
    if (!glfwInit()) {
        LOG_ERR << "Failed to initialize GLFW";
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
        std::stringstream ss;
        ss << "OpenGL error " << error  << " ";
        if (error == GL_INVALID_VALUE) ss << "GL_INVALID_VALUE";
        if (error == GL_INVALID_OPERATION) ss << "GL_INVALID_OPERATION";
        LOG_ERR << ss.str();
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
        LOG_DEBUG << infoLog;
        memset(infoLog, 0, 1024);
    }
    glGetShaderiv(vs, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        LOG_ERR << "Vertex shader compilation failed";
    }

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, fs_text, NULL);
    glCompileShader(fs);
    glGetShaderInfoLog(fs, 1024, NULL, infoLog);
    if (*infoLog) {
        LOG_DEBUG << infoLog;
        memset(infoLog, 0, 1024);
    }
    glGetShaderiv(fs, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        LOG_ERR << "Fragment shader compilation failed";
    }

    GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);
    glValidateProgram(program);
    glGetProgramInfoLog(program, 1024, NULL, infoLog);
    if (*infoLog) {
        LOG_DEBUG << infoLog;
    }
    glGetProgramiv(program, GL_LINK_STATUS, &status);
    if (status == GL_FALSE) {
        LOG_ERR << "Shader program link failed";
    }

    glDeleteShader(vs);
    glDeleteShader(fs);

    CheckGLError();

    return program;
}

void Mirror(QImage& img)
{
    int i = 0;
    while (true) {
        QRgb *line0 = (QRgb *) img.scanLine(i);
        QRgb *line1 = (QRgb *) img.scanLine(img.height() - 1 - i);
        i++;
        for (int j = 0; j < img.width(); ++j)
            std::swap(line0[j], line1[j]);
        if (i > img.height() / 2) break;
    }
}


// TextureObject class implementation
// ==================================

TextureObject::TextureObject()
{
}

TextureObject::~TextureObject()
{
    for (std::size_t i = 0; i < texNameVec.size(); ++i)
        Release(i);
}

void TextureObject::AddImage(std::shared_ptr<QImage> img)
{
    imgVec.push_back(img);
    texNameVec.push_back(0);
}

/// !!! NOTE !!!
/// using GL_BGRA when reading from-to QImages as the QImage format is AA RR GG BB
/// and on little endian machines this is stored as [BB] [GG] [RR] [AA] when reinterpreted
/// as a pointer to the first byte of the four components
void TextureObject::Bind(int i)
{
    assert(i >= 0 && i < (int)imgVec.size());
    // load texture from qimage on first use
    if (texNameVec[i] == 0) {

        /*
        glGenTextures(1, &texNameVec[i]);

        // generate 3 mip levels
        std::vector<unsigned> l0(64*64, 0xffffffff);
        std::vector<unsigned> l1(1024 , 0xff0000ff);
        std::vector<unsigned> l2(256  , 0xffff0000);
        std::vector<unsigned> l3(64   , 0xff00ff00);
        std::vector<unsigned> l4(16   , 0xffff00ff);
        std::vector<unsigned> l5(4    , 0xff00ffff);
        std::vector<unsigned> l6(1    , 0xffffff00);


        GLenum format, channels, type;
        format = GL_RGBA8; channels = GL_BGRA; type = GL_UNSIGNED_BYTE;
        glBindTexture(GL_TEXTURE_2D, texNameVec[i]);
        glTexStorage2D(GL_TEXTURE_2D, 7, format, 64, 64);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 64, 64, channels, type, l0.data());
        glTexSubImage2D(GL_TEXTURE_2D, 1, 0, 0, 32, 32, channels, type, l1.data());
        glTexSubImage2D(GL_TEXTURE_2D, 2, 0, 0, 16, 16, channels, type, l2.data());
        glTexSubImage2D(GL_TEXTURE_2D, 3, 0, 0, 8, 8, channels, type, l3.data());
        glTexSubImage2D(GL_TEXTURE_2D, 4, 0, 0, 4, 4, channels, type, l4.data());
        glTexSubImage2D(GL_TEXTURE_2D, 5, 0, 0, 2, 2, channels, type, l5.data());
        glTexSubImage2D(GL_TEXTURE_2D, 6, 0, 0, 1, 1, channels, type, l6.data());


        CheckGLError();
        */





        QImage& img = *imgVec[i];
        glGenTextures(1, &texNameVec[i]);

        GLenum format, channels, type;

        switch(img.format()) {
        case QImage::Format_RGB32:
        case QImage::Format_ARGB32:
            format = GL_RGBA8; channels = GL_BGRA; type = GL_UNSIGNED_BYTE; break;
        default:
            LOG_ERR << "Unsupported texture format " << img.format();
            std::exit(-1);
        }
        Mirror(img);
        //QImage mirrored = img.mirrored(); // mirror to match opengl convention
        int miplevels = std::log2((float) img.width());
        glBindTexture(GL_TEXTURE_2D, texNameVec[i]);
        glTexStorage2D(GL_TEXTURE_2D, miplevels, format, img.width(), img.height());
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width(), img.height(), channels, type, img.constBits());
        glGenerateMipmap(GL_TEXTURE_2D);
        CheckGLError();
        Mirror(img); // Mirror the QImage back so that its data view remains consistent
    }
    else {
        glBindTexture(GL_TEXTURE_2D, texNameVec[i]);
        CheckGLError();
    }
}

void TextureObject::Release(int i)
{
    assert(i >= 0 && i < (int)imgVec.size());
    if (texNameVec[i]) {
        glDeleteTextures(1, &texNameVec[i]);
        texNameVec[i] = 0;
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

int64_t TextureObject::TextureArea(std::size_t i)
{
    assert(i < imgVec.size());
    return ((int64_t) TextureWidth(i)) * TextureHeight(i);
}

std::size_t TextureObject::ArraySize()
{
    return imgVec.size();
}

