#ifndef GL_UTIL_H
#define GL_UTIL_H

#include <iostream>

#include <GL/glew.h>
#include <QImage>

static void ErrorCallback(int error, const char *message)
{
    std::cout << "GLFW error: " << message << " (code " << error << ")" << std::endl;
}

static void CheckGLError()
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

/// !!! NOTE !!!
/// using GL_BGRA when reading from-to QImages as the QImage format is AA RR GG BB
/// and on little endian machines this is stored as [BB] [GG] [RR] [AA] when reinterpreted
/// as a pointer to the first byte of the four components
static void LoadTexture2DFromQImage(const QImage& img, GLuint texture)
{
    GLenum format, channels, type;

    switch(img.format()) {
    case QImage::Format_RGB32:
    case QImage::Format_ARGB32:
        format = GL_RGBA8; channels = GL_BGRA; type = GL_UNSIGNED_BYTE; break;
    default:
        cout << "Unupported texture format" << endl; exit(-1);
    }

    QImage mirrored = img.mirrored();

    glBindTexture(GL_TEXTURE_2D, texture);
    glTexStorage2D(GL_TEXTURE_2D, 1, format, img.width(), img.height());
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width(), img.height(), channels, type, mirrored.constBits());

    CheckGLError();
}

#endif // GL_UTIL_H

