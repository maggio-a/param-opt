#ifndef GL_UTILS_H
#define GL_UTILS_H

#include <vector>
#include <memory>

#include <GL/glew.h>

class QImage;
class TextureObject;

using TextureObjectHandle = std::shared_ptr<TextureObject>;

/* Initializes glfw */
void GLInit();

/* Terminates glfw */
void GLTerminate();

/* Prints the last OpenGL error code */
void CheckGLError();

/* Compiles a vertex shader source and a fragment shader source into a program */
GLuint CompileShaders(const GLchar **vs_text, const GLchar **fs_text);

/* wrapper to an (eventually) array texture */
struct TextureObject {

    std::vector<std::shared_ptr<QImage>> imgVec;
    GLuint _texture;

    TextureObject();
    ~TextureObject();

    TextureObject(const TextureObject &) = delete;
    TextureObject &operator=(const TextureObject &) = delete;

    /* Add QImage ref to the texture object */
    void AddImage(std::shared_ptr<QImage> img);

    /* Binds the texture object, creating a new opengl texture from the imgVec array */
    void Bind();

    /* Deletes the current opengl texture, without unbinding it if it is bound */
    void Release();

    int TextureWidth(std::size_t i);
    int TextureHeight(std::size_t i);

    std::size_t ArraySize();
};


#endif // GL_UTIL_H

