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

/* Vertically mirrors a QImage in-place, useful to match the OpenGL convention
 * for texture data storage */
void Mirror(QImage& img);

/* wrapper to an (eventually) array of textures */
struct TextureObject {

    std::vector<std::shared_ptr<QImage>> imgVec;
    std::vector<GLuint> texNameVec;

    TextureObject();
    ~TextureObject();

    TextureObject(const TextureObject &) = delete;
    TextureObject &operator=(const TextureObject &) = delete;

    /* Add QImage ref to the texture object */
    void AddImage(std::shared_ptr<QImage> img);

    /* Binds the texture at index i */
    void Bind(int i);

    /* Releases the texture i, without unbinding it if it is bound */
    void Release(int i);

    int TextureWidth(std::size_t i);
    int TextureHeight(std::size_t i);

    int64_t TextureArea(std::size_t i);

    std::size_t ArraySize();
};


#endif // GL_UTIL_H

