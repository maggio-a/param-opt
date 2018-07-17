#ifndef PARAMETERIZATION_CHECKER_H
#define PARAMETERIZATION_CHECKER_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <fstream>

#include <QImage>

#include "mesh.h"
#include "mesh_graph.h"
#include "gl_utils.h"

struct RasterizedParameterizationStats {
    int rw; // raster width
    int rh; // raster height
    int totalFragments;
    int totalFragments_bilinear;
    int overwrittenFragments; // number of fragments that were written more than once
    int lostFragments; // number of fragments lost due to overwrites: if fragment f has fw>1 writes, than lostFragmens += (fw-1)
    int boundaryFragments;
};

static const char *vs_text_checker[] = {
    "#version 430 core                                           \n"
    "                                                            \n"
    "in vec2 position;                                           \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);               \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                        \n"
    "}                                                           \n"
};

static const char *fs_text_checker[] = {
    "#version 430 core                                                \n"
    "                                                                 \n"
    "layout (r32ui) uniform uimage2D imgbuf;                          \n"
    "out vec4 color;                                                  \n"
    "                                                                 \n"
    "void main(void)                                                  \n"
    "{                                                                \n"
    "    color = vec4(1.0, 1.0, 1.0, 1.0);                            \n"
    "    imageAtomicAdd(imgbuf, ivec2(gl_FragCoord.xy), 1);           \n"
    "}                                                                \n"
};

static bool CheckLocalInjectivity(Mesh& m)
{
    std::cout << "TODO" << std::endl;
    return true;
}

static bool CheckUVConnectivity(Mesh& m)
{
    std::cout << "TODO" << std::endl;
    return true;
}

static RasterizedParameterizationStats GetRasterizationStats(const std::vector<Mesh::FacePointer>& faces, int width, int height, const vcg::Box2d& resizeBox);

static RasterizedParameterizationStats GetRasterizationStats(Mesh& m, int width, int height)
{
    vcg::Box2d nullbox;
    nullbox.SetNull();
    std::vector<Mesh::FacePointer> faces;
    for (auto & f : m.face) { faces.push_back(&f); }
    return GetRasterizationStats(faces, width, height, nullbox);
}

static RasterizedParameterizationStats GetRasterizationStats(ChartHandle chart, int width, int height)
{
    return GetRasterizationStats(chart->fpVec, width, height, chart->UVBox());
}

static RasterizedParameterizationStats GetRasterizationStats(const std::vector<Mesh::FacePointer>& faces, int width, int height, const vcg::Box2d& resizeBox)
{
    std::cout << "RASTER STATS RESOLUTION = " << width << "x" << height << std::endl;
    // Create a hidden window
    GLFWwindow *parentWindow = glfwGetCurrentContext();
    bool sharedContext = (parentWindow != nullptr);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    //glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);

    //GLFWwindow *window = glfwCreateWindow(width, height, "Window", nullptr, parentWindow);
    GLFWwindow *window = glfwCreateWindow(512, 512, "Window", nullptr, parentWindow);
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
    glGetError();

    int fbw, fbh;
    glfwGetFramebufferSize(window, &fbw, &fbh);
    //assert(fbw == width && fbh == height);

    // OpenGL setup

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLint program = CompileShaders(vs_text_checker, fs_text_checker);
    glUseProgram(program);

    // Allocate vertex data

    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    //vcg::Box2d box = chart->UVBox();

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, faces.size()*6*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto fptr : faces) {
        for (int i = 0; i < 3; ++i) {
            if (! resizeBox.IsNull()) {
                // normalize coordinates
                *p++ = (fptr->cWT(i).U() - resizeBox.min.X()) / resizeBox.DimX();
                *p++ = (fptr->cWT(i).V() - resizeBox.min.Y()) / resizeBox.DimY();
            } else {
                *p++ = fptr->cWT(i).U();
                *p++ = fptr->cWT(i).V();
            }
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint pos_location = glGetAttribLocation(program, "position");
    glVertexAttribPointer(pos_location, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(pos_location);

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Create texture to store fragment writes
    constexpr int imgbuf_unit = 0;
    GLint loc_imgbuf = glGetUniformLocation(program, "imgbuf");
    glUniform1i(loc_imgbuf, imgbuf_unit);

    assert(sizeof(unsigned) == 4);
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32UI, width, height);

    // clear the texure
    unsigned *sb = new unsigned[width*height];
    memset(sb, 0x00, width*height*sizeof(unsigned));
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RED_INTEGER, GL_UNSIGNED_INT, sb);

    glBindImageTexture(imgbuf_unit, tex, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32UI);

    // bugfix (?) setup an offscreen context of the appropriate size to make sure
    // the data is fully rendered
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    GLuint renderTarget;
    glGenTextures(1, &renderTarget);
    glBindTexture(GL_TEXTURE_2D, renderTarget);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, width, height);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, renderTarget, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glClearColor(0.0f, 1.0f, 0.0f, 1.0f);

    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glClear(GL_COLOR_BUFFER_BIT);
    glfwPollEvents();

    glDrawArrays(GL_TRIANGLES, 0, faces.size()*3);

    glBindTexture(GL_TEXTURE_2D, tex);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, sb);

    CheckGLError();

    //glReadBuffer(GL_BACK);
    //glReadPixels(0, 0, width, height, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, sb);

    RasterizedParameterizationStats stats{0, 0, 0, 0, 0, 0, 0};
    stats.rw = width;
    stats.rh = height;
    bool outside;
    for (int j = 0; j < height; ++j) {
        outside = true;
        for (int i = 0; i < width; ++i) {
            int k = j * width + i;
            //std::cout << k << std::endl;
            int n = sb[k]; // stencil value
            if (n > 0) {
                stats.totalFragments += n;
                stats.totalFragments_bilinear += 1;
                if (n > 1) {
                    stats.overwrittenFragments++;
                    stats.lostFragments += n - 1;
                }
                if (outside     // entering a segment
                    || j == 0   // first line in the raster
                    || j == (height-1) // last line in the raster
                    || i == (width-1) // at the right boundary of the raster
                    || sb[(j+1)*width+i] == 0 // pixel above is empty
                    || sb[(j-1)*width+i] == 0 // pixel below is empty
                    || sb[k+1] == 0 // pixel to the right is empty
                        ) { // then it is a boundary fragment
                    stats.boundaryFragments++;
                }
                outside = false;
            } else {
                outside = true;

                // inspect a 3x3 kernel and add pixel to bilinear area if it overlaps a fragment
                bool bilinear = false;
                for(int ker_j = std::max(j-1, 0); ker_j < std::min(j+1, height); ++ker_j) {
                    for(int ker_i = std::max(i-1, 0); ker_i < std::min(i+1, width); ++ker_i) {
                        if (sb[ker_j * width + ker_i] > 0) {
                            bilinear = true;
                        }
                    }
                }
                if (bilinear) stats.totalFragments_bilinear++;
            }
        }
    }

    /*
    ofstream ofs("overlap.ppm", ios_base::binary | ios_base::out | ios_base::trunc);
    if (!ofs) {
        assert(0);
    } else {
        ofs << "P6 " << width << " " << height << " 255\n";
        for (int i = 0; i < width*height; ++i) {
            if (sb[i] > 0) {
                ofs << (unsigned char) 255 << (unsigned char) 255 << (unsigned char) 255;
            } else {
                ofs << (unsigned char) 255 << (unsigned char) 0 << (unsigned char) 0;
            }
        }
        ofs.close();
    }
    */

    delete sb;
    sb = nullptr;

    glfwPollEvents();

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteBuffers(1, &vertexbuf);
    glDeleteTextures(1, &tex);
    glDeleteTextures(1, &renderTarget);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);

    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }

    return stats;
}


#endif // PARAMETERIZATION_CHECKER_H

