#ifndef PARAMETERIZATION_CHECKER_H
#define PARAMETERIZATION_CHECKER_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>

#include <QImage>

#include "mesh.h"
#include "mesh_graph.h"
#include "gl_utils.h"

struct RasterizedParameterizationStats {
    int totalFragments;
    int overwrittenFragments; // number of fragments that were written more than once
    int lostFragments; // number of fragments lost due to overwrites: if fragment f has fw>1 writes, than lostFragmens += (fw-1)
    int boundaryFragments;
};

static const char *vs_text_checker[] = {
    "#version 410 core                                           \n"
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
    "#version 410 core                                           \n"
    "                                                            \n"
    "out vec4 color;                                             \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    color = vec4(1.0, 1.0, 1.0, 1.0);                       \n"
    "}                                                           \n"
};

static RasterizedParameterizationStats GetRasterizationStats(ChartHandle chart, int width, int height)
{
    // Create a hidden window
    GLFWwindow *parentWindow = glfwGetCurrentContext();
    bool sharedContext = (parentWindow != nullptr);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    GLFWwindow *window = glfwCreateWindow(width, height, "Window", nullptr, parentWindow);
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

    glfwSwapInterval(1);

    int fbw, fbh;
    glfwGetFramebufferSize(window, &fbw, &fbh);
    assert(fbw == width && fbh == height);

    // OpenGL setup

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLint program = CompileShaders(vs_text_checker, fs_text_checker);
    glUseProgram(program);

    CheckGLError();

    // Allocate vertex data

    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    vcg::Box2d box = chart->UVBox();

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, chart->FN()*6*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto fptr : chart->fpVec) {
        for (int i = 0; i < 3; ++i) {
            // normalize coordinates
            *p++ = (fptr->cWT(i).U() - box.min.X()) / box.DimX();
            *p++ = (fptr->cWT(i).V() - box.min.Y()) / box.DimY();
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint pos_location = glGetAttribLocation(program, "position");
    glVertexAttribPointer(pos_location, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(pos_location);

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glViewport(0, 0, width, height);

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);

    glStencilFunc(GL_ALWAYS, 1, 0xFF);
    glStencilOp(GL_KEEP, GL_INCR, GL_INCR); // increment stencil regardless of depth pass|fail
    glStencilMask(0xFF);

    glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 0, chart->FN()*3);

    CheckGLError();

    //glReadBuffer(GL_COLOR_ATTACHMENT0);
    unsigned char sb[width*height];
    glReadPixels(0, 0, width, height, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, sb);

    RasterizedParameterizationStats stats{0, 0, 0, 0};

    bool outside;
    for (int j = 0; j < height; ++j) {
        outside = true;
        for (int i = 0; i < width; ++i) {
            int k = j * width + i;
            int n = sb[k]; // stencil value
            if (n > 0) {
                stats.totalFragments += n;
                if (n > 1) stats.overwrittenFragments++;
                stats.lostFragments += n - 1;
                if (outside     // enering a segment
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

    glfwPollEvents();

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    //glDeleteTextures(1, &renderTarget);
    //glDeleteFramebuffers(1, &fbo);
    glDeleteBuffers(1, &vertexbuf);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);

    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }

    return stats;
}


#endif // PARAMETERIZATION_CHECKER_H

