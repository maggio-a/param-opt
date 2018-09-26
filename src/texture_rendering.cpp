/*
 * References for the bicubic interpolated texture lookup:
 *  - GPU gems 2 ch 20 (Sigg and Hadwiger 2005)
 *  - Efficient GPU-Based Texture Interpolation using Uniform B-Splines  (Ruijters et al. 2009)
 * */

#include "texture_rendering.h"
#include "mesh.h"
#include "pushpull.h"
#include "uv.h"
#include "mesh_attribute.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <algorithm>

#include <QImage>


/* === Shaders ============================================================== */

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
    "    //if (uv.s < 0) uv = vec2(0.0, 0.0);                    \n"
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
    "    if (use_cubic_interpolation == 0) {                               \n"
    "        if (uv.s < 0)                                                 \n"
    "            color = vec4(0, 1, 0, 1);                                 \n"
    "        else                                                          \n"
    "            color = texture2D(img0, uv);                              \n"
    "    } else {                                                          \n"
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


static const char *vs_text_checker[] = {
    "#version 430 core                                           \n"
    "                                                            \n"
    "in vec2 position;                                           \n"
    "in uint chart_id;                                           \n"
    "out flat uint chartId;                                      \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);               \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                        \n"
    "    chartId = chart_id;                                     \n"
    "}                                                           \n"
};

static const char *fs_text_checker[] = {
    "#version 430 core                                                                \n"
    "                                                                                 \n"
    "layout (r32ui) uniform uimage2D imgbuf;                                          \n"
    "layout (r32ui) uniform uimage2D idbuf;                                           \n"
    "in flat uint chartId;                                                            \n"
    "out vec4 color;                                                                  \n"
    "                                                                                 \n"
    "void main(void)                                                                  \n"
    "{                                                                                \n"
    "    color = vec4(1.0, 1.0, 1.0, 1.0);                                            \n"
    "    imageAtomicAdd(imgbuf, ivec2(gl_FragCoord.xy), 1);                           \n"
    "    imageAtomicExchange(idbuf, ivec2(gl_FragCoord.xy), chartId);                 \n"
    #if 0
    "    uint val = imageAtomicCompSwap(idbuf, ivec2(gl_FragCoord.xy), 0, chartId);   \n"
    "    if (val == 0 || val == 0xffffffff)                                           \n" // no one set the vat for this fragment yet
    "        return;                                                                  \n" // or its a 'marked' pixel already
    "    else if (val != chartId)                                                     \n"
    "        imageAtomicExchange(idbuf, ivec2(gl_FragCoord.xy), 0xffffffff);          \n" // pixel was set to a different chart, mark it
    #endif
    "}                                                                                \n"
};


static const char *vs_text_geom_image[] = {
    "#version 430 core                                           \n"
    "                                                            \n"
    "in vec2 position;                                           \n"
    "in vec3 position_3d;                                        \n"
    "out vec3 p3d;                                               \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);               \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                        \n"
    "    p3d = position_3d;                                      \n"
    "}                                                           \n"
};

static const char *fs_text_geom_image[] = {
    "#version 430 core                                                       \n"
    "                                                                        \n"
    "layout (r32f) uniform image2D p_x;                                      \n"
    "layout (r32f) uniform image2D p_y;                                      \n"
    "layout (r32f) uniform image2D p_z;                                      \n"
    "layout (r32ui) uniform uimage2D img_fault;                              \n"
    "in vec3 p3d;                                                            \n"
    "uniform float magic_threshold;                                          \n"
    "uniform int pass;                                                       \n"
    "                                                                        \n"
    "void main(void)                                                         \n"
    "{                                                                       \n"
    "    if (pass == 1) {                                                    \n"
    "        imageStore(p_x, ivec2(gl_FragCoord.xy), vec4(p3d.x, 0, 0, 0));  \n"
    "        imageStore(p_y, ivec2(gl_FragCoord.xy), vec4(p3d.y, 0, 0, 0));  \n"
    "        imageStore(p_z, ivec2(gl_FragCoord.xy), vec4(p3d.z, 0, 0, 0));  \n"
    "    } else if (pass == 2) {                                             \n"
    "        vec4 x = imageLoad(p_x, ivec2(gl_FragCoord.xy));                \n"
    "        vec4 y = imageLoad(p_y, ivec2(gl_FragCoord.xy));                \n"
    "        vec4 z = imageLoad(p_z, ivec2(gl_FragCoord.xy));                \n"
    "        vec3 gimage_pos = vec3(x.x, y.x, z.x);                          \n"
    "        if (distance(p3d, gimage_pos) > magic_threshold) {              \n"
    "            imageAtomicAdd(img_fault, ivec2(gl_FragCoord.xy), 1);       \n"
    "        }                                                               \n"
    "    }                                                                   \n"
    "}                                                                       \n"
};


/* === Functions ============================================================ */

static RasterizedParameterizationStats GetRasterizationStats(Mesh& m, const std::vector<Mesh::FacePointer>& faces, int width, int height);

static int FacesByTextureIndex(Mesh& m, std::vector<std::vector<Mesh::FacePointer>>& fv);

static std::shared_ptr<QImage> RenderTexture(std::vector<Mesh::FacePointer>& fvec,
                                             Mesh &m, TextureObjectHandle textureObject,
                                             bool filter, InterpolationMode imode,
                                             int textureWidth, int textureHeight, GLFWwindow *parentWindow);



#include <wrap/io_trimesh/export.h>
static GeometryImageStats GetGeometryImageStats(Mesh& m, const std::vector<Mesh::FacePointer>& faces, int width, int height)
{
    assert(sizeof(unsigned) == 4);

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

    // OpenGL setup

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLint program = CompileShaders(vs_text_geom_image, fs_text_geom_image);
    glUseProgram(program);

    // Allocate vertex data

    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);


    // TODO scale

    //vcg::Box2d box = chart->UVBox();

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, faces.size()*15*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto fptr : faces) {
        for (int i = 0; i < 3; ++i) {
            *p++ = (float) fptr->cWT(i).U();
            *p++ = (float) fptr->cWT(i).V();
            *p++ = (float) fptr->P(i).X();
            *p++ = (float) fptr->P(i).Y();
            *p++ = (float) fptr->P(i).Z();
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint pos_location = glGetAttribLocation(program, "position");
    glVertexAttribPointer(pos_location, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
    glEnableVertexAttribArray(pos_location);

    GLint pos_id = glGetAttribLocation(program, "position_3d");
    glVertexAttribPointer(pos_id, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), (const GLvoid *) (2*sizeof(float)));
    glEnableVertexAttribArray(pos_id);

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0);


    vcg::Box3d bb;
    for (auto fptr : faces) {
        bb.Add(fptr->P(0));
        bb.Add(fptr->P(1));
        bb.Add(fptr->P(2));
    }

    GLint loc_inj_threshold = glGetUniformLocation(program, "magic_threshold");
    glUniform1f(loc_inj_threshold, bb.Diag()*0.001);

    unsigned *sb = new unsigned[width*height](); // zero initialize


    // textures for position components

    const int imgunit_p[] = {0, 1, 2};

    GLint loc_imgunit_p;
    loc_imgunit_p = glGetUniformLocation(program, "p_x");
    glUniform1i(loc_imgunit_p, imgunit_p[0]);

    loc_imgunit_p = glGetUniformLocation(program, "p_y");
    glUniform1i(loc_imgunit_p, imgunit_p[1]);

    loc_imgunit_p = glGetUniformLocation(program, "p_z");
    glUniform1i(loc_imgunit_p, imgunit_p[2]);

    GLuint tex_p[3];
    glGenTextures(3, tex_p);

    for (int i = 0; i < 3; ++i) {
        glBindTexture(GL_TEXTURE_2D, tex_p[i]);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32F, width, height);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RED, GL_FLOAT, sb);
        glBindImageTexture(imgunit_p[i], tex_p[i], 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32F);
    }

    // texture to store interpolation faults
    constexpr int faultbuf_unit = 3;
    GLint loc_faultbuf = glGetUniformLocation(program, "img_fault");
    glUniform1i(loc_faultbuf, faultbuf_unit);
    GLuint tex_fault;
    glGenTextures(1, &tex_fault);
    glBindTexture(GL_TEXTURE_2D, tex_fault);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32UI, width, height);

    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RED_INTEGER, GL_UNSIGNED_INT, sb);
    glBindImageTexture(faultbuf_unit, tex_fault, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32UI);


    // bugfix (?) setup an offscreen context of the appropriate size to make sure
    // the data is fully rendered
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glFramebufferParameteri(GL_FRAMEBUFFER, GL_FRAMEBUFFER_DEFAULT_WIDTH, width);
    glFramebufferParameteri(GL_FRAMEBUFFER, GL_FRAMEBUFFER_DEFAULT_HEIGHT, height);

    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glfwPollEvents();

    GLint loc_pass = glGetUniformLocation(program, "pass");

    // First pass
    glUniform1i(loc_pass, 1);
    glDrawArrays(GL_TRIANGLES, 0, faces.size()*3);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    // Second pass
    glUniform1i(loc_pass, 2);
    glDrawArrays(GL_TRIANGLES, 0, faces.size()*3);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glBindTexture(GL_TEXTURE_2D, tex_fault);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, sb);

#if 0

    float * x = new float[width * height]();
    float * y = new float[width * height]();
    float * z = new float[width * height]();
    glBindTexture(GL_TEXTURE_2D, tex_p[0]);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, x);

    glBindTexture(GL_TEXTURE_2D, tex_p[1]);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, y);

    glBindTexture(GL_TEXTURE_2D, tex_p[2]);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, z);


    std::cout << "REMOVE MEEEEEEEEEEEEEEEEEEEE" << std::endl;
    Mesh pc;
    Box3f bb2;
    int u = 0;
    for (int i = 0; i < height; ++i) {
        for (int j =  0; j < width; ++j) {
            int k = i * width + j;
            if (x[k] != 0) {
                u++;
                Point3f pf(x[k], y[k], z[k]);
                if(k%2) tri::Allocator<Mesh>::AddVertex(pc, Point3d(pf.X(), pf.Y(), pf.Z()));
                bb2.Add(pf);
            }
        }
    }
    tri::io::Exporter<Mesh>::Save(pc, "pc.obj", tri::io::Mask::IOM_VERTCOORD);
    std::cout << "############################### " << bb2.Diag() << "           " << u <<std::endl;
    delete[] x;
    delete[] y;
    delete[] z;
#endif

    CheckGLError();

    //glReadBuffer(GL_BACK);
    //glReadPixels(0, 0, width, height, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, sb);

    GeometryImageStats stats = {};

    stats.rw = width;
    stats.rh = height;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int k = i * width + j;
            int n = sb[k]; // fault value
            if (n > 0) {
                stats.uniqueFaults++;
                stats.totalFaults += (n + 1);
            }
        }
    }

    glfwPollEvents();

    delete [] sb;

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteBuffers(1, &vertexbuf);
    glDeleteTextures(1, &tex_fault);
    glDeleteTextures(3, tex_p);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);

    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }

    return stats;
}

std::vector<GeometryImageStats> GetGeometryImageStats(Mesh& m, TextureObjectHandle textureObject)
{
    std::vector<std::vector<Mesh::FacePointer>> facesByTexture;
    int ntex = FacesByTextureIndex(m, facesByTexture);

    std::vector<GeometryImageStats> statsVec;
    for (int i = 0; i < ntex; ++i) {
        int tw = textureObject->TextureWidth(i);
        int th = textureObject->TextureHeight(i);
        while (std::min(tw, th) > 4096) {
            tw /= 2;
            th /= 2;
        }
        GeometryImageStats stats = GetGeometryImageStats(m, facesByTexture[i], tw, th);
        statsVec.push_back(stats);
    }

    return statsVec;
}

static int FacesByTextureIndex(Mesh& m, std::vector<std::vector<Mesh::FacePointer>>& fv)
{
    fv.clear();

    // Detect the number of required textures
    int nTex = 1;
    for (auto&f : m.face) {
        nTex = std::max(nTex, f.cWT(0).N() + 1);
    }

    fv.resize(nTex);

    for (auto& f : m.face) {
        int ti = f.cWT(0).N();
        assert(ti < nTex);
        fv[ti].push_back(&f);
    }

    return fv.size();
}


RasterizedParameterizationStats GetRasterizationStats(ChartHandle chart, int width, int height)
{
    Box2d uvBox = chart->UVBox();
    std::vector<TexCoordStorage> tcVec(chart->FN());
    for (unsigned k = 0; k < chart->FN(); ++k) {
        for (int i = 0; i < 3; ++i) {
            tcVec[k].tc[i] = chart->fpVec[k]->WT(i);
            chart->fpVec[k]->WT(i).U() = (chart->fpVec[k]->WT(i).U() - uvBox.min.X()) / uvBox.DimX();
            chart->fpVec[k]->WT(i).V() = (chart->fpVec[k]->WT(i).V() - uvBox.min.Y()) / uvBox.DimY();
        }
    }

    RasterizedParameterizationStats stats = GetRasterizationStats(chart->mesh, chart->fpVec, width, height);

    for (unsigned k = 0; k < chart->FN(); ++k)
        for (int i = 0; i < 3; ++i)
            chart->fpVec[k]->WT(i) = tcVec[k].tc[i];

    return stats;
}

std::vector<RasterizedParameterizationStats> GetRasterizationStats(Mesh& m, TextureObjectHandle textureObject)
{
    std::vector<std::vector<Mesh::FacePointer>> facesByTexture;
    int ntex = FacesByTextureIndex(m, facesByTexture);

    std::vector<RasterizedParameterizationStats> statsVec;
    for (int i = 0; i < ntex; ++i) {
        RasterizedParameterizationStats stats = GetRasterizationStats(m, facesByTexture[i], textureObject->TextureWidth(i), textureObject->TextureHeight(i));
        statsVec.push_back(stats);
    }

    return statsVec;
}

std::vector<std::vector<RasterizedParameterizationStats>> GetRasterizationStatsAtMipmapLevels(Mesh& m, TextureObjectHandle textureObject)
{
    constexpr int MIN_DIM = 4;
    std::vector<std::vector<Mesh::FacePointer>> facesByTexture;
    int ntex = FacesByTextureIndex(m, facesByTexture);

    std::vector<std::vector<RasterizedParameterizationStats>> levelStats;
    for (int i = 0; i < ntex; ++i) {
        std::vector<RasterizedParameterizationStats> textureStats;
        int tw = textureObject->TextureWidth(i);
        int th = textureObject->TextureHeight(i);
        while (std::max(tw, th) >= MIN_DIM) {
            RasterizedParameterizationStats stats = GetRasterizationStats(m, facesByTexture[i], tw, th);
            textureStats.push_back(stats);
            tw /= 2;
            th /= 2;
        }
        levelStats.push_back(textureStats);
    }

    return levelStats;
}

enum KernelBit {
    Unset = 1,
    Set = 2
};

inline int Inspect3x3Kernel(unsigned *buffer, int row, int col, int width, int height)
{
    int mask = 0;
    if (row == 0 || row == height - 1)
        mask |= Unset;
    if (col == 0 || col == width - 1)
        mask |= Unset;
    for (int i = std::max(row-1, 0); i < std::min(row+2, height); ++i) {
        for (int j = std::max(col-1, 0); j < std::min(col+2, width); ++j) {
            if (buffer[i * width + j] == 0)
                mask |= Unset;
            else
                mask |= Set;
        }
    }
    return mask;
}

inline bool Inspect3x3KernelClash(unsigned *buffer, int row, int col, int width, int height)
{
    unsigned id = buffer[row*width + col];
    for (int i = std::max(row-1, 0); i < std::min(row+2, height); ++i) {
        for (int j = std::max(col-1, 0); j < std::min(col+2, width); ++j) {
            unsigned buff_id = buffer[i * width + j];
            if (i != row && j != col && (buff_id != 0 && buff_id != id))
                return true;
        }
    }
    return false;
}

static RasterizedParameterizationStats GetRasterizationStats(Mesh& m, const std::vector<Mesh::FacePointer>& faces, int width, int height)
{
    assert(sizeof(unsigned) == 4);

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

    assert(HasConnectedComponentIDAttribute(m));
    auto CCIDh = GetConnectedComponentIDAttribute(m);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, faces.size()*9*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto fptr : faces) {
        for (int i = 0; i < 3; ++i) {
            *p++ = fptr->cWT(i).U();
            *p++ = fptr->cWT(i).V();

            unsigned int *up = (unsigned int *) p;
            *up = (unsigned int) CCIDh[fptr] + 1;

            p++;
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint pos_location = glGetAttribLocation(program, "position");
    glVertexAttribPointer(pos_location, 2, GL_FLOAT, GL_FALSE, 3*sizeof(float), 0);
    glEnableVertexAttribArray(pos_location);

    GLint pos_id = glGetAttribLocation(program, "chart_id");
    glVertexAttribPointer(pos_id, 1, GL_UNSIGNED_INT, GL_FALSE, 3*sizeof(float), (const GLvoid *) (2*sizeof(float)));
    glEnableVertexAttribArray(pos_id);

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0);


    unsigned *sb = new unsigned[width*height](); // zero initialize

    // Create textures to store fragment writes

    // Texture to count overdraw
    constexpr int imgbuf_unit = 0;
    GLint loc_imgbuf = glGetUniformLocation(program, "imgbuf");
    glUniform1i(loc_imgbuf, imgbuf_unit);

    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32UI, width, height);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RED_INTEGER, GL_UNSIGNED_INT, sb);

    glBindImageTexture(imgbuf_unit, tex, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32UI);

    // Texture to store marks. The texture is initially zero-initialized, and each fragment
    // write checks if the stored pixel value is either zero, an 'invalid' value (0xffffffff), or a non-zero chart id.
    // If it reads zero, it writes its value to the pixel (with atomic compare and swap). Otherwise, if the comparison failed
    // and the returned value is invalid, nothing happens. If the returned value is another chart id instead, then it is checked
    // against the chart id of the fragment, and if it is different it means that two different regions overlap on the pizxel, and
    // is therefore set to the invalid value 0xffffffff
    // Note that since 0 is a perfectly valid id, the id of each face as a vertex attribute is incremented by one
    constexpr int idbuf_unit = 1;
    GLint loc_idbuf = glGetUniformLocation(program, "idbuf");
    glUniform1i(loc_idbuf, idbuf_unit);
    GLuint tex_id;
    glGenTextures(1, &tex_id);
    glBindTexture(GL_TEXTURE_2D, tex_id);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32UI, width, height);

    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RED_INTEGER, GL_UNSIGNED_INT, sb);
    glBindImageTexture(idbuf_unit, tex_id, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32UI);


    // bugfix (?) setup an offscreen context of the appropriate size to make sure
    // the data is fully rendered
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glFramebufferParameteri(GL_FRAMEBUFFER, GL_FRAMEBUFFER_DEFAULT_WIDTH, width);
    glFramebufferParameteri(GL_FRAMEBUFFER, GL_FRAMEBUFFER_DEFAULT_HEIGHT, height);

    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glClearColor(0.0f, 1.0f, 0.0f, 1.0f);

    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glClear(GL_COLOR_BUFFER_BIT);
    glfwPollEvents();

    glDrawArrays(GL_TRIANGLES, 0, faces.size()*3);

    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glBindTexture(GL_TEXTURE_2D, tex);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, sb);

    unsigned *mb = new unsigned[width*height];

    glBindTexture(GL_TEXTURE_2D, tex_id);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RED_INTEGER, GL_UNSIGNED_INT, mb);

    CheckGLError();

    //glReadBuffer(GL_BACK);
    //glReadPixels(0, 0, width, height, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, sb);

    RasterizedParameterizationStats stats = {};
    stats.rw = width;
    stats.rh = height;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int k = i * width + j;
            int n = sb[k]; // stencil value
            int kernelMask = Inspect3x3Kernel(sb, i, j, width, height);
            if (n > 0) {
                stats.totalFragments += n;
                stats.totalFragments_bilinear++;
                if ((n > 1)) {
                    stats.overwrittenFragments++;
                    stats.lostFragments += (n - 1);
                }

                if (Inspect3x3KernelClash(mb, i, j, width, height))
                    stats.fragmentClashes++;

                if (kernelMask &= KernelBit::Unset)
                    stats.boundaryFragments++;
            } else {
                if (kernelMask &= KernelBit::Set)
                    stats.totalFragments_bilinear++;
            }
        }
    }

    glfwPollEvents();

    delete [] sb;
    delete [] mb;

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteBuffers(1, &vertexbuf);
    glDeleteTextures(1, &tex);
    glDeleteTextures(1, &tex_id);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);

    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }

    return stats;
}

/* There are two possibilities here. If we packed into the same number of containers
 * as the source containers the container sizes should match (note that this case
 * includes the one from a single source texture to a single destination texture).
 * If we packed into a single container, and the source containers were more than
 * one, then we should ensure that the size is smallest power of two that is
 * equal or larger than the input texture storage. */

std::vector<TextureSize> ComputeSizes(int ntex, TextureObjectHandle inputTexture)
{
    std::vector<TextureSize> textureSizes;
    int ntex_in = inputTexture->ArraySize();
    if (ntex == ntex_in) {
        for (int i = 0; i < ntex_in; ++i) {
            textureSizes.push_back({inputTexture->TextureWidth(i), inputTexture->TextureHeight(i)});
        }
        return textureSizes;
    } else {
        assert(ntex == 1);
        int64_t textureArea = 0;
        for (int i = 0; i < ntex_in; ++i) {
            textureArea += inputTexture->TextureArea(i);
        }
        if (textureArea <= 4096 * 4096)
            return { {4096, 4096} };
        else if (textureArea <= 8192 * 8192)
            return { {8192, 8192} };
        else if (textureArea <= 16384 * 16384)
            return { {16384, 16384} };
        else
            assert(0 && "Unable to find a texture large enough to store all the data");
    }
}

TextureObjectHandle RenderTexture(Mesh &m, TextureObjectHandle textureObject, bool filter, InterpolationMode imode, GLFWwindow *parentWindow)
{
    std::vector<std::vector<Mesh::FacePointer>> facesByTexture;
    int nTex = FacesByTextureIndex(m, facesByTexture);

    std::vector<TextureSize> texSizes = ComputeSizes(nTex, textureObject);
    TextureObjectHandle newTexture = std::make_shared<TextureObject>();
    for (int i = 0; i < nTex; ++i) {
        std::shared_ptr<QImage> teximg = RenderTexture(facesByTexture[i], m, textureObject, filter, imode, texSizes[i].w, texSizes[i].h, parentWindow);
        newTexture->AddImage(teximg);
    }

    return newTexture;
}

static std::shared_ptr<QImage> RenderTexture(std::vector<Mesh::FacePointer>& fvec,
                                      Mesh &m, TextureObjectHandle textureObject,
                                      bool filter, InterpolationMode imode,
                                      int textureWidth, int textureHeight, GLFWwindow *parentWindow)
{
    bool sharedContext = (parentWindow != nullptr);
    auto WTCSh = GetWedgeTexCoordStorageAttribute(m);

    // sort the faces in increasing order of input texture unit
    auto FaceComparatorByInputTexIndex = [&WTCSh](const Mesh::FacePointer& f1, const Mesh::FacePointer& f2) {
        return WTCSh[f1].tc[0].N() < WTCSh[f2].tc[0].N();
    };

    std::sort(fvec.begin(), fvec.end(), FaceComparatorByInputTexIndex);

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
    if (err) {
        cout << "glew init error " << glewGetErrorString(err) << endl;
    }

    glGetError(); // suppress possible error on glew init

    glfwSwapInterval(1);

    // OpenGL setup

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLint program = CompileShaders(vs_text, fs_text);
    glUseProgram(program);

    CheckGLError();

    // Allocate vertex data

    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    std::vector<TextureSize> inTexSizes;
    for (std::size_t i = 0; i < textureObject->ArraySize(); ++i) {
        int iw = textureObject->TextureWidth(i);
        int ih = textureObject->TextureHeight(i);
        inTexSizes.push_back({iw, ih});
    }

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, m.FN()*12*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (auto fptr : fvec) {
        int ti = WTCSh[fptr].tc[0].N();
        for (int i = 0; i < 3; ++i) {
            *p++ = fptr->cWT(i).U();
            *p++ = fptr->cWT(i).V();
            vcg::Point2d uv = WTCSh[fptr].tc[i].P();
            *p++ = uv.X() / inTexSizes[ti].w;
            *p++ = uv.Y() / inTexSizes[ti].h;
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

    int renderedTexWidth = textureWidth;
    int renderedTexHeight = textureHeight;

    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glViewport(0, 0, renderedTexWidth, renderedTexHeight);

    CheckGLError();

    GLuint renderTarget;
    glGenTextures(1, &renderTarget);
    glBindTexture(GL_TEXTURE_2D, renderTarget);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, renderedTexWidth, renderedTexHeight);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, renderTarget, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    CheckGLError();

    std::cout << "[LOG] Using interpolation mode " << imode << std::endl;

    GLint loc_cubic_flag = glGetUniformLocation(program, "use_cubic_interpolation");
    glUniform1i(loc_cubic_flag, 0);
    switch (imode) {
    case Cubic:
        glUniform1i(loc_cubic_flag, 1);
        // fall through
    case Linear:
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        break;
    case Nearest:
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    }
    CheckGLError();

    std::shared_ptr<QImage> textureImage = std::make_shared<QImage>(renderedTexWidth, renderedTexHeight, QImage::Format_ARGB32);

    // disable depth and stencil test (if they were enabled) as the render target does not have the buffers attached
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    auto f0 = fvec.begin();
    auto fbase = f0;
    while (fbase != fvec.end()) {
        auto fcurr = fbase;
        int currTexIndex = WTCSh[*fcurr].tc[0].N();
        while (fcurr != fvec.end() && WTCSh[*fcurr].tc[0].N() == currTexIndex)
            fcurr++;
        int baseIndex = std::distance(f0, fbase) * 3;
        int count = std::distance(fbase, fcurr) * 3;

        // Load texture image
        glActiveTexture(GL_TEXTURE0);
        std::cout << "Binding texture unit " << currTexIndex << std::endl;
        textureObject->Bind(currTexIndex);

        GLint loc_img0 = glGetUniformLocation(program, "img0");
        glUniform1i(loc_img0, 0);
        GLint loc_texture_size = glGetUniformLocation(program, "texture_size");
        glUniform2f(loc_texture_size, float(textureObject->TextureWidth(currTexIndex)), float(textureObject->TextureHeight(currTexIndex)));

        glDrawArrays(GL_TRIANGLES, baseIndex, count);
        CheckGLError();

        // if the context is not shared, I simply assume all the texture names are used somewhere else...
        if (sharedContext == false)
            textureObject->Release(currTexIndex);

        fbase = fcurr;
    }

    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, renderedTexWidth, renderedTexHeight, GL_BGRA, GL_UNSIGNED_BYTE, textureImage->bits());

    glfwPollEvents();

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteTextures(1, &renderTarget);
    glDeleteFramebuffers(1, &fbo);
    glDeleteBuffers(1, &vertexbuf);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);
    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }

    if (filter) vcg::PullPush(*textureImage, qRgba(0, 255, 0, 255));

    Mirror(*textureImage);
    return textureImage;
}
