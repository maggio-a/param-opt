#include "mesh.h"
#include "mesh_graph.h"
#include "uv.h"
#include "texture_optimization.h"
#include "texture_rendering.h"
#include "timer.h"
#include "gl_utils.h"
#include "texture.h"
#include "mesh_viewer.h"
#include "logging.h"
#include "utils.h"

#include <wrap/io_trimesh/io_mask.h>

#include <string>
#include <vector>
#include <iostream>

#include <QCoreApplication>
#include <QImage>
#include <QDir>
#include <QFileInfo>
#include <QString>

#include <GL/glew.h>
#include <GLFW/glfw3.h>


static const char *vs_text[] = {
    "#version 410 core                                                             \n"
    "                                                                              \n"
    "uniform mat4 modelViewMatrix;                                                 \n"
    "uniform mat4 projectionMatrix;                                                \n"
    "                                                                              \n"
    "in vec3 position;                                                             \n"
    "in vec2 texcoord;                                                             \n"
    "out vec2 uv;                                                                  \n"
    "                                                                              \n"
    "void main(void)                                                               \n"
    "{                                                                             \n"
    "    uv = texcoord;                                                            \n"
    "    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0f);  \n"
    "}                                                                             \n"
};

static const char *fs_text[] = {
    "#version 410 core                                              \n"
    "                                                               \n"
    "uniform sampler2D tex0;                                        \n"
    "                                                               \n"
    "uniform int zero_mask = 0;                                     \n"
    "                                                               \n"
    "in vec2 uv;                                                    \n"
    "out vec4 color;                                                \n"
    "                                                               \n"
    "void main(void)                                                \n"
    "{                                                              \n"
    "    if (zero_mask == 0) {                                      \n"
    "        color = texture(tex0, uv);                             \n"
    "    } else {                                                   \n"
    "        if (uv == vec2(0.0, 0.0))                              \n"
    "            color = vec4(1.0, 1.0, 1.0, 1.0);                  \n"
    "        else                                                   \n"
    "            color = vec4(0.0, 0.0, 0.0, 1.0);                  \n"
    "    }                                                          \n"
    "}                                                              \n"
};


using namespace vcg;

void RenderStills(Mesh &m, TextureObjectHandle textureObject, bool halfres, float lodbias, float max_aniso, bool mask)
{
    tri::UpdateBounding<Mesh>::Box(m);

    bool scaled = false;
    double scaleFactor = 0;

    if (m.bbox.Diag() > 100) {
        scaled = true;
        scaleFactor = 1.0 / m.bbox.Diag();
        tri::UpdatePosition<Mesh>::Scale(m, scaleFactor);
        tri::UpdateBounding<Mesh>::Box(m);
    }

    int w = 2048;
    int h = 2048;
    if (halfres) {
        w /= 2;
        h /= 2;
    }
    float aspect = w / (float) h;

    int ntex = textureObject->ArraySize();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    GLFWwindow *window = glfwCreateWindow(512, 512, "Window", NULL, NULL);
    if (!window)
    {
        LOG_ERR << "Failed to create window or context";
    }

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err) {
        LOG_ERR << "glew init error " << glewGetErrorString(err);
    }

    glGetError(); // suppress possible error on glew init

    glfwSwapInterval(1);

    GLint program = CompileShaders(vs_text, fs_text);
    glUseProgram(program);

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Allocate vertex data

    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, m.FN()*15*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

    /* iterate ntex times over the mesh, each time recording in the buffer the
     * faces that refer to a different texture unit. The facesPerTexUnit vector
     * stores at index i the number of primitives that must be rendered when
     * texture i is bound. This avoids loading all the textures at once and
     * running out of gpu memory if the mesh has many high resolution textues */

    std::vector<unsigned> facesPerTexUnit;


    for (int i = 0; i < ntex; ++i) {
        unsigned nfi = 0;
        for (const auto& f : m.face) {
            int ti = f.cWT(0).N();
            if (ti == i) {
                for (int j = 0; j < 3; ++j) {
                    *p++ = f.cP(j).X();
                    *p++ = f.cP(j).Y();
                    *p++ = f.cP(j).Z();
                    *p++ = f.cWT(j).U();
                    *p++ = f.cWT(j).V();
                }
                nfi++;
            }
        }
        facesPerTexUnit.push_back(nfi);
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint pos_location = glGetAttribLocation(program, "position");
    glVertexAttribPointer(pos_location, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
    glEnableVertexAttribArray(pos_location);

    GLint tc_location = glGetAttribLocation(program, "texcoord");
    glVertexAttribPointer(tc_location, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void *) (3*sizeof(float)));
    glEnableVertexAttribArray(tc_location);

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0); // done, unbind

    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glViewport(0, 0, w, h);

    GLuint colorTexture;
    glGenTextures(1, &colorTexture);
    glBindTexture(GL_TEXTURE_2D, colorTexture);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, w, h);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, colorTexture, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    GLuint depthTexture;
    glGenTextures(1, &depthTexture);
    glBindTexture(GL_TEXTURE_2D, depthTexture);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_DEPTH_COMPONENT32F, w, h);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depthTexture, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    CheckGLError();


    /* To compute the viewpoints distance we take the radius of the bounding
     * sphere incremented by a small amount */

    vcg::Point3d center = m.bbox.Center();
    double radius = 0;
    for (auto& v : m.vert)
        radius = std::max(radius, vcg::Distance(center, v.P()));

    mat4x4 model;
    {
        vcg::Matrix44f translate;
        translate.SetTranslate(-center.X(), -center.Y(), -center.Z());

        /*
        vcg::Matrix44f scale;
        scale.SetScale(s, s, s);

        vcg::Matrix44f mod = scale * translate;
        */
        vcg::Matrix44f mod = translate;
        mod.transposeInPlace();

        memcpy(&model, mod.V(), 16*sizeof(float));
    }

    mat4x4 proj;
    mat4x4_perspective(proj, 45.0f * M_PI / 180.0f, aspect, 0.001f, 2000.0f);
    //mat4x4_perspective(proj, 0.1f * M_PI / 180.0f, aspect, 0.001f, 2000.0f);

    // projection matrix is fixed
    GLint loc_projectionMatrix = glGetUniformLocation(program, "projectionMatrix");
    glUniformMatrix4fv(loc_projectionMatrix, 1, GL_FALSE, (const GLfloat *) proj);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_SCISSOR_TEST);

    Mesh viewPoints;
    tri::Icosahedron(viewPoints);
    tri::UpdatePosition<Mesh>::Scale(viewPoints, radius);
    int n = 0;
    for (const auto& vp : viewPoints.vert) {

        LOG_DEBUG << "Rendering image " << n;

        Point3f viewPoint;
        viewPoint.Import(vp.cP());

        float target[] = { 0.0f, 0.0f, 0.0f };
        float up[] = { 0.0f, 1.0f, 0.0f };
        float eye[] = { viewPoint.X(), viewPoint.Y(), viewPoint.Z() };
        //float eye[] = { 0.0f, 0.0f, 3.0f };

        mat4x4 view;
        mat4x4_identity(view);
        mat4x4_look_at(view, eye, target, up);

        mat4x4 modelView;
        mat4x4_mul(modelView, view, model);

        GLint loc_modelViewMatrix = glGetUniformLocation(program, "modelViewMatrix");
        glUniformMatrix4fv(loc_modelViewMatrix, 1, GL_FALSE, (const GLfloat *) modelView);


        // render in multiple passes

        if (mask == false)
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        else
            glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        int base = 0;
        for (int i = 0; i < ntex; ++i) {
            int vertexCount = facesPerTexUnit[i] * 3;

            // Load texture image
            glActiveTexture(GL_TEXTURE0);
            textureObject->Bind(i);

            // set filtering
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, lodbias);
            glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, max_aniso);

            CheckGLError();

            GLint loc_tex0 = glGetUniformLocation(program, "tex0");
            glUniform1i(loc_tex0, 0);

            GLint loc_zero_mask = glGetUniformLocation(program, "zero_mask");
            glUniform1i(loc_zero_mask, (mask ? 1 : 0));

            glDrawArrays(GL_TRIANGLES, base, vertexCount);
            CheckGLError();
            textureObject->Release(i);

            base += vertexCount;
        }


        glMemoryBarrier(GL_ALL_BARRIER_BITS);

        QImage img(w, h, QImage::Format_ARGB32);
        glReadBuffer(GL_COLOR_ATTACHMENT0);
        glReadPixels(0, 0, w, h, GL_BGRA, GL_UNSIGNED_BYTE, img.bits());

        Mirror(img);

        char name[256] = {};
        std::snprintf(name, 255, "img%3d.png", n++);
        img.save(name, "png", 66);
        glfwPollEvents();
    }


    glfwPollEvents();

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteTextures(1, &colorTexture);
    glDeleteTextures(1, &depthTexture);
    glDeleteFramebuffers(1, &fbo);
    glDeleteBuffers(1, &vertexbuf);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);


    if (scaled) {
        tri::UpdatePosition<Mesh>::Scale(m, (1.0 / scaleFactor));
        tri::UpdateBounding<Mesh>::Box(m);
    }
}


int main(int argc, char *argv[])
{
    // Make sure the executable directory is added to Qt's library path
    {
        QFileInfo executableInfo(argv[0]);
        QCoreApplication::addLibraryPath(executableInfo.dir().absolutePath());
    }

    ensure_condition(argc > 1 && "Mesh argument missing");

    bool downsample1 = false;
    bool downsample2 = false;
    bool halfres = false;
    bool mask = false;

    float lodbias = 0.0f;
    float max_aniso = 1.0f;

    for (int i = 2; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--downsample1") {
            downsample1 = true;
        } else if (arg == "--downsample2") {
            downsample2 = true;
        } else if (arg == "--halfres") {
            halfres = true;
        } else if (arg == "--mask") {
            mask = true;
        } else if (arg == "--lodbias") {
            i++;
            ensure_condition(i < argc);
            lodbias = std::stof(std::string(argv[i]));
        } else if (arg == "--max_anisotropy") {
            i++;
            ensure_condition(i < argc);
            max_aniso = std::stof(std::string(argv[i]));
        } else {
            LOG_WARN << "Unrecognized optional argument '" << arg << "'";
        }
    }

    Mesh m;
    TextureObjectHandle textureObject;
    int loadMask;
    if (LoadMesh(m, argv[1], textureObject, loadMask) == false) {
        LOG_ERR << "Failed to open mesh";
        std::exit(-1);
    }

    tri::Clean<Mesh>::RemoveUnreferencedVertex(m);
    tri::Allocator<Mesh>::CompactEveryVector(m);

    for (auto& f : m.face) {
        if (DistortionMetric::AreaUV(f) == 0) {
            f.WT(0).P() = vcg::Point2d(0, 0);
            f.WT(1).P() = vcg::Point2d(0, 0);
            f.WT(2).P() = vcg::Point2d(0, 0);
        }
    }

    if (downsample1 || downsample2) {
        LOG_INFO << "Downsampling textures...";
        for (auto &qimg : textureObject->imgVec) {
            std::shared_ptr<QImage> qimgdown = std::make_shared<QImage>();
            int div = 0;
            if (downsample1)
                div = 2;
            if (downsample2)
                div = 4;
            assert(div > 0);
            *qimgdown = qimg->scaledToHeight(qimg->height() / div, Qt::SmoothTransformation);
            *qimg = qimgdown->convertToFormat(QImage::Format_ARGB32); // required because scaled() can change the internal format of the QImage...
        }
    }

    ensure_condition(loadMask & tri::io::Mask::IOM_WEDGTEXCOORD);

    GLInit();

    LOG_INFO << "Rendering images...";
    RenderStills(m, textureObject, halfres, lodbias, max_aniso, mask);

    GLTerminate();

    return 0;
}
