#include "texture.h"

#include "mesh_graph.h"
#include "gl_utils.h"
#include "metric.h"
#include "mesh_attribute.h"

#include <vcg/complex/algorithms/update/color.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

static const char *vs_text_distortion[] = {
    "#version 430 core                                           \n"
    "                                                            \n"
    "in vec2 position;                                           \n"
    "in vec4 color_area;                                              \n"
    "in vec4 color_angle;                                              \n"
    "out vec4 c1;                                                 \n"
    "out vec4 c2;                                                 \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);               \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                        \n"
    "    c1 = color_area;                                              \n"
    "    c2 = color_angle;                                              \n"
    "}                                                           \n"
};

static const char *fs_text_distortion[] = {
    "#version 430 core                                                     \n"
    "                                                                      \n"
    "in vec4 c1;                                                            \n"
    "in vec4 c2;                                                            \n"
    "layout(location = 0) out vec4 color1;                                                       \n"
    "layout(location = 1) out vec4 color2;                                                       \n"
    "                                                                      \n"
    "void main(void)                                                       \n"
    "{                                                                     \n"
    "    color1 = c1;                                                        \n"
    "    color2 = c2;                                                        \n"
    "}                                                                     \n"
};

static const char *vs_text_pack[] = {
    "#version 430 core                                   \n"
    "                                                    \n"
    "in vec2 position;                                   \n"
    "in vec4 color;                                      \n"
    "out vec4 c1;                                        \n"
    "                                                    \n"
    "void main(void)                                     \n"
    "{                                                   \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);       \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                \n"
    "    c1 = color;                                     \n"
    "}                                                   \n"
};

static const char *fs_text_pack[] = {
    "#version 430 core                                  \n"
    "                                                   \n"
    "in vec4 c1;                                        \n"
    "layout(location = 0) out vec4 color1;              \n"
    "                                                   \n"
    "void main(void)                                    \n"
    "{                                                  \n"
    "    color1 = c1;                                   \n"
    "}                                                  \n"
};

void GeneratePackingQualityTexture(GraphHandle graph, TextureObjectHandle textureObject, std::unordered_map<RegionID,float>& qualityMap)
{
    assert(0 && "TODO add multitexture support");
    // Create a hidden window
    GLFWwindow *parentWindow = glfwGetCurrentContext();
    bool sharedContext = (parentWindow != nullptr);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

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

    GLint program = CompileShaders(vs_text_pack, fs_text_pack);
    glUseProgram(program);


    Mesh& m = graph->mesh;

    for (const auto& entry : graph->charts) {
        for (auto fptr : entry.second->fpVec) {
            fptr->Q() = qualityMap[entry.first];
        }
    }
    tri::UpdateColor<Mesh>::PerFaceQualityGray(m, 0, 1);
    //tri::UpdateColor<Mesh>::PerFaceQualityGray(m);

    // Allocate vertex data
    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, m.face.size()*9*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (const auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            *p++ = f.cWT(i).U();
            *p++ = f.cWT(i).V();
            vcg::Color4b color = f.cC();
            unsigned char *colorptr = (unsigned char *) p;
            *colorptr++ = color[0];
            *colorptr++ = color[1];
            *colorptr++ = color[2];
            *colorptr++ = color[3];
            p++;
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint loc_position = glGetAttribLocation(program, "position");
    glVertexAttribPointer(loc_position, 2, GL_FLOAT, GL_FALSE, 3*sizeof(float), 0);
    glEnableVertexAttribArray(loc_position);

    GLint loc_color_area = glGetAttribLocation(program, "color");
    glVertexAttribPointer(loc_color_area, 4, GL_UNSIGNED_BYTE, GL_TRUE, 3*sizeof(float), (const GLvoid *) (2*sizeof(float)));
    glEnableVertexAttribArray(loc_color_area);

    CheckGLError();

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    std::size_t width = textureObject->TextureWidth(0);
    std::size_t height = textureObject->TextureHeight(0);

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

    glClearColor(0.0f, 0.65f, 0.0f, 1.0f);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawArrays(GL_TRIANGLES, 0, m.face.size()*3);

    glfwPollEvents();

    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    QImage packabilityTexture(width, height, QImage::Format_ARGB32);
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, packabilityTexture.bits());

    std::stringstream ss;
    ss << m.name << "_" << "_packability.png";

    packabilityTexture.save(ss.str().c_str(), "png", 66);

    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteBuffers(1, &vertexbuf);
    glDeleteTextures(1, &renderTarget);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);

    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }
}









void GenerateDistortionTextures(Mesh& m, TextureObjectHandle textureObject)
{
    assert(0 && "TODO add multitexture support");
    // Create a hidden window
    GLFWwindow *parentWindow = glfwGetCurrentContext();
    bool sharedContext = (parentWindow != nullptr);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

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

    GLint program = CompileShaders(vs_text_distortion, fs_text_distortion);
    glUseProgram(program);


    // Compute distortion stats

    double texArea = textureObject->TextureWidth(0) * (double) textureObject->TextureHeight(0);
    double uvRatio = textureObject->TextureWidth(0) / (double) textureObject->TextureHeight(0);
    assert(HasWedgeTexCoordStorageAttribute(m));
    auto WTCSh = GetWedgeTexCoordStorageAttribute(m);

    double c_min = std::numeric_limits<double>::max();
    double c_max = std::numeric_limits<double>::min();
    for (auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            if (uvRatio > 1)
                f.WT(i).U() *= uvRatio;
            else if (uvRatio < 1)
                f.WT(i).V() /= uvRatio;
        }

        double q = DistortionMetric::AngleDistortion(m, f, ParameterizationGeometry::Texture);
        c_min = std::min(c_min, q);
        c_max = std::max(c_max, q);
        f.Q() = q;

        for (int i = 0; i < 3; ++i) {
            if (uvRatio > 1)
                f.WT(i).U() /= uvRatio;
            else if (uvRatio < 1)
                f.WT(i).V() *= uvRatio;
        }


        /*
        const Point2d& u0 = WTCSh[f].tc[0].P();
        const Point2d& u1 = WTCSh[f].tc[1].P();
        const Point2d& u2 = WTCSh[f].tc[2].P();
        Point2d u10 = u1 - u0;
        Point2d u20 = u2 - u0;
        double area = std::abs(u10 ^ u20) / 2.0;
        // Compute the isometry to the model coordinates and scale them to match the texture area
        Point2d x10;
        Point2d x20;
        LocalIsometry(f.P(1) - f.P(0), f.P(2) - f.P(0), x10, x20);
        double areaModel = std::abs(x10 ^ x20) / 2.0;
        x10 *= std::sqrt(area / areaModel);
        x20 *= std::sqrt(area / areaModel);
        // Evaluate the angle difference between the existing parameterization and the model
        // if the distortion is low, we use the triangle defined by the parameterization, otherwise
        // we linearly interpolate between the texture and model shape to try and correct the distortion
        double angleDist = 0;
        angleDist += std::abs(VecAngle(u10, u20) - VecAngle(x10, x20));
        angleDist += std::abs(VecAngle(u10 + (u20 - u10), -u10) - VecAngle(x10 + (x20 - x10), -x10));
        angleDist += std::abs(VecAngle(u20 - u10, -u20) - VecAngle(x20 - x10, -x20));
        double interpolationFactor = angleDist / (2.0 * M_PI);

        if (angleDist > 0.5 * M_PI)
            interpolationFactor = 1.0;
        else
            interpolationFactor = 0.2;
        */

        // compute the singular values of the transformation matrix, s2 > s1
        // ref for the formula: smith&schaefer 2015 bijective,
        /*
        Eigen::Matrix2d phi = ComputeTransformationMatrix(x10, x20, u10, u20);
        double bcplus  = std::pow(phi(0, 1) + phi(1, 0), 2.0);
        double bcminus = std::pow(phi(0, 1) - phi(1, 0), 2.0);
        double adplus  = std::pow(phi(0, 0) + phi(1, 1), 2.0);
        double adminus = std::pow(phi(0, 0) - phi(1, 1), 2.0);
        double s1 = 0.5 * std::abs(std::sqrt(bcplus + adminus) - std::sqrt(bcminus + adplus));
        double s2 = 0.5 * (std::sqrt(bcplus + adminus) + std::sqrt(bcminus + adplus));

        interpolationFactor = 1.0 - (s1 / s2);
        assert(interpolationFactor > 0);
        assert(interpolationFactor <= 1);


        f.Q() = interpolationFactor;
        */
    }

    tri::UpdateColor<Mesh>::PerFaceQualityGray(m);
    std::vector<vcg::Color4b> color_angle;
    for (auto& f : m.face)
        color_angle.push_back(f.C());

    std::cout << "[LOG] Texture angle distortion range: (" << c_min << ", " << c_max << ")" << std::endl;

    float qa_min = std::numeric_limits<float>::max();
    float qa_max = std::numeric_limits<float>::min();

    texArea = 1;
    for (std::size_t i = 0; i < m.face.size(); ++i) {
        auto &f = m.face[i];
        vcg::Point2d vuv[3];
        for (int i = 0; i < 3; ++i) {
            vcg::Point2d uv = WTCSh[f].tc[i].P();
            if (uvRatio > 1)
                uv.X() /= uvRatio;
            else if (uvRatio < 1)
                uv.Y() *= uvRatio;
            vuv[i] = uv;
        }
        double initialArea = (0.5 * ((vuv[1] - vuv[0]) ^ (vuv[2] - vuv[0]))) * texArea;
        double finalArea = (0.5 * (f.cWT(1).P() - f.cWT(0).P()) ^ (f.cWT(2).P() - f.cWT(0).P())) * texArea;
        double qa = (finalArea - initialArea) / (0.5 * (f.cWT(1).P() - f.cWT(0).P()) ^ (f.cWT(2).P() - f.cWT(0).P()));
        if (qa < qa_min) qa_min = qa;
        if (qa > qa_max) qa_max = qa;

        f.Q() = qa;
    }

    std::cout << "[LOG] Texture area distortion range: (" << qa_min << ", " << qa_max << ")" << std::endl;

    for (auto& f : m.face) {
        float q = f.Q();
        if (q < 0) {
            float v = 1.0f + (q / (-qa_min));
            f.C().Import(Color4f{1.0f, v, v, 1.0f});
        } else {
            //float v = 1.0f - (q / range.second);
            float v = 1.0f - (q / qa_max);
            f.C().Import(Color4f{v, v, 1.0f, 1.0f});
        }
    }

    //tri::UpdateColor<Mesh>::PerFaceQualityRamp(m);

    // Allocate vertex data
    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, m.face.size()*12*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (const auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            *p++ = WTCSh[f].tc[i].U() / ((uvRatio > 1) ? uvRatio : 1.0);
            *p++ = WTCSh[f].tc[i].V() * ((uvRatio < 1) ? uvRatio : 1.0);
            //*p++ = f.cWT(i).U() / ((uvRatio > 1) ? uvRatio : 1.0);
            //*p++ = f.cWT(i).V() * ((uvRatio < 1) ? uvRatio : 1.0);
            vcg::Color4b color = f.cC();
            vcg::Color4b angleColor = color_angle[tri::Index(m, f)];
            unsigned char *colorptr = (unsigned char *) p;
            *colorptr++ = color[0];
            *colorptr++ = color[1];
            *colorptr++ = color[2];
            *colorptr++ = color[3];
            *colorptr++ = angleColor[0];
            *colorptr++ = angleColor[1];
            *colorptr++ = angleColor[2];
            *colorptr++ = angleColor[3];
            p++;
            p++;
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint loc_position = glGetAttribLocation(program, "position");
    glVertexAttribPointer(loc_position, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), 0);
    glEnableVertexAttribArray(loc_position);

    GLint loc_color_area = glGetAttribLocation(program, "color_area");
    glVertexAttribPointer(loc_color_area, 4, GL_UNSIGNED_BYTE, GL_TRUE, 4*sizeof(float), (const GLvoid *) (2*sizeof(float)));
    glEnableVertexAttribArray(loc_color_area);

    GLint loc_color_angle = glGetAttribLocation(program, "color_angle");
    glVertexAttribPointer(loc_color_angle, 4, GL_UNSIGNED_BYTE, GL_TRUE, 4*sizeof(float), (const GLvoid *) (3*sizeof(float)));
    glEnableVertexAttribArray(loc_color_angle);

    CheckGLError();

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    std::size_t width = textureObject->TextureWidth(0);
    std::size_t height = textureObject->TextureHeight(0);

    GLuint renderTarget[2];
    glGenTextures(2, renderTarget);

    GLenum drawBuffer[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
    for (int i = 0; i < 2; ++i) {
        glBindTexture(GL_TEXTURE_2D, renderTarget[i]);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, width, height);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture(GL_FRAMEBUFFER, drawBuffer[i], renderTarget[i], 0);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glClearColor(0.0f, 0.2f, 0.0f, 1.0f);
    glDrawBuffers(2, drawBuffer);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawArrays(GL_TRIANGLES, 0, m.face.size()*3);

    glfwPollEvents();

    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    QImage distortionMap(width, height, QImage::Format_ARGB32);
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, distortionMap.bits());

    std::stringstream ss;
    ss << m.name << "_" << qa_min << "_" << qa_max << "_texture_area_distortion.png";

    distortionMap.save(ss.str().c_str(), "png", 66);

    glReadBuffer(GL_COLOR_ATTACHMENT1);
    glReadPixels(0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, distortionMap.bits());
    ss.str("");
    ss.clear();
    ss << m.name << "_" << c_min << "_" << c_max << "_texture_angle_distortion.png";
    distortionMap.save(ss.str().c_str(), "png", 66);


    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteBuffers(1, &vertexbuf);
    glDeleteTextures(2, renderTarget);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);

    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }
}

static const char *vs_text_gutter[] = {
    "#version 410 core                                           \n"
    "                                                            \n"
    "in vec2 position;                                           \n"
    "in vec3 coord;                                              \n"
    "out vec2 vPos;                                              \n"
    "out vec3 vCoord;                                            \n"
    "                                                            \n"
    "void main(void)                                             \n"
    "{                                                           \n"
    "    vec2 p = 2.0 * position - vec2(1.0, 1.0);               \n"
    "    gl_Position = vec4(p, 0.5, 1.0);                        \n"
    "    vCoord = coord;                                         \n"
    "    vPos = position;                                        \n"
    "}                                                           \n"
};

static const char *fs_text_gutter_pass1[] = {
    "#version 410 core                                                     \n"
    "                                                                      \n"
    "in vec3 vCoord;                                                       \n"
    "out vec4 color;                                                       \n"
    "                                                                      \n"
    "void main(void)                                                       \n"
    "{                                                                     \n"
    "    color = vec4(vCoord, 1.0);                                        \n"
    "}                                                                     \n"
};

static const char *fs_text_gutter_pass2[] = {
    "#version 410 core                                                     \n"
    "                                                                      \n"
    "uniform vec2 texture_size;                                            \n"
    "uniform sampler2D geometryTexture;                                    \n"
    "                                                                      \n"
    "in vec3 vCoord;                                                       \n"
    "in vec2 vPos;                                                         \n"
    "out vec4 color2;                                                      \n"
    "                                                                      \n"
    "void main(void)                                                       \n"
    "{                                                                     \n"
    "    vec3 p = texture2D(geometryTexture, vPos);                        \n"
    "    color2 = vec4(distance(p, vPos), 1.0);                            \n"
    "}                                                                     \n"
};

/*
 * TODO
 * L'idea e' che mi creo una geometry texture, e gli do una passata di push-pull
 * Dopodiche' faccio un nuovo rendering della texture dove misuro la differenza
 * tra il lookup della posizione nella geometry texture e la posizione interpolata
 * In teoria piu' questa differenza e' grande, piu' e' probabile che il lookup
 * nella geometry texture abbia interpolato valori appartenenti a due patch diverse
 * Il rendering lo faccio come immagine di intensita' normalizzata (magari con un
 * gradiente non lineare se occorre, per favorire la visualizzazione).
 *
 * Bisogna vedere se il push-pull funziona senza colori
 */

void EvaluateGutterResistance(Mesh& m, TextureObjectHandle textureObject)
{
#if false
    // Create a hidden window
    GLFWwindow *parentWindow = glfwGetCurrentContext();
    bool sharedContext = (parentWindow != nullptr);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

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

    GLint program1 = CompileShaders(vs_text_gutter, fs_text_gutter_pass1);
    GLint program2 = CompileShaders(vs_text_gutter, fs_text_gutter_pass2);

    glUseProgram(program1);

    // Allocate vertex data
    GLuint vertexbuf;
    glGenBuffers(1, &vertexbuf);

    glBindBuffer(GL_ARRAY_BUFFER, vertexbuf);
    glBufferData(GL_ARRAY_BUFFER, m.face.size()*15*sizeof(float), NULL, GL_STATIC_DRAW);
    float *p = (float *) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (const auto& f : m.face) {
        for (int i = 0; i < 3; ++i) {
            *p++ = f.cWT(i).U();
            *p++ = f.cWT(i).V();
            *p++ = f.cP(i)[0];
            *p++ = f.cP(i)[1];
            *p++ = f.cP(i)[2];
        }
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLint loc_position = glGetAttribLocation(program1, "position");
    glVertexAttribPointer(loc_position, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
    glEnableVertexAttribArray(loc_position);

    GLint loc_color_area = glGetAttribLocation(program1, "coord");
    glVertexAttribPointer(loc_color_area, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), (const GLvoid *) (2*sizeof(float)));
    glEnableVertexAttribArray(loc_color_area);

    CheckGLError();

    p = nullptr;
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    std::size_t width = textureObject->TextureWidth(0);
    std::size_t height = textureObject->TextureHeight(0);

    GLuint renderTarget[2];
    glGenTextures(2, renderTarget);

    for (int i = 0; i < 2; ++i) {
        glBindTexture(GL_TEXTURE_2D, renderTarget[i]);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, renderTarget[0], 0);

    glViewport(0, 0, width, height);
    glScissor(0, 0, width, height);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_STENCIL_TEST);

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawArrays(GL_TRIANGLES, 0, m.face.size()*3);

    glfwPollEvents();

    glMemoryBarrier(GL_ALL_BARRIER_BITS);


    // Pass 2

    /*
    QImage distortionMap(width, height, QImage::Format_ARGB32);
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, distortionMap.bits());

    std::stringstream ss;
    ss << m.name << "_" << qa_min << "_" << qa_max << "_texture_area_distortion.png";

    distortionMap.save(ss.str().c_str(), "png", 66);

    glReadBuffer(GL_COLOR_ATTACHMENT1);
    glReadPixels(0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, distortionMap.bits());
    ss.str("");
    ss.clear();
    ss << m.name << "_" << c_min << "_" << c_max << "_texture_angle_distortion.png";
    distortionMap.save(ss.str().c_str(), "png", 66);
    */


    // clean up
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);

    glDeleteBuffers(1, &vertexbuf);
    glDeleteTextures(2, renderTarget);
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(window);

    if (sharedContext) {
        glfwMakeContextCurrent(parentWindow);
    }
#endif
}

