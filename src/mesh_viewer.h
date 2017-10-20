#ifndef MESH_VIEWER_H
#define MESH_VIEWER_H

#include <memory>
#include <GL/glew.h>

#include "mesh_graph.h"
#include "linmath.h"

class GLFWwindow;

class MeshViewer {

public:
    struct Camera2D {
        // center
        float x;
        float y;
        // dimension of the view rectangle
        float viewSize;

        Camera2D() : x{0.5f}, y{0.5f}, viewSize{1.0f} {}

        void Reset() {
            x = y = 0.5f; viewSize = 1.0f;
        }

        void MoveX(float dx) {
            x += dx * viewSize;
            ClipView();
        }

        void MoveY(float dy) {
            y += dy * viewSize;
            ClipView();
        }

        void ClipView() {
            x += std::max(0.0f - (x - viewSize/2.0f), 0.0f);
            x -= std::max((x + viewSize/2.0f) - 1.0f, 0.0f);
            y += std::max(0.0f - (y - viewSize/2.0f), 0.0f);
            y -= std::max((y + viewSize/2.0f) - 1.0f, 0.0f);
        }

        void ZoomIn(float mouse_offset_x, float mouse_offset_y) {
            float delta = (viewSize - 0.9f*viewSize)/2.0f;
            x += delta*2*mouse_offset_x;
            y += delta*2*mouse_offset_y;
            viewSize *= 0.9f;
            //ClipView();
        }

        void ZoomOut(float mouse_offset_x, float mouse_offset_y) {
            float delta = (viewSize - viewSize*1.1f)/2.0f;
            x += delta*2*mouse_offset_x;
            y += delta*2*mouse_offset_y;
            viewSize = std::min(viewSize*1.1f, 1.0f);
            ClipView();
        }
    };

    enum DragMode { DISABLED, PERSPECTIVE, TEXTURE, DETAIL };

    static void MouseButtonCallback(GLFWwindow *, int, int, int);
    static void CursorPositionCallback(GLFWwindow *, double, double);
    static void ScrollCallback(GLFWwindow *, double, double);
    static void KeyCallback(GLFWwindow *, int, int, int, int);
    static void FramebufferSizeCallback(GLFWwindow *, int, int);

private:
    std::shared_ptr<MeshGraph> meshParamData;

    // GUI related variables///////////////////////////////////////// TODO move into info
    GLFWwindow *_window = nullptr;
    double _xpos = 0, _ypos = 0;
    DragMode _dragMode = DISABLED;

    std::shared_ptr<FaceGroup> _selected;

    float _dragX = 0.0f, _dragY = 0.0f;

    struct {
        GLuint mesh = 0;
        GLuint selection = 0;
    } _vertexBuffers;

    GLuint _texture = 0;

    struct {
        GLuint program = 0;
        GLuint vao = 0;
        struct {
            GLint loc_position;
            GLint loc_texcoord;
        } attributes;
        struct {
            GLint loc_modelView;
            GLint loc_projection;
            GLint loc_weight;
        } uniforms;

        struct {
            GLuint program = 0;
            GLuint vao = 0;
            struct {
                GLint loc_position;
                GLint loc_texcoord;
            } attributes;
            struct {
                GLint loc_modelView;
                GLint loc_projection;
                GLint loc_weight;
            } uniforms;
        } selection;
    } _perspectiveView;

    struct {
        GLuint program = 0;
        GLuint vao = 0;
        struct {
            GLint loc_texcoord;
        } attributes;
        struct {
            GLint loc_projection;
        } uniforms;
    } _textureView;

    struct {
        GLuint program = 0;
        GLuint vao = 0;
        struct {
            GLint loc_texcoord;
        } attributes;
        struct {
            GLint loc_projection;
        } uniforms;
        bool wireframe = false;
    } _detailView;

    struct {
        vec3 eye = {0.0f, 0.0f, 15.0f};
        vec3 target = { 0.0f, 0.0f, 0.0f};
        vec3 up = {0.0f, 1.0f, 0.0f};
        float near = 0.1f;
        float far = 2000.0f;
    } _perspectiveCamera;

    Camera2D _textureCamera;
    Camera2D _detailCamera;

    struct {
        mat4x4 positionMatrix;
        mat4x4 orientationMatrix;
        mat4x4 scaleMatrix;
        mat4x4 viewMatrix;
        mat4x4 projectionMatrix;
    } _meshTransform;

    struct {
        mat4x4 modelMatrix;
    } _detailTransform;

    struct {
        int width;
        int height;
        int xSplit;
        float perspectiveViewAspect;
    } info;

public:
    MeshViewer(std::shared_ptr<MeshGraph> meshParamData_);
    void Run();
    GLuint CompileShaders(const GLchar **vs_text, const GLchar **fs_text);
    void InitBuffers();
    void SetupViews();
    void SetupDetailView(const RegionID id);

    void UpdateTransforms();
    void DrawViews();
    void Draw3DView();
    void DrawTextureView();
    void DrawDetailView();

    bool InPerspectiveView();
    bool InTextureView();
    bool InDetailView();
    bool IntersectionMouseRayModel(Mesh::ConstFacePointer *fp, float &u, float &v);
    void CenterPerspectiveViewFromMouse();
    void PickRegion();
    void PerspectivePick();
    void TexturePick();
    void InitializeSelection(const RegionID id);
    void ClearSelection();
};


#endif // MESH_VIEWER_H

