#ifndef MESH_VIEWER_H
#define MESH_VIEWER_H

#include <memory>
#include <vector>
#include <string>

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

    // TODO the mesh graph (meshParamData) should be encapsulated by the GraphManager object
    // and queried through it
    std::shared_ptr<MeshGraph> meshParamData;
    std::shared_ptr<GraphManager> gm;

    std::string fileName;

    // TODO move the optimizer parameters somewhere else
    std::size_t minRegionSize;

    // GUI related variables///////////////////////////////////////// TODO move into info
    GLFWwindow *_window = nullptr;
    double _xpos = 0, _ypos = 0;
    DragMode _dragMode = DISABLED;

    float _dragX = 0.0f, _dragY = 0.0f;

    struct SelectedRegionInfo {
        std::shared_ptr<FaceGroup> chart;
        GLint first;
        GLsizei count;
        vcg::Color4f color;
    };

    enum SelectionType { None, Chart, Edge };

    SelectionType selectionType = None;
    std::vector<SelectedRegionInfo> selectionVector;

    struct {
        GLuint mesh = 0;
        GLuint selection = 0;
    } _vertexBuffers;

    TextureObjectHandle _currentTexture;

    struct {
        int colorSource = 0; // see shader code
        GLuint program = 0;
        GLuint vao = 0;
        struct {
            GLint loc_position;
            GLint loc_texcoord;
            GLint loc_distcolor;
        } attributes;
        struct {
            GLint loc_modelView;
            GLint loc_projection;
            GLint loc_colorSource;
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
    MeshViewer(std::shared_ptr<MeshGraph> meshParamData_, std::size_t minRegionSize_, const std::string &fileName_);
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
    void ClearSelection();
    void Select(const RegionID id);
    void Select(const GraphManager::Edge& e);
    void InitializeSelection(const std::vector<std::pair<RegionID,vcg::Color4f>>& vsel);

    void ManageImGuiState();

    // GraphManager interface
    void gmClose();
    bool gmHasNextEdge();
    const std::pair<GraphManager::Edge,double>& gmPeekNextEdge();
    void gmRemoveNextEdge();
    GraphManager::ChartHandle gmCollapse(const GraphManager::Edge& e);

};

#endif // MESH_VIEWER_H
