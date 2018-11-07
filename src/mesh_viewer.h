#ifndef MESH_VIEWER_H
#define MESH_VIEWER_H

#include <memory>
#include <vector>
#include <string>
#include <algorithm>

#include <GL/glew.h>

#include "mesh.h"
#include "mesh_graph.h"
#include "texture_optimization.h"
#include "linmath.h"
#include "parameterization.h"
#include "utils.h"

struct GLFWwindow;

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

    enum ColorMode { VERTEX, FACE, NONE };

    static void MouseButtonCallback(GLFWwindow *, int, int, int);
    static void CursorPositionCallback(GLFWwindow *, double, double);
    static void ScrollCallback(GLFWwindow *, double, double);
    static void KeyCallback(GLFWwindow *, int, int, int, int);
    static void FramebufferSizeCallback(GLFWwindow *, int, int);

    static constexpr int ColorMask_EMPTY      = 0;
    static constexpr int ColorMask_TEXTURE    = 1;
    static constexpr int ColorMask_PRIMITIVE  = 2;
    static constexpr int ColorMask_CHECKBOARD = 4;


private:

    GraphHandle graph;
    TextureObjectHandle _currentTexture;
    std::shared_ptr<GraphManager> gm;

    std::unordered_map<RegionID, vcg::Color4f> regionColors;

    std::string fileName;


    ParameterizationStrategy strategy;
    int nw; // number of workers for the graph parameterization


    // TODO move the optimizer parameters somewhere else
    int regionCount;

    // GUI related variables///////////////////////////////////////// TODO move into info
    GLFWwindow *_window = nullptr;
    double _xpos = 0, _ypos = 0;
    DragMode _dragMode = DISABLED;

    float _dragX = 0.0f, _dragY = 0.0f;

    // Selection related information
    std::shared_ptr<FaceGroup> shellGroup;
    std::shared_ptr<ParameterizerObject> parameterizer;
    ColorMode shellColorMode = FACE;

    struct SelectionBufferInfo {
        ChartHandle chart;
        std::size_t bufferIndex;
        GLint first;
        GLsizei count;
        GLint first_highlight;
        GLuint texicon;
        int referenceCount;
    };

    std::unordered_map<RegionID, SelectionBufferInfo> selectedRegions;
    std::unordered_map<RegionID, int> primaryCharts;  // charts selected and ready for merge

    std::vector<GLuint> selectionVao;
    //std::vector<GLuint> detailVao;
    std::vector<GLuint> highlightVao;

    struct {
        GLuint mesh = 0;
        GLuint detail = 0;
        GLuint detailBorder = 0;
        std::vector<GLuint> selection;
        std::vector<GLuint> highlight;
    } _vertexBuffers;

    struct {
        int colorMask = 1; // see  ColorMask enum and shader code
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
            GLint loc_colorMask;
            GLint loc_weight;
        } uniforms;

        struct {
            GLuint program = 0;
            //GLuint vao = 0;
            struct {
                GLint loc_position;
                GLint loc_texcoord;
                GLint loc_distcolor;
            } attributes;
            struct {
                GLint loc_modelView;
                GLint loc_projection;
                GLint loc_colorMask;
                GLint loc_weight;
            } uniforms;
        } selection;
    } _perspectiveView;

    struct {
        GLuint program = 0;
        GLuint vao = 0;
        struct {
            GLint loc_position;
            GLint loc_texcoord;
        } attributes;
        struct {
            GLint loc_projection;
        } uniforms;

        struct {
            GLuint program;
            GLuint vao;
            struct {
                GLint loc_position;
                GLint loc_texcoord;
            } attributes;
            struct {
                GLint loc_primitiveColor;
                GLint loc_projection;
            } uniforms;
        } highlight;

    } _textureView;

    struct {
        GLuint program = 0;
        GLuint vao = 0;
        GLuint borderVao = 0;
        struct {
            GLint loc_position;
            GLint loc_texcoord;
            GLint loc_color;
        } attributes;
        struct {
            GLint loc_projection;
        } uniforms;
        bool wireframe = false;
        GLsizei count = 0;
        GLsizei borderCount = 0;
        GLsizei seamCount = 0;
    } _detailView;

    struct {
        vec3 eye = {0.0f, 0.0f, 3.0f};
        vec3 target = { 0.0f, 0.0f, 0.0f};
        vec3 up = {0.0f, 1.0f, 0.0f};
        float nearClip = 0.1f;
        float farClip = 2000.0f;
    } _perspectiveCamera;

    Camera2D _textureCamera;
    Camera2D _detailCamera;

    struct {
        mat4x4 positionMatrix;
        mat4x4 trackballMatrix;
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
        int perspectiveViewport[4];
        int textureViewport[4];
        int detailViewport[4];
    } info;

    /*
    // ImGui state
    // =======================================================

    bool showInfoArea = true;

    // distortion display
    bool distortionFromTexture = false;
    DistortionMetric::Type distortion[] = {
        DistortionMetric::Area,
        DistortionMetric::Angle,
    };
    int distortionIndex = 0;
    int activeDistIndex = -1;

    // parameterization strategy
    DirectParameterizer parameterizer[] = { DCP, FixedBorderBijective };
    TexCoordOptimizer optimizer[] = { AreaPreserving, SymmetricDirichletOpt, MIPS };
    ParameterizationGeometry geometry[] = { Model, Texture };

    int parameterizerInUse = 0;
    int optimizerInUse = 0;
    int geometryInUse = 1;
    */

public:
    MeshViewer(GraphHandle gh, const Args& args);
    void Run();
    void InitBuffers();
    void UpdateDetailBuffers();
    void SetupViews();
    //void SetupDetailView(ChartHandle chart);
    void SetupDetailView(const Mesh& detailMesh);

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
    void UpdateSelection(const RegionID id);

    void ManageImGuiState();

    // GraphManager interface
    void gmClose();
    bool gmHasNextEdge();
    const std::pair<GraphManager::Edge,double>& gmPeekNextEdge();
    void gmRemoveNextEdge();
    GraphManager::ChartHandle gmCollapse(const GraphManager::Edge& e);

    template <typename ChartInputIterator>
    std::pair<int,GraphManager::ChartHandle> gmCollapse(ChartInputIterator first, ChartInputIterator last);
};

#endif // MESH_VIEWER_H
