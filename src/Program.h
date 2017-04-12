#pragma once

#include "Camera.h"
#include "Utils.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <AntTweakBar.h>
#include <PolyVoxCore/MarchingCubesSurfaceExtractor.h>
#include <PolyVoxCore/SurfaceMesh.h>
#include <PolyVoxCore/SimpleVolume.h>
#include <vector>

// Callback function called by AntTweakBar when the single-step button is clicked.
void TW_CALL stepButtonCallback(void* clientData);

// Callback function called by AntTweakBar when the particle reset button is clicked.
void TW_CALL particleResetButtonCallback(void* clientData);

// Returns the density value above which a voxel is considered solid.
// Template specialization of PolyVox::DefaultMarchingCubesController<float>::getThreshold()
template<> inline PolyVox::DefaultMarchingCubesController<float>::DensityType
    PolyVox::DefaultMarchingCubesController<float>::getThreshold()
{
    return 10.0f;
}

struct GridCellNeighborhood
{
    int minX, maxX;
    int minY, maxY;
    int minZ, maxZ;
};

class Program
{
public:
    Program(GLFWwindow* window);
    ~Program();

    // Updates the state of the program.
    void update();
    // Draws a new frame.
    void draw();
    // Resets all particles to their initial position, with zero velocity.
    void resetParticles();

    // Called when the mouse cursor is moved.
    void onMouseMoved(float dxPos, float dyPos);
    // Called when the mouse wheel is scrolled.
    void onMouseScrolled(float yOffset);
    // Called when a key is pressed or released.
    void onKeypress(int key, int action);

    // Delta time step; the simulation advances exactly dt seconds every update() call. 
    float dt = 0.01f;
    // Is the simulation paused?
    bool paused = true;

    bool holdForward = false;
    bool holdBackward = false;
    bool holdRight = false;
    bool holdLeft = false;
    bool holdUp = false;
    bool holdDown = false;
    bool holdShift = false;

private:
    Camera camera;
    GLFWwindow* window;
    TwBar* antTweakBar;

    GLuint meshVAO, pointsVAO;
    GLuint meshVBO, pointsVBO;
    GLuint meshEBO;
    GLuint simpleVertexShader, simpleFragmentShader, simpleShaderProgram;
    GLuint waterVertexShader, waterFragmentShader, waterShaderProgram;
    GLint simpleViewUniform, waterViewUniform;
    GLint waterCamUniform;
    int windowSizeX;
    int windowSizeY;

    /* DEBUG */
    // The particle positions array ('r'), for now, contains (x, y, z) coordinates for a cube with sides of size cubeSize
    static const int cubeSize = 9;
    static const int particleCount = cubeSize * cubeSize * cubeSize;
    /* END DEBUG */

    // All particles reside inside a box with these dimensions
    static constexpr float minPosX = -2.0f;
    static constexpr float minPosY = -2.0f;
    static constexpr float minPosZ = -2.0f;
    static constexpr float maxPosX =  2.0f;
    static constexpr float maxPosY =  2.0f;
    static constexpr float maxPosZ =  2.0f;
    glm::vec3 minPos{minPosX, minPosY, minPosZ};
    glm::vec3 maxPos{maxPosX, maxPosY, maxPosZ};
    glm::vec3 gridOffset = { 0.0f, 0.0f, 0.0f };

    glm::vec3 worldBoundsVertices[16] = {
        {minPos.x, minPos.y, minPos.z},
        {maxPos.x, minPos.y, minPos.z},
        {maxPos.x, maxPos.y, minPos.z},
        {minPos.x, maxPos.y, minPos.z},
        {minPos.x, minPos.y, minPos.z},
        {minPos.x, minPos.y, maxPos.z},
        {maxPos.x, minPos.y, maxPos.z},
        {maxPos.x, maxPos.y, maxPos.z},
        {minPos.x, maxPos.y, maxPos.z},
        {minPos.x, minPos.y, maxPos.z},
        {minPos.x, maxPos.y, maxPos.z},
        {minPos.x, maxPos.y, minPos.z},
        {maxPos.x, maxPos.y, minPos.z},
        {maxPos.x, maxPos.y, maxPos.z},
        {maxPos.x, minPos.y, maxPos.z},
        {maxPos.x, minPos.y, minPos.z},
    };

    // Radius of influence
    static constexpr float h = 0.5f;
    // Gas constant
    float k = 1000.0f;
    // Rest density
    float rho0 = 20.0f;
    // Mass of each particle
    float m = 1.0f;
    // Fluid viscosity
    float mu = 3.0f;
    // Surface tension coefficient
    float sigma = 0.01f;
    // Surface tension is only evaluated if |n| exceeds this threshold
    // (where n is the gradient field of the smoothed color field).
    float csNormThreshold = 1.0f;
    // Gravity acceleration
    glm::vec3 gravity{0.0f, -10.0f, 0.0f};

    // Particle positions
    glm::vec3 r[particleCount];
    // Particle velocities
    glm::vec3 v[particleCount];

    // Particle grid data structure: changes O(n^2) to O(nm): we don't have to
    // check all other particles, but only particles in adjacent grid cells.
    // Possible TODO: Store entire particle data for better cache usage
    static constexpr int gridSizeX = ceiling((maxPosX - minPosX) / h);
    static constexpr int gridSizeY = ceiling((maxPosY - minPosY) / h);
    static constexpr int gridSizeZ = ceiling((maxPosZ - minPosZ) / h);
    std::vector<size_t> particleGrid[gridSizeX][gridSizeY][gridSizeZ];

    void fillParticleGrid();

    glm::vec3 calcPressureForce(size_t particleId, float* rho, float* p);
    glm::vec3 calcViscosityForce(size_t particleId, float* rho);
    glm::vec3 calcSurfaceForce(size_t particleId, float* rho);

    static constexpr int lookupTableSize = (int)(h*h * 1e4f); // e.g. (0.5 * 0.5) * 1e4 = 2500

    float poly6LookupTable[lookupTableSize];

    void fillKernelLookupTables();

    float calcDensity(size_t particleId) const;
    float calcDensity(const glm::vec3& pos) const;

    GridCellNeighborhood getAdjacentCells(const glm::vec3& pos) const;

    PolyVox::Vector3DInt32 worldPosToVoxelIndex(const glm::vec3& worldPos) const;
    void fillVoxelVolume();
    // World positions will be multiplied by this scale for the purposes of voxel indexing.
    // Example: if this scale is 10, a world pos of (-1, 0, 2.5) will map to the voxel at (-10, 0, 25).
    // Higher values result in a voxel grid of a higher resolution; individual voxels would be smaller.
    float voxelVolumeResolutionScale = 10.0f;
    PolyVox::SimpleVolume<float> voxelVolume{{
            worldPosToVoxelIndex(minPos) - PolyVox::Vector3DInt32{
                (int)std::ceil(h * voxelVolumeResolutionScale),
                (int)std::ceil(h * voxelVolumeResolutionScale),
                (int)std::ceil(h * voxelVolumeResolutionScale)
            },
            worldPosToVoxelIndex(maxPos) + PolyVox::Vector3DInt32{
                (int)std::ceil(h * voxelVolumeResolutionScale),
                (int)std::ceil(h * voxelVolumeResolutionScale),
                (int)std::ceil(h * voxelVolumeResolutionScale)
            },
        }};
    PolyVox::SurfaceMesh<PolyVox::PositionMaterialNormal> surfaceMesh;
    PolyVox::MarchingCubesSurfaceExtractor<PolyVox::SimpleVolume<float>> surfaceExtractor{
        &voxelVolume, voxelVolume.getEnclosingRegion(), &surfaceMesh
    };

};
