#pragma once

#include "FluidBase.h"
#include "Utils.h"
#include <vector>

class FluidSimulator : public FluidBase
{
public:
    explicit FluidSimulator(GLFWwindow* window);

    // Updates the state of the program.
    void update();
    // Resets all particles to their initial position, with zero velocity.
    void resetParticles();

    // Delta time step; the simulation advances exactly dt seconds every update() call. 
    float dt = 0.01f;
    // Is the simulation paused?
    bool paused = true;

protected:
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
    // Particle positions in the previous frame
    glm::vec3 rPrev[particleCount];

    float calcDensity(size_t particleId) const;
    float calcDensity(const glm::vec3& pos) const;

    glm::vec3 calcPressureForce(size_t particleId, const float* const rho, const float* const p) const;
    glm::vec3 calcViscosityForce(size_t particleId, const float* const rho) const;
    glm::vec3 calcSurfaceForce(size_t particleId, const float* const rho) const;

    static constexpr int lookupTableSize = (int)(h*h * 1e4f); // e.g. (0.5 * 0.5) * 1e4 = 2500
    float poly6LookupTable[lookupTableSize];

    void fillKernelLookupTables();

    // Particle grid data structure: changes O(n^2) to O(nm): we don't have to
    // check all other particles, but only particles in adjacent grid cells.
    // Possible TODO: Store entire particle data for better cache usage
    static constexpr int gridSizeX = ceiling((maxPosX - minPosX) / h);
    static constexpr int gridSizeY = ceiling((maxPosY - minPosY) / h);
    static constexpr int gridSizeZ = ceiling((maxPosZ - minPosZ) / h);
    std::vector<size_t> particleGrid[gridSizeX][gridSizeY][gridSizeZ];

    void fillParticleGrid();

    struct GridCellNeighborhood
    {
        int minX, maxX;
        int minY, maxY;
        int minZ, maxZ;
    };
    GridCellNeighborhood getAdjacentCells(const glm::vec3& pos) const;
};
