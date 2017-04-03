#include "FluidSimulator.h"
#include "Kernels.h"
#include "Definitions.h"

FluidSimulator::FluidSimulator(GLFWwindow* window)
    : FluidBase(window)
{
    //TwAddVarRW(antTweakBar, "h",             TW_TYPE_FLOAT,   &h,               "min=  0.1    max=   5     step=  0.1  ");
    TwAddVarRW(antTweakBar, "k",             TW_TYPE_FLOAT,   &k,               "min=100      max=3000     step=100    ");
    TwAddVarRW(antTweakBar, "rho0",          TW_TYPE_FLOAT,   &rho0,            "min=  1      max=  50     step=  1    ");
    TwAddVarRW(antTweakBar, "m",             TW_TYPE_FLOAT,   &m,               "min=  0.1    max=   5     step=  0.1  ");
    TwAddVarRW(antTweakBar, "mu",            TW_TYPE_FLOAT,   &mu,              "min=  0.1    max=   5     step=  0.1  ");
    TwAddVarRW(antTweakBar, "sigma",         TW_TYPE_FLOAT,   &sigma,           "min=  0.001  max=   0.05  step=  0.001");
    TwAddVarRW(antTweakBar, "csNormThresh.", TW_TYPE_FLOAT,   &csNormThreshold, "min=  0.1    max=   5     step=  0.1  ");
    TwAddVarRW(antTweakBar, "gravityY",      TW_TYPE_FLOAT,   &gravity.y,       "min=-50      max=  -1     step=  1    ");
    TwAddVarRW(antTweakBar, "dt",            TW_TYPE_FLOAT,   &dt,              "min=  0.001  max=   0.05  step=  0.001");
    TwAddVarRW(antTweakBar, "Paused",        TW_TYPE_BOOLCPP, &paused,          "");
    TwAddButton(antTweakBar, "Single step", [](void* clientData)
    {
        FluidSimulator* program = static_cast<FluidSimulator*>(clientData);
        program->paused = false;
        program->update();
        program->paused = true;
    }, this, "");
    TwAddButton(antTweakBar, "Reset particles", [](void* clientData)
    {
        FluidSimulator* program = static_cast<FluidSimulator*>(clientData);
        program->resetParticles();
        program->paused = true;
    }, this, "");

    // Create cube positions
    resetParticles();

    // Fill kernel lookup tables
    fillKernelLookupTables();
}

// Updates the state of the program.
void FluidSimulator::update()
{
    if (paused)
        return;

    float add = holdShift ? 0.06f : 0.03f;

    if (holdForward)
        gridOffset.z -= add;
    else if (holdBackward)
        gridOffset.z += add;

    if (holdRight)
        gridOffset.x += add;
    else if (holdLeft)
        gridOffset.x -= add;

    if (holdUp)
        gridOffset.y += add;
    else if (holdDown)
        gridOffset.y -= add;

    worldBoundsVertices[0] = { minPos.x + gridOffset.x, minPos.y + gridOffset.y, minPos.z + gridOffset.z };
    worldBoundsVertices[1] = { maxPos.x + gridOffset.x, minPos.y + gridOffset.y, minPos.z + gridOffset.z };
    worldBoundsVertices[2] = { maxPos.x + gridOffset.x, maxPos.y + gridOffset.y, minPos.z + gridOffset.z };
    worldBoundsVertices[3] = { minPos.x + gridOffset.x, maxPos.y + gridOffset.y, minPos.z + gridOffset.z };
    worldBoundsVertices[4] = { minPos.x + gridOffset.x, minPos.y + gridOffset.y, minPos.z + gridOffset.z };
    worldBoundsVertices[5] = { minPos.x + gridOffset.x, minPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[6] = { maxPos.x + gridOffset.x, minPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[7] = { maxPos.x + gridOffset.x, maxPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[8] = { minPos.x + gridOffset.x, maxPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[9] = { minPos.x + gridOffset.x, minPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[10] = { minPos.x + gridOffset.x, maxPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[11] = { minPos.x + gridOffset.x, maxPos.y + gridOffset.y, minPos.z + gridOffset.z };
    worldBoundsVertices[12] = { maxPos.x + gridOffset.x, maxPos.y + gridOffset.y, minPos.z + gridOffset.z };
    worldBoundsVertices[13] = { maxPos.x + gridOffset.x, maxPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[14] = { maxPos.x + gridOffset.x, minPos.y + gridOffset.y, maxPos.z + gridOffset.z };
    worldBoundsVertices[15] = { maxPos.x + gridOffset.x, minPos.y + gridOffset.y, minPos.z + gridOffset.z };

#ifdef USEPARTICLEGRID
    fillParticleGrid();
#endif

    // Copy current positions to previous positions array
    std::copy(std::begin(r), std::end(r), std::begin(rPrev));

    float rho[particleCount] = {}; // Particle densities
    float p[particleCount] = {}; // Particle pressure
    // First, calculate density and pressure at each particle position
#pragma omp parallel for
    for (int i = 0; i < particleCount; ++i)
    {
        rho[i] = calcDensity(i);
        p[i] = std::max(0.0f, k * (rho[i] - rho0));
    }

    // Calculate the force acting on each particle
    glm::vec3 forces[particleCount];
#pragma omp parallel for
    for (int i = 0; i < particleCount; ++i)
    {
        glm::vec3 pressureForce = calcPressureForce(i, rho, p);

        glm::vec3 viscosityForce = calcViscosityForce(i, rho);

        glm::vec3 surfaceForce = calcSurfaceForce(i, rho);

        glm::vec3 gravityForce = gravity * rho[i];

        forces[i] = pressureForce + viscosityForce + surfaceForce + gravityForce;
    }

    // TODO: add external forces (for example, forces created by user input)

    // Move each particle
#pragma omp parallel for
    for (int i = 0; i < particleCount; ++i)
    {
        glm::vec3 a = forces[i] / rho[i]; // Acceleration
        // Semi-implicit Euler integration (TODO: better integration?)
        v[i] += a * dt;
        r[i] += v[i] * dt;

        // Rudimentary collision; the particles reside inside a hard-coded AABB
        // Upon collision with bounds, push particles out of objects, and reflect their velocity vector
        for (int dim = 0; dim != 3; ++dim) // Loop over x, y and z components
        {
            if (r[i][dim] < minPos[dim] + gridOffset[dim])
            {
                r[i][dim] = minPos[dim] + gridOffset[dim];
                v[i][dim] = -v[i][dim];
            }
            else if (r[i][dim] > maxPos[dim] + gridOffset[dim])
            {
                r[i][dim] = maxPos[dim] + gridOffset[dim];
                v[i][dim] = -v[i][dim];
            }
        }
    }
}

void FluidSimulator::resetParticles()
{
    // Create cube positions
    for (int z = 0; z != cubeSize; ++z)
        for (int y = 0; y != cubeSize; ++y)
            for (int x = 0; x != cubeSize; ++x)
                r[z * cubeSize * cubeSize + y * cubeSize + x] =
                    glm::vec3(0.3f * (x - cubeSize / 2), 0.3f * (y - cubeSize / 2), 0.3f * (z - cubeSize / 2));

    // Copy current positions to previous positions array
    std::copy(std::begin(r), std::end(r), std::begin(rPrev));

    for (auto& vel : v)
        vel = glm::vec3{};

#ifdef USEPARTICLEGRID
    fillParticleGrid();
#endif
}

void FluidSimulator::fillKernelLookupTables()
{
    // poly6 lookup table
#pragma omp parallel for
    for (int i = 0; i < lookupTableSize; ++i)
    {
        float squaredDistance = i * 1e-4f;
        poly6LookupTable[i] = poly6(squaredDistance, h);
    }
}


float FluidSimulator::calcDensity(size_t particleId) const
{
    return calcDensity(r[particleId]);
}

float FluidSimulator::calcDensity(const glm::vec3& pos) const
{
    float rho = 0.0f;
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        rho += m * poly6(pos - r[j], h);
#else
    auto neighborhood = getAdjacentCells(pos);

    for (size_t x = neighborhood.minX; x <= neighborhood.maxX; ++x)
        for (size_t y = neighborhood.minY; y <= neighborhood.maxY; ++y)
            for (size_t z = neighborhood.minZ; z <= neighborhood.maxZ; ++z)
                for (size_t j : particleGrid[x][y][z])
                {
                    {
                        float sqLen = glm::length2(pos - r[j]);
                        if (sqLen < h*h)
                        {
                            int lookupIndex = (int)(1e4f * glm::length2(pos - r[j]));
                            rho += m * poly6LookupTable[lookupIndex];

                            // Non-lookup table version:
                            //rho += m * poly6(sqLen, h);
                        }
                    }
                }
#endif
    return rho;
}

glm::vec3 FluidSimulator::calcPressureForce(size_t particleId, const float* const rho, const float* const p) const
{
    glm::vec3 pressureForce;
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        if (particleId != j && r[particleId] != r[j]) // spikyGradient will return NaN when r_i == r_j
            pressureForce += -m * ((p[particleId] + p[j]) / (2 * rho[j])) * spikyGradient(r[particleId] - r[j], h);
#else
    auto neighborhood = getAdjacentCells(r[particleId]);

    for (size_t x = neighborhood.minX; x <= neighborhood.maxX; ++x)
        for (size_t y = neighborhood.minY; y <= neighborhood.maxY; ++y)
            for (size_t z = neighborhood.minZ; z <= neighborhood.maxZ; ++z)
                for (size_t j : particleGrid[x][y][z])
                    if (particleId != j && r[particleId] != r[j]) // spikyGradient will return NaN when r_i == r_j
                        pressureForce += -m * ((p[particleId] + p[j]) / (2 * rho[j])) * spikyGradient(r[particleId] - r[j], h);
#endif
    return pressureForce;
}

glm::vec3 FluidSimulator::calcViscosityForce(size_t particleId, const float* const rho) const
{
    glm::vec3 viscosityForce;
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        if (particleId != j)
            viscosityForce += mu * m * ((v[j] - v[particleId]) / rho[j]) * viscosityLaplacian(r[particleId] - r[j], h);
#else
    auto neighborhood = getAdjacentCells(r[particleId]);

    for (size_t x = neighborhood.minX; x <= neighborhood.maxX; ++x)
        for (size_t y = neighborhood.minY; y <= neighborhood.maxY; ++y)
            for (size_t z = neighborhood.minZ; z <= neighborhood.maxZ; ++z)
                for (size_t j : particleGrid[x][y][z])
                    if (particleId != j)
                        viscosityForce += mu * m * ((v[j] - v[particleId]) / rho[j]) * viscosityLaplacian(r[particleId] - r[j], h);
#endif

    return viscosityForce;
}

glm::vec3 FluidSimulator::calcSurfaceForce(size_t particleId, const float* const rho) const
{
    glm::vec3 n; // Gradient field of the color field
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        n += m * (1 / rho[j]) * poly6Gradient(r[particleId] - r[j], h);
#else
    auto neighborhood = getAdjacentCells(r[particleId]);

    for (size_t x = neighborhood.minX; x <= neighborhood.maxX; ++x)
        for (size_t y = neighborhood.minY; y <= neighborhood.maxY; ++y)
            for (size_t z = neighborhood.minZ; z <= neighborhood.maxZ; ++z)
                for (size_t j : particleGrid[x][y][z])
                    n += m * (1 / rho[j]) * poly6Gradient(r[particleId] - r[j], h);
#endif

    if (glm::length(n) < csNormThreshold)
        return glm::vec3{};

    float csLaplacian = 0.0f; // Laplacian of the color field
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        csLaplacian += m * (1 / rho[j]) * poly6Laplacian(r[particleId] - r[j], h);
#else
    for (size_t x = neighborhood.minX; x <= neighborhood.maxX; ++x)
        for (size_t y = neighborhood.minY; y <= neighborhood.maxY; ++y)
            for (size_t z = neighborhood.minZ; z <= neighborhood.maxZ; ++z)
                 for (size_t j : particleGrid[x][y][z])
                     csLaplacian += m * (1 / rho[j]) * poly6Laplacian(r[particleId] - r[j], h);
#endif
     return -sigma * csLaplacian * glm::normalize(n);
}

void FluidSimulator::fillParticleGrid()
{
    // Possible TODO: don't use vectors for fast clearing using memset
    // However, a fixed size could possibly cost a lot of memory for larger grids
    for (size_t x = 0; x != gridSizeX; ++x)
        for (size_t y = 0; y != gridSizeY; ++y)
            for (size_t z = 0; z != gridSizeZ; ++z)
                particleGrid[x][y][z].clear();

    for (size_t i = 0; i != particleCount; ++i)
    {
        const glm::vec3& pos = r[i];

        // Subtract EPSILON to make sure values exactly at grid edge don't lead to incorrect array slot
        particleGrid
            [clamp(0, gridSizeX - 1, (int)((pos.x - (minPos.x + gridOffset.x) - EPSILON) / h))]
        [clamp(0, gridSizeY - 1, (int)((pos.y - (minPos.y + gridOffset.y) - EPSILON) / h))]
        [clamp(0, gridSizeZ - 1, (int)((pos.z - (minPos.z + gridOffset.z) - EPSILON) / h))].push_back(i);
    }
}

FluidSimulator::GridCellNeighborhood FluidSimulator::getAdjacentCells(const glm::vec3& pos) const
{
    int gridX = (int)((pos.x - (minPos.x + gridOffset.x) - EPSILON) / h);
    int gridY = (int)((pos.y - (minPos.y + gridOffset.y) - EPSILON) / h);
    int gridZ = (int)((pos.z - (minPos.z + gridOffset.z) - EPSILON) / h);
    return {
        std::max(gridX - 1, 0), std::min(gridX + 1, gridSizeX - 1),
        std::max(gridY - 1, 0), std::min(gridY + 1, gridSizeY - 1),
        std::max(gridZ - 1, 0), std::min(gridZ + 1, gridSizeZ - 1)
    };
}
