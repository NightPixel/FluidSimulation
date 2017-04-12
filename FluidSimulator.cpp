#include "FluidSimulator.h"
#include "Kernels.h"
#include "Definitions.h"
#include <tuple>

FluidSimulator::FluidSimulator(GLFWwindow* window)
    : FluidBase(window)
{
    //TwAddVarRW(antTweakBar, "h",             TW_TYPE_FLOAT,   &h,               "min=  0.1    max=   5     step=  0.1  ");
    TwAddVarRW(antTweakBar, "k",               TW_TYPE_FLOAT,   &k,               "min=100      max=3000     step=100    ");
    TwAddVarRW(antTweakBar, "rho0",            TW_TYPE_FLOAT,   &rho0,            "min=  1      max=  50     step=  1    ");
    TwAddVarRW(antTweakBar, "m",               TW_TYPE_FLOAT,   &m,               "min=  0.1    max=   5     step=  0.1  ");
    TwAddVarRW(antTweakBar, "mu",              TW_TYPE_FLOAT,   &mu,              "min=  0.1    max=   15     step=  0.1  ");
    TwAddVarRW(antTweakBar, "sigma",           TW_TYPE_FLOAT,   &sigma,           "min=  0.001  max=   0.05  step=  0.001");
    TwAddVarRW(antTweakBar, "csNormThresh.",   TW_TYPE_FLOAT,   &csNormThreshold, "min=  0.1    max=   5     step=  0.1  ");
    TwAddVarRW(antTweakBar, "gravityY",        TW_TYPE_FLOAT,   &gravity.y,       "min=-50      max=  -1     step=  1    ");
    TwAddVarRW(antTweakBar, "restitution cf.", TW_TYPE_FLOAT,   &cr,              "min=  0      max=   1     step=  0.1  ");
    TwAddVarRW(antTweakBar, "dt",              TW_TYPE_FLOAT,   &dt,              "min=  0.001  max=   0.05  step=  0.001");
    TwAddVarRW(antTweakBar, "Paused",          TW_TYPE_BOOLCPP, &paused,          "");
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

    fillTriangleGrid();
}

// Updates the state of the program.
void FluidSimulator::update()
{
    if (paused)
        return;

    float add = holdShift ? 0.06f : 0.03f;
    glm::vec3 sceneMovement = glm::vec3();

    if (holdForward)
        sceneMovement.z -= add;
    else if (holdBackward)
        sceneMovement.z += add;

    if (holdRight)
        sceneMovement.x += add;
    else if (holdLeft)
        sceneMovement.x -= add;

    if (holdUp)
        sceneMovement.y += add;
    else if (holdDown)
        sceneMovement.y -= add;

    if (sceneMovement != glm::vec3())
    {
        sceneOffset += sceneMovement;
        // Handle collisions caused by scenemovement for each water particle

        for (int i = 0; i < particleCount; ++i)
        {
            handleCollisions(i, r[i] + sceneMovement);
        }
    }

#ifdef USEPARTICLEGRID
    fillParticleGrid();
#endif

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
        const glm::vec3 a = forces[i] / rho[i]; // Acceleration
        // Semi-implicit Euler integration (TODO: better integration?)
        v[i] += a * dt;
        r[i] += v[i] * dt;

        handleCollisions(i, r[i] - v[i] * dt);

        // Rudimentary collision; the particles reside inside a hard-coded AABB
        // Upon collision with bounds, push particles out of objects, and reflect their velocity vector
        for (int dim = 0; dim != 3; ++dim) // Loop over x, y and z components
        {
            if (r[i][dim] < minPos[dim] + sceneOffset[dim])
            {
                r[i][dim] = minPos[dim] + sceneOffset[dim];
                v[i][dim] = -v[i][dim];
            }
            else if (r[i][dim] > maxPos[dim] + sceneOffset[dim])
            {
                r[i][dim] = maxPos[dim] + sceneOffset[dim];
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
                    glm::vec3(0.2f * (x - cubeSize / 2), 0.2f * (y - cubeSize / 2), 0.2f * (z - cubeSize / 2)) + sceneOffset;

    for (auto& vel : v)
        vel = glm::vec3{};

#ifdef USEPARTICLEGRID
    fillParticleGrid();
#endif
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

    __m128 posQF = _mm_loadu_ps(&pos.x);
    float hSq = h * h;
    for (size_t x = neighborhood.minX; x <= neighborhood.maxX; ++x)
        for (size_t y = neighborhood.minY; y <= neighborhood.maxY; ++y)
            for (size_t z = neighborhood.minZ; z <= neighborhood.maxZ; ++z)
                for (size_t j : particleGrid[x][y][z])
                {
                    __m128 rQF = _mm_loadu_ps(&r[j].x);
                    __m128 rMinusPosQF = _mm_sub_ps(posQF, rQF);
                    __m128 multPosRQF = _mm_mul_ps(rMinusPosQF, rMinusPosQF);
                    float* multPosRArray = (float*)&multPosRQF;
                    float sqLen = multPosRArray[0] + multPosRArray[1] + multPosRArray[2];
                    if (sqLen < hSq)
                    {
                        int lookupIndex = (int)(1e4f * sqLen);
                        rho += poly6LookupTable[lookupIndex];

                        // Non-lookup table version:
                        //rho += m * poly6(sqLen, h);
                    }
                }
    rho *= m;
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

void FluidSimulator::handleCollisions(size_t particleId, glm::vec3 oldPos)
{
const glm::ivec3 oldGridIndex = glm::clamp({0, 0, 0}, worldPosToGridIndex(oldPos), {gridSizeX - 1, gridSizeY - 1, gridSizeZ - 1});
    const glm::ivec3 newGridIndex = glm::clamp({0, 0, 0}, worldPosToGridIndex(r[particleId]), {gridSizeX - 1, gridSizeY - 1, gridSizeZ - 1});

    std::pair<bool, float> closestIntersection(false, std::numeric_limits<float>::max());
    Triangle* closestIntersectionTriangle = nullptr;

    for (int x = std::min(oldGridIndex.x, newGridIndex.x); x <= std::max(oldGridIndex.x, newGridIndex.x); ++x)
        for (int y = std::min(oldGridIndex.y, newGridIndex.y); y <= std::max(oldGridIndex.y, newGridIndex.y); ++y)
            for (int z = std::min(oldGridIndex.z, newGridIndex.z); z <= std::max(oldGridIndex.z, newGridIndex.z); ++z)
            {
                std::vector<Triangle*>& trianglesInCell = triangleGrid[x][y][z];

                for (Triangle* triangle : trianglesInCell)
                {
                    auto intersection = triangleLineSegmentIntersection(*triangle, oldPos, r[particleId], sceneOffset);
                    if (!intersection.first)
                        continue; // No collision

                    if (intersection.second < closestIntersection.second)
                    {
                        closestIntersection = intersection;
                        closestIntersectionTriangle = triangle;
                    }
                }
            }

    if (closestIntersection.first)
    {
        // One (or maybe even more) collisions were found: resolve using the closest one

        const glm::vec3 newPos = r[particleId];
        const glm::vec3 delta = (r[particleId] - oldPos);
        float deltaNorm = glm::length(delta);
        // In case of collision with a triangle,
        // we project the particle back on the triangle...
        // We project it to be just a little bit outside of the triangle to prevent the particle from ending up "in the face" and not detecting a collision in the next frame
        r[particleId] = r[particleId] - closestIntersectionTriangle->normal * (glm::dot(r[particleId] - (closestIntersectionTriangle->positions[0] + sceneOffset), closestIntersectionTriangle->normal) - 1e-3f);

        // ... and reflect the velocity vector along the surface normal,
        // scaled by the coefficient of restitution and the penetration depth (between 0 and 1)
        // v = v - (1 + cR * penetration depth) * (v . n) * n
        v[particleId] = v[particleId] - (1.0f + cr * closestIntersection.second) * glm::dot(v[particleId], closestIntersectionTriangle->normal) * closestIntersectionTriangle->normal;
    }
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
        const glm::ivec3 gridPos = worldPosToGridIndex(r[i]);
        particleGrid
            [clamp(0, gridSizeX - 1, gridPos.x)]
            [clamp(0, gridSizeY - 1, gridPos.y)]
            [clamp(0, gridSizeZ - 1, gridPos.z)].push_back(i);
    }
}

void FluidSimulator::fillTriangleGrid()
{
    // Pre-compute the bounding boxes for each grid cell.
    // Each grid cell's bounding box is described by a (grid center, half edge lengths) vector pair.
    std::pair<glm::vec3, glm::vec3> gridCellBoundingBoxes[gridSizeX][gridSizeY][gridSizeZ];

    for (size_t x = 0; x != gridSizeX; ++x)
        for (size_t y = 0; y != gridSizeY; ++y)
            for (size_t z = 0; z != gridSizeZ; ++z)
            {
                triangleGrid[x][y][z].clear();

                const glm::vec3 lowerCorner{
                    x * h + (minPos.x + sceneOffset.x),
                    y * h + (minPos.y + sceneOffset.y),
                    z * h + (minPos.z + sceneOffset.z)
                };
                const glm::vec3 upperCorner{
                    (x + 1) * h + (minPos.x + sceneOffset.x),
                    (y + 1) * h + (minPos.y + sceneOffset.y),
                    (z + 1) * h + (minPos.z + sceneOffset.z)
                };

                const glm::vec3 boxCenter = (lowerCorner + upperCorner) / 2.0f;
                const glm::vec3 halfLengths = (upperCorner - lowerCorner) / 2.0f;

                gridCellBoundingBoxes[x][y][z] = std::make_pair(boxCenter, halfLengths);
            }

    // For each triangle, we find the grid cells in which its bounding box resides.
    // For each of these grid cells, we do a triangle-box intersection test;
    // if the grid cell's box intersects the triangle, we store a pointer to the triangle in that grid cell.
    for (auto& model : models)
    {
        for (auto& triangle : model.triangles)
        {
            // Each triangle's bounding box is described by a (lower corner, upper corner) vector pair.
            // Note that this is a different description than the grid cell's bounding boxes!
            const auto bounds = triangle.getBoundingBox();

            // Convert the triangle bounding box to grid cells
            const glm::ivec3 lowerGridPos = glm::clamp({0, 0, 0}, worldPosToGridIndex(bounds.first), {gridSizeX - 1, gridSizeY - 1, gridSizeZ - 1});
            const glm::ivec3 upperGridPos = glm::clamp({0, 0, 0}, worldPosToGridIndex(bounds.second), {gridSizeX - 1, gridSizeY - 1, gridSizeZ - 1});

            // For each grid cell that overlaps with the triangle's bounding box...
            for (size_t x = lowerGridPos.x; x <= upperGridPos.x; ++x)
                for (size_t y = lowerGridPos.y; y <= upperGridPos.y; ++y)
                    for (size_t z = lowerGridPos.z; z <= upperGridPos.z; ++z)
                    {
                        const glm::vec3& gridCellCenter = gridCellBoundingBoxes[x][y][z].first;
                        const glm::vec3& gridCellHalfLengths = gridCellBoundingBoxes[x][y][z].second;

                        // ...if the triangle actually intersects with the grid cell...
                        if (triangleBoxIntersection(triangle, gridCellCenter, gridCellHalfLengths))
                        {
                            // ...store the triangle in that cell.
                            triangleGrid[x][y][z].push_back(&triangle);
                            printf("(%i, %i, %i) (%i, %i, %i)\n", lowerGridPos.x, lowerGridPos.y, lowerGridPos.z, upperGridPos.x, upperGridPos.y, upperGridPos.z);
                        }
                    }
        }
    }
}

glm::ivec3 FluidSimulator::worldPosToGridIndex(const glm::vec3& pos) const
{
    // Subtract EPSILON to make sure values exactly at grid edge don't lead to incorrect array slot
    return {
        (int)((pos.x - (minPos.x + sceneOffset.x) - EPSILON) / h),
        (int)((pos.y - (minPos.y + sceneOffset.y) - EPSILON) / h),
        (int)((pos.z - (minPos.z + sceneOffset.z) - EPSILON) / h)
    };
}

FluidSimulator::GridCellNeighborhood FluidSimulator::getAdjacentCells(const glm::vec3& pos) const
{
    const glm::ivec3 gridPos = worldPosToGridIndex(pos);
    return {
        std::max(gridPos.x - 1, 0), std::min(gridPos.x + 1, gridSizeX - 1),
        std::max(gridPos.y - 1, 0), std::min(gridPos.y + 1, gridSizeY - 1),
        std::max(gridPos.z - 1, 0), std::min(gridPos.z + 1, gridSizeZ - 1)
    };
}
