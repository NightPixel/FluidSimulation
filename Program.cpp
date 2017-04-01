#include "Program.h"
#include "Kernels.h"
#include "OpenGLUtils.h"
#include "Randomizer.h"
#include "Definitions.h"
#include "Utils.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Callback function called by AntTweakBar when the single-step button is clicked.
void TW_CALL stepButtonCallback(void* clientData)
{
    Program* program = static_cast<Program*>(clientData);

    program->paused = false;
    program->update();
    program->paused = true;
}

// Callback function called by AntTweakBar when the particle reset button is clicked.
void TW_CALL particleResetButtonCallback(void* clientData)
{
    Program* program = static_cast<Program*>(clientData);

    program->resetParticles();
    program->paused = true;
}

Program::Program(GLFWwindow* window)
    : ProgramBase(window)
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
    TwAddButton(antTweakBar, "Single step", stepButtonCallback, this, "");
    TwAddButton(antTweakBar, "Reset particles", particleResetButtonCallback, this, "");


    // Set up fragment shader uniforms
    glUseProgram(simpleShaderProgram);
    glUniform3fv(glGetUniformLocation(simpleShaderProgram, "minWorldPos"), 1, glm::value_ptr(minPos));
    glUniform3fv(glGetUniformLocation(simpleShaderProgram, "maxWorldPos"), 1, glm::value_ptr(maxPos));
    
    // Create cube positions
    resetParticles();
}

// Updates the state of the program.
void Program::update()
{
    if (paused)
        return;
    /* DEBUG: Randomize positions to test CPU-to-GPU vertex streaming */
    //for (auto& vert : r)
    //    vert += Randomizer::random(-0.001f, 0.001f);
    /* END DEBUG */

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

// Draws a new frame.
void Program::draw()
{
    // Run marching cubes using PolyVox, and retrieve the vertex and index buffers
    fillVoxelVolume();
    surfaceExtractor.execute();
    const auto& lowerCorner = voxelVolume.getEnclosingRegion().getLowerCorner();
    surfaceMesh.translateVertices({(float)lowerCorner.getX(), (float)lowerCorner.getY(), (float)lowerCorner.getZ()});
    surfaceMesh.scaleVertices(1.0f / voxelVolumeResolutionScale);
    const std::vector<uint32_t>& indices = surfaceMesh.getIndices();
    std::vector<PolyVox::PositionMaterialNormal>& vertices = surfaceMesh.getRawVertexData();
    for (auto& vert : vertices) // Clamp vertex locations to world boundaries
    {
        vert.position.setElements(
            clamp(minPos.x + gridOffset.x, maxPos.x + gridOffset.x, vert.position.getX() + gridOffset.x),
            clamp(minPos.y + gridOffset.y, maxPos.y + gridOffset.y, vert.position.getY() + gridOffset.y),
            clamp(minPos.z + gridOffset.z, maxPos.z + gridOffset.z, vert.position.getZ() + gridOffset.z)
        );
    }

    const glm::mat4 view = camera.getViewMatrix();
    glUseProgram(waterShaderProgram);
    glUniformMatrix4fv(waterViewUniform, 1, GL_FALSE, glm::value_ptr(view));
    glUniform3fv(waterCamUniform, 1, glm::value_ptr(camera.getPosition()));

    // Draw surface mesh
    glBindVertexArray(meshVAO);
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(PolyVox::PositionMaterialNormal), vertices.data(), GL_STREAM_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(uint32_t), indices.data(), GL_STREAM_DRAW);
    glDrawElements(GL_TRIANGLES, (GLsizei)indices.size(), GL_UNSIGNED_INT, nullptr);

    glUseProgram(simpleShaderProgram);
    glUniformMatrix4fv(simpleViewUniform, 1, GL_FALSE, glm::value_ptr(view));

    // Draw points
    glBindVertexArray(pointsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(r), r, GL_STREAM_DRAW);
    glDrawArrays(GL_POINTS, 0, (GLsizei)std::size(r));

    // Draw world bounds
    glBufferData(GL_ARRAY_BUFFER, sizeof(worldBoundsVertices), worldBoundsVertices, GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)std::size(worldBoundsVertices));
}

void Program::resetParticles()
{
    // Create cube positions
    for (int z = 0; z != cubeSize; ++z)
        for (int y = 0; y != cubeSize; ++y)
            for (int x = 0; x != cubeSize; ++x)
                r[z * cubeSize * cubeSize + y * cubeSize + x] =
                    glm::vec3(0.3f * (x - cubeSize / 2), 0.3f * (y - cubeSize / 2), 0.3f * (z - cubeSize / 2));

    for (auto& vel : v)
        vel = glm::vec3{};

#ifdef USEPARTICLEGRID
    fillParticleGrid();
#endif
}

void Program::fillParticleGrid()
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

float Program::calcDensity(size_t particleId) const
{
    return calcDensity(r[particleId]);
}

float Program::calcDensity(const glm::vec3& pos) const
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
                    rho += m * poly6(pos - r[j], h);
#endif
    return rho;
}

glm::vec3 Program::calcPressureForce(size_t particleId, float* rho, float* p)
{
    glm::vec3 pressureForce;
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        if (particleId != j)
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

glm::vec3 Program::calcViscosityForce(size_t particleId, float* rho)
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

glm::vec3 Program::calcSurfaceForce(size_t particleId, float* rho)
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

Program::GridCellNeighborhood Program::getAdjacentCells(const glm::vec3& pos) const
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

PolyVox::Vector3DInt32 Program::worldPosToVoxelIndex(const glm::vec3& worldPos) const
{
    return {
        (int)((worldPos.x - gridOffset.x) * voxelVolumeResolutionScale),
        (int)((worldPos.y - gridOffset.y) * voxelVolumeResolutionScale),
        (int)((worldPos.z - gridOffset.z) * voxelVolumeResolutionScale)
    };
}

void Program::fillVoxelVolume()
{
    const PolyVox::Region volumeRegion = voxelVolume.getEnclosingRegion();
    const PolyVox::Vector3DInt32& lowerCorner = volumeRegion.getLowerCorner();
    const PolyVox::Vector3DInt32& upperCorner = volumeRegion.getUpperCorner();

    float invVoxelVolumeResolutionScale = 1.0f / voxelVolumeResolutionScale;
#pragma omp parallel for
    for (int z = lowerCorner.getZ(); z <= upperCorner.getZ(); z++)
        for (int y = lowerCorner.getY(); y <= upperCorner.getY(); y++)
            for (int x = lowerCorner.getX(); x <= upperCorner.getX(); x++)
                voxelVolume.setVoxelAt(x, y, z, calcDensity(
                    { x * invVoxelVolumeResolutionScale + gridOffset.x, y * invVoxelVolumeResolutionScale + gridOffset.y, z * invVoxelVolumeResolutionScale + gridOffset.z }));

    /* DEBUG: Only fill voxels that contain an actual particle */
    /*for (int z = lowerCorner.getZ(); z <= upperCorner.getZ(); z++)
        for (int y = lowerCorner.getY(); y <= upperCorner.getY(); y++)
            for (int x = lowerCorner.getX(); x <= upperCorner.getX(); x++)
                voxelVolume.setVoxelAt(x, y, z, 0);
    for (const auto& pos : r)
        voxelVolume.setVoxelAt(worldPosToVoxelIndex(pos), 100.0f);
    /* END DEBUG */
}