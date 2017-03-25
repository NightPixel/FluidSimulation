#include "Program.h"
#include "Kernels.h"
#include "OpenGLUtils.h"
#include "Randomizer.h"
#include "Definitions.h"
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
    : window(window)
    , voxelVolume(PolyVox::Region(worldPosToVoxelIndex(minPos), worldPosToVoxelIndex(maxPos)))
    , surfaceExtractor(&voxelVolume, voxelVolume.getEnclosingRegion(), &surfaceMesh)
{
    glfwGetWindowSize(window, &windowSizeX, &windowSizeY);
    antTweakBar = TwNewBar("Simulation settings");
    TwDefine("GLOBAL fontsize=3");
    TwAddVarRW(antTweakBar, "h",             TW_TYPE_FLOAT,   &h,               "min=  0.1    max=   5     step=  0.1  ");
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

    // Initialize OpenGL
    glEnable(GL_DEPTH_TEST);
    glPointSize(5.0f);

    // Create Vertex Array Object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Create cube positions
    resetParticles();

    // Create a Vertex Buffer Object
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
#ifdef DRAWPOINTS
    // Copy the vertex data
    glBufferData(GL_ARRAY_BUFFER, sizeof(r), r, GL_STREAM_DRAW);
#endif

    // Create an Element Buffer Object
    glGenBuffers(1, &ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

    // Create and compile the vertex shader
    vertexShader = createShaderFromSource("Simple.vert", GL_VERTEX_SHADER);
    auto shaderInfo = checkShaderCompilation(vertexShader);
    if (!shaderInfo.first)
        printf("Vertex shader failed to compile!\n%s\n", shaderInfo.second.c_str());
    else if (!shaderInfo.second.empty())
        printf("Vertex shader compiled with warnings.\n%s\n", shaderInfo.second.c_str());

    // Create and compile the fragment shader
    fragmentShader = createShaderFromSource("Simple.frag", GL_FRAGMENT_SHADER);
    shaderInfo = checkShaderCompilation(fragmentShader);
    if (!shaderInfo.first)
        printf("Vertex shader failed to compile!\n%s\n", shaderInfo.second.c_str());
    else if (!shaderInfo.second.empty())
        printf("Vertex shader compiled with warnings.\n%s\n", shaderInfo.second.c_str());

    // Link the vertex and fragment shader into a shader program
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glBindFragDataLocation(shaderProgram, 0, "outColor");
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    // Specify the layout of the vertex data
#ifdef DRAWPOINTS
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<void*>(0 * sizeof(float)));
#else
    // PolyVox::PositionMaterialNormal layout:
    //    (x, y, z)-position (3 floats)
    //    (x, y, z)-normal   (3 floats)
    //    material           (1 float)
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(0 * sizeof(float)));
    GLint normAttrib = glGetAttribLocation(shaderProgram, "normal");
    glEnableVertexAttribArray(normAttrib);
    glVertexAttribPointer(normAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(3 * sizeof(float)));
    GLint matAttrib = glGetAttribLocation(shaderProgram, "material");
    glEnableVertexAttribArray(matAttrib);
    glVertexAttribPointer(matAttrib, 1, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(6 * sizeof(float)));
#endif

    // Set up model, view, projection matrices
    glm::mat4 model{}; // Identity matrix
    GLint uniModel = glGetUniformLocation(shaderProgram, "model");
    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

    uniView = glGetUniformLocation(shaderProgram, "view");

    glm::mat4 proj = glm::perspective(glm::radians(45.0f), (float)windowSizeX / windowSizeY, 1.0f, 25.0f);
    GLint uniProj = glGetUniformLocation(shaderProgram, "proj");
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));
}

Program::~Program()
{
    glDeleteProgram(shaderProgram);
    glDeleteShader(fragmentShader);
    glDeleteShader(vertexShader);

    glDeleteBuffers(1, &ebo);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
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

#ifdef USEPARTICLEGRID
    fillParticleGrid();
#endif

    float rho[particleCount] = {}; // Particle densities
    float p[particleCount] = {}; // Particle pressure
    // First, calculate density and pressure at each particle position
    for (size_t i = 0; i != particleCount; ++i)
    {
        rho[i] = calcDensity(i);
        p[i] = std::max(0.0f, k * (rho[i] - rho0));
    }

    // Calculate the force acting on each particle
    glm::vec3 forces[particleCount];
    for (size_t i = 0; i != particleCount; ++i)
    {
        glm::vec3 pressureForce = calcPressureForce(i, rho, p);

        glm::vec3 viscosityForce = calcViscosityForce(i, rho);

        glm::vec3 surfaceForce = calcSurfaceForce(i, rho);

        glm::vec3 gravityForce = gravity * rho[i];

        forces[i] = pressureForce + viscosityForce + surfaceForce + gravityForce;
    }

    // TODO: add external forces (for example, forces created by user input)

    // Move each particle
    for (size_t i = 0; i != particleCount; ++i)
    {
        glm::vec3 a = forces[i] / rho[i]; // Acceleration
        // Semi-implicit Euler integration (TODO: better integration?)
        v[i] += a * dt;
        r[i] += v[i] * dt;

        // Rudimentary collision; the particles reside inside a hard-coded AABB
        // Upon collision with bounds, push particles out of objects, and reflect their velocity vector
        for (int dim = 0; dim != 3; ++dim) // Loop over x, y and z components
        {
            if (r[i][dim] < minPos[dim])
            {
                r[i][dim] = minPos[dim];
                v[i][dim] = -v[i][dim];
            }
            else if (r[i][dim] > maxPos[dim])
            {
                r[i][dim] = maxPos[dim];
                v[i][dim] = -v[i][dim];
            }
        }
    }

}

// Draws a new frame.
void Program::draw()
{
    // Re-bind vertex buffer; AntTweakBar changed it.
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Update view matrix
    glm::mat4 view = camera.getViewMatrix();
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

#ifdef DRAWPOINTS
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(r), r);
    glDrawArrays(GL_POINTS, 0, cubeSize * cubeSize * cubeSize);
#else
    fillVoxelVolume();
    surfaceExtractor.execute();
    surfaceMesh.scaleVertices(1.0f / (2.0f * voxelVolumeResolutionScale));
    const std::vector<uint32_t>& indices = surfaceMesh.getIndices();
    const std::vector<PolyVox::PositionMaterialNormal>& vertices = surfaceMesh.getVertices();

    // Upload new vertex and index data
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(PolyVox::PositionMaterialNormal), vertices.data(), GL_STREAM_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(uint32_t), indices.data(), GL_STREAM_DRAW);

    glDrawElements(GL_TRIANGLES, (GLsizei)indices.size(), GL_UNSIGNED_INT, nullptr);
#endif
}

void Program::resetParticles()
{
    // Create cube positions
    for (int z = 0; z != cubeSize; ++z)
        for (int y = 0; y != cubeSize; ++y)
            for (int x = 0; x != cubeSize; ++x)
            {
                r[z * cubeSize * cubeSize + y * cubeSize + x] =
                    glm::vec3(0.3f * (x - cubeSize / 2), 0.3f * (y - cubeSize / 2), 0.3f * (z - cubeSize / 2));
                printf("[%i] = (%f, %f, %f)\n",
                    z * cubeSize * cubeSize + y * cubeSize + x,
                    0.3f * (x - cubeSize / 2), 0.3f * (y - cubeSize / 2), 0.3f * (z - cubeSize / 2));
            }

    for (auto& vel : v)
        vel = glm::vec3{};
}

// Called when the mouse cursor is moved.
void Program::onMouseMoved(float dxPos, float dyPos)
{
    const bool leftMouseButtonPressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    const bool rightMouseButtonPressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;

    if (leftMouseButtonPressed == rightMouseButtonPressed)
        return; // Neither is pressed, or both are pressed

    if (leftMouseButtonPressed)
    {
        float dTheta = -dxPos / (0.5f * windowSizeX);
        float dPhi   = -dyPos / (0.5f * windowSizeY);
        camera.rotate(dTheta, dPhi);
    }
    if (rightMouseButtonPressed)
    {
        float dx =  2.0f * dxPos / (0.5f * windowSizeX);
        float dy = -2.0f * dyPos / (0.5f * windowSizeY);
        camera.pan(dx, dy);
    }

}

// Called when the mouse wheel is scrolled.
void Program::onMouseScrolled(float yOffset)
{
    camera.zoom(yOffset / 30.0f);
}

void Program::fillParticleGrid()
{
    // Possible TODO: don't use vectors for fast clearing using memset
    // However, a fixed size could possibly cost a lot of memory for larger grids
    for(size_t x = 0; x != 4; ++x)
        for (size_t y = 0; y != 4; ++y)
            for (size_t z = 0; z != 4; ++z)
            {
                particleGrid[x][y][z].clear();
            }
    

    for (size_t i = 0; i != particleCount; ++i)
    {
        glm::vec3 pos = r[i];
        
        // Make sure values exactly at grid edge don't lead to incorrect array slot
        particleGrid[(int)(pos.x + 2.0f - EPSILON)][(int)(pos.y + 2.0f - EPSILON)][(int)(pos.z + 2.0f - EPSILON)].push_back(i);
    }
}

glm::vec3 Program::calcPressureForce(size_t particleId, float* rho, float* p)
{
    glm::vec3 pressureForce;
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        pressureForce += -m * ((p[particleId] + p[j]) / (2 * rho[j])) * spikyGradient(r[particleId] - r[j], h);
#else
    glm::vec3 pos = r[particleId];
    int gridX = (int)(pos.x + 2.0f - EPSILON); int minX, maxX;
    int gridY = (int)(pos.y + 2.0f - EPSILON); int minY, maxY;
    int gridZ = (int)(pos.z + 2.0f - EPSILON); int minZ, maxZ;
    getAdjacentCells(gridX, gridY, gridZ, minX, maxX, minY, maxY, minZ, maxZ);

    for (size_t x = minX; x <= maxX; ++x) for (size_t y = minY; y <= maxY; ++y) for (size_t z = minZ; z <= maxZ; ++z)
        for (size_t j : particleGrid[x][y][z])
        {
            pressureForce += -m * ((p[particleId] + p[j]) / (2 * rho[j])) * spikyGradient(r[particleId] - r[j], h);
            
            float dot = glm::dot(spikyGradient(r[j] - r[particleId], h), r[j] - r[particleId]);
            if (dot > 0.0f)
                printf("Pos\n");

            float dot2 = glm::dot(-m * ((p[particleId] + p[j]) / (2 * rho[j])) * spikyGradient(r[particleId] - r[j], h), r[j] - r[particleId]);
            if (dot2 > 0.0f)
                printf("Pos2\n");
        }
#endif
    return pressureForce;
}

glm::vec3 Program::calcViscosityForce(size_t particleId, float* rho)
{
    glm::vec3 viscosityForce;
#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        viscosityForce += mu * m * ((v[j] - v[particleId]) / rho[j]) * viscosityLaplacian(r[particleId] - r[j], h);
#else
    glm::vec3 pos = r[particleId];
    int gridX = (int)(pos.x + 2.0f - EPSILON); int minX, maxX;
    int gridY = (int)(pos.y + 2.0f - EPSILON); int minY, maxY;
    int gridZ = (int)(pos.z + 2.0f - EPSILON); int minZ, maxZ;
    getAdjacentCells(gridX, gridY, gridZ, minX, maxX, minY, maxY, minZ, maxZ);

    for (size_t x = minX; x <= maxX; ++x) for (size_t y = minY; y <= maxY; ++y) for (size_t z = minZ; z <= maxZ; ++z)
        for (size_t j : particleGrid[x][y][z])
            viscosityForce += mu * m * ((v[j] - v[particleId]) / rho[j]) * viscosityLaplacian(r[particleId] - r[j], h);
#endif

    return viscosityForce;
}

glm::vec3 Program::calcSurfaceForce(size_t particleId, float* rho)
{
    glm::vec3 surfaceForce;
    glm::vec3 n; // Gradient field of the color field

#ifndef USEPARTICLEGRID
    for (size_t j = 0; j != particleCount; ++j)
        n += m * (1 / rho[j]) * poly6Gradient(r[particleId] - r[j], h);
#else
    glm::vec3 pos = r[particleId];
    int gridX = (int)(pos.x + 2.0f - EPSILON); int minX, maxX;
    int gridY = (int)(pos.y + 2.0f - EPSILON); int minY, maxY;
    int gridZ = (int)(pos.z + 2.0f - EPSILON); int minZ, maxZ;
    getAdjacentCells(gridX, gridY, gridZ, minX, maxX, minY, maxY, minZ, maxZ);

    for (size_t x = minX; x <= maxX; ++x) for (size_t y = minY; y <= maxY; ++y) for (size_t z = minZ; z <= maxZ; ++z)
        for (size_t j : particleGrid[x][y][z])
        {
            n += m * (1 / rho[j]) * poly6Gradient(r[particleId] - r[j], h);
        }
#endif

    if (glm::length(n) > csNormThreshold)
    {
        float csLaplacian = 0.0f; // Laplacian of the color field

#ifndef USEPARTICLEGRID
        for (size_t j = 0; j != particleCount; ++j)
            csLaplacian += m * (1 / rho[j]) * poly6Laplacian(r[particleId] - r[j], h);
#else
        for (size_t x = minX; x <= maxX; ++x) for (size_t y = minY; y <= maxY; ++y) for (size_t z = minZ; z <= maxZ; ++z)
            for (size_t j : particleGrid[x][y][z])
            {
                csLaplacian += m * (1 / rho[j]) * poly6Laplacian(r[particleId] - r[j], h);
            }
#endif
        surfaceForce = -sigma * csLaplacian * glm::normalize(n);
    }

    return surfaceForce;
}

float Program::calcDensity(size_t particleId) const
{
    float rho = 0.0f;
    for (size_t j = 0; j != particleCount; ++j)
        rho += m * poly6(r[particleId] - r[j], h);
    return rho;
}

float Program::calcDensity(const glm::vec3& position) const
{
    float rho = 0.0f;
    for (size_t j = 0; j != particleCount; ++j)
        rho += m * poly6(position - r[j], h);
    return rho;
}

void Program::getAdjacentCells(int gridX, int gridY, int gridZ, int & minXOut, int & maxXOut, int & minYOut, int & maxYOut, int & minZOut, int & maxZOut)
{
    minXOut = std::max(gridX - 1, 0); maxXOut = std::min(gridX + 1, 3);
    minYOut = std::max(gridY - 1, 0); maxYOut = std::min(gridY + 1, 3);
    minZOut = std::max(gridZ - 1, 0); maxZOut = std::min(gridZ + 1, 3);
}

PolyVox::Vector3DInt32 Program::worldPosToVoxelIndex(const glm::vec3& worldPos) const
{
    return {
        (int)(worldPos.x * voxelVolumeResolutionScale),
        (int)(worldPos.y * voxelVolumeResolutionScale),
        (int)(worldPos.z * voxelVolumeResolutionScale)
    };
}

void Program::fillVoxelVolume()
{
    const PolyVox::Region volumeRegion = voxelVolume.getEnclosingRegion();
    const PolyVox::Vector3DInt32& lowerCorner = volumeRegion.getLowerCorner();
    const PolyVox::Vector3DInt32& upperCorner = volumeRegion.getUpperCorner();

    for (int z = lowerCorner.getZ(); z <= upperCorner.getZ(); z++)
        for (int y = lowerCorner.getY(); y <= upperCorner.getY(); y++)
            for (int x = lowerCorner.getX(); x <= upperCorner.getX(); x++)
                voxelVolume.setVoxelAt(x, y, z, calcDensity({ x, y, z }));

    /* DEBUG: Only fill voxels that contain an actual particle */
    /*for (int z = lowerCorner.getZ(); z <= upperCorner.getZ(); z++)
        for (int y = lowerCorner.getY(); y <= upperCorner.getY(); y++)
            for (int x = lowerCorner.getX(); x <= upperCorner.getX(); x++)
                voxelVolume.setVoxelAt(x, y, z, 0);
    for (const auto& pos : r)
        voxelVolume.setVoxelAt(worldPosToVoxelIndex(pos), 100.0f);
    /* END DEBUG */
}
