#pragma once

#include "Camera.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>

class Program
{
public:
    Program(GLFWwindow* window);
    ~Program();

    // Updates the state of the program, and draws a new frame.
    // The dt parameter is the elapsed frame time since last frame, in ms.
    void update(float dt);

    // Called when the mouse cursor is moved.
    void onMouseMoved(float dxPos, float dyPos);
    // Called when the mouse wheel is scrolled.
    void onMouseScrolled(float yOffset);

private:
    Camera camera;
    GLFWwindow* window;

    GLuint vao;
    GLuint vbo;
    GLuint vertexShader;
    GLuint fragmentShader;
    GLuint shaderProgram;
    GLint uniView;
    int windowSizeX;
    int windowSizeY;

    /* DEBUG */
    // The particle positions array ('r'), for now, contains (x, y, z) coordinates for a cube with sides of size cubeSize
    static const int cubeSize = 7;
    static const int particleCount = cubeSize * cubeSize * cubeSize;
    // For now, all particles reside inside a larger cube with these dimensions
    glm::vec3 minPos = glm::vec3(-2.0f, -2.0f, -2.0f);
    glm::vec3 maxPos = glm::vec3(2.0f, 2.0f, 2.0f);
    /* END DEBUG */

    // Radius of influence
    float h = 0.5f;
    // Gas constant
    float k = 1000.0f;
    // Rest density
    float rho0 = 20.0f;
    // Mass of each particle
    float m = 1.0f;
    // Fluid viscosity
    float mu = 3.0f; // 3.0f
    // Surface tension coefficient
    float sigma = 0.01f; // 0.01f
    // Surface tension is only evaluated if |n| exceeds this threshold
    // (where n is the gradient field of the smoothed color field).
    float csNormThreshold = 1.0f;
    // Gravity acceleration
    glm::vec3 gravity = glm::vec3(0.0f, -9.81f, 0.0f);

    // Particle positions
    glm::vec3 r[particleCount];
    // Particle velocities
    glm::vec3 v[particleCount];

    // Particle grid data structure: changes O(n^2) to O(nm): we don't have to check all other particles, but
    // only particles in adjacent grid cells.
    // Possible TODO: Store entire particle data for better cache usage
    // Possible TODO: Particle grid sizes are 4.0f / 1.0f: make this dependent on actual h and minPos/maxPos with dyn alloc (slow?)
    std::vector<int> particleGrid[4][4][4];

    void fillParticleGrid();

    glm::vec3 calcPressureForce(int particleId, float* rho, float* p);
    glm::vec3 calcViscosityForce(int particleId, float* rho);
    glm::vec3 calcSurfaceForce(int particleId, float* rho);

    void calcDensity(int particleId, float& rho);

    void getAdjacentCells(int gridX, int gridY, int gridZ,
        int& minXOut, int& maxXOut,
        int& minYOut, int& maxYOut,
        int& minZOut, int& maxZOut);
};
