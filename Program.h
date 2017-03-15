#pragma once

#include "Camera.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>

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
    // positions array, for now, contains (x, y, z) coordinates for a cube with sides of size cubeSize
    static const int cubeSize = 7;
    static const int particleCount = cubeSize * cubeSize * cubeSize;
    /* END DEBUG */

    // Radius of influence
    static const float h;
    // Gas constant
    static const float k;
    // Rest density
    static const float rho0;
    // Mass of each particle
    static const float m;
    // Fluid viscosity
    static const float mu;

    // Calculates the density at the given position.
    float getDensity(const glm::vec3& position) const;

    // Calculates the pressure for the given density.
    float getPressure(float density) const;

    // Particle positions
    glm::vec3 r[particleCount];
    // Particle velocities
    glm::vec3 v[particleCount];
};
