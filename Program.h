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

    // vertices, for now, contains (x, y, z) coordinates for a cube with sides of size cubeSize
    static const int cubeSize = 7;
    GLfloat vertices[cubeSize * cubeSize * cubeSize * 3];
};
