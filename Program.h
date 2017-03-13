#pragma once

#include <GL/glew.h>

class Program
{
public:
    Program();
    ~Program();

    // Updates the state of the program, and draws a new frame.
    // The dt parameter is the elapsed frame time since last frame, in ms.
    void update(float dt);

private:
    GLuint vao;
    GLuint vbo;
    GLuint vertexShader;
    GLuint fragmentShader;
    GLuint shaderProgram;
    GLint uniModel;
    float rotationAngle = 0.0f;

    // vertices, for now, contains (x, y, z) coordinates for a cube with sides of size cubeSize
    static const int cubeSize = 7;
    GLfloat vertices[cubeSize * cubeSize * cubeSize * 3];
};
