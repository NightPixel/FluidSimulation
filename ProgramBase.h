#pragma once

#include "Camera.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <AntTweakBar.h>

class ProgramBase
{
public:
    explicit ProgramBase(GLFWwindow* window);

    // Called when the mouse cursor is moved.
    void onMouseMoved(float dxPos, float dyPos);
    // Called when the mouse wheel is scrolled.
    void onMouseScrolled(float yOffset);
    // Called when a key is pressed or released.
    void onKeypress(int key, int action);

protected:
    ~ProgramBase();

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

    bool holdForward = false;
    bool holdBackward = false;
    bool holdRight = false;
    bool holdLeft = false;
    bool holdUp = false;
    bool holdDown = false;
    bool holdShift = false;
};
