#pragma once

#include "Camera.h"
#include "Geometry.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <AntTweakBar.h>
#include <vector>

class FluidBase
{
public:
    explicit FluidBase(GLFWwindow* window);

    void addModelsFromOBJFile(const std::string& fileName,
        const glm::vec3& offset = {}, const glm::vec3& rotation = {}, const glm::vec3& scale = glm::vec3{1.0f, 1.0f, 1.0f});

    // Called when the mouse cursor is moved.
    void onMouseMoved(float dxPos, float dyPos);
    // Called when the mouse wheel is scrolled.
    void onMouseScrolled(float yOffset);
    // Called when a key is pressed or released.
    void onKeypress(int key, int action);

protected:
    Camera camera;
    GLFWwindow* window;
    TwBar* antTweakBar;
    int windowSizeX;
    int windowSizeY;

    bool holdForward = false;
    bool holdBackward = false;
    bool holdRight = false;
    bool holdLeft = false;
    bool holdUp = false;
    bool holdDown = false;
    bool holdShift = false;

    std::vector<Model> models;
};
