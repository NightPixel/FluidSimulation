#include "FluidBase.h"

FluidBase::FluidBase(GLFWwindow* window)
    : window(window)
{
    glfwGetWindowSize(window, &windowSizeX, &windowSizeY);
    antTweakBar = TwNewBar("Simulation settings");
    TwDefine("GLOBAL fontsize=3");

    //addModelsFromOBJFile("cube.obj", glm::vec3{-1.0f, -1.0f, -1.0f}, glm::vec3{45.0f, 75.0f, 0.0f}, glm::vec3{1.5f, 0.5f, 0.5f});
    //addModelsFromOBJFile("cube.obj", glm::vec3{0.3f, 1.0f, 0.3f}, glm::vec3{-25.0f, -30.0f, 100.0f}, glm::vec3{0.2f, 0.3f, 0.4f});
    //addModelsFromOBJFile("cube.obj", glm::vec3{0.0f, -1.75f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{3.0f, 0.2f, 3.0f});
    //addModelsFromOBJFile("lowresbunny.obj", glm::vec3{0.0f, -3.3f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.75f, 0.75f, 0.75f});
}

void FluidBase::addModelsFromOBJFile(const std::string& fileName, const glm::vec3& offset, const glm::vec3& rotation, const glm::vec3& scale)
{
    auto newModels = loadOBJFile(fileName, offset, rotation, scale);
    models.insert(models.end(), newModels.begin(), newModels.end());
}

void FluidBase::loadScene(int sceneNumber)
{
    models.clear();
    switch (sceneNumber)
    {
    case 0:
        break;
    case 1:
        addModelsFromOBJFile("cube.obj", glm::vec3{ 0.0f, -1.75f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.0f }, glm::vec3{ 3.0f, 0.2f, 3.0f });
        break;
    case 2:
        addModelsFromOBJFile("cube.obj", glm::vec3{ 0.0f, 0.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 90.0f }, glm::vec3{ 3.0f, 0.2f, 3.0f });
        break;
    case 3:
        addModelsFromOBJFile("cube.obj", glm::vec3{ -1.5f, -1.75f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.0f }, glm::vec3{ 3.0f, 0.2f, 3.0f });
        break;
    case 4:
        addModelsFromOBJFile("cube.obj", glm::vec3{ 1.5f, -1.0f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.0f }, glm::vec3{ 3.0f, 0.2f, 3.0f });
        addModelsFromOBJFile("cube.obj", glm::vec3{ -1.5f, -1.75f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.0f }, glm::vec3{ 3.0f, 0.2f, 3.0f });
        break;
    case 5:
        addModelsFromOBJFile("lowresbunny.obj", glm::vec3{ 0.0f, -3.3f, 0.0f }, glm::vec3{ 0.0f, 0.0f, 0.0f }, glm::vec3{ 0.75f, 0.75f, 0.75f });
        break;
    }
}

void FluidBase::onMouseMoved(float dxPos, float dyPos)
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

void FluidBase::onMouseScrolled(float yOffset)
{
    camera.zoom(yOffset / 30.0f);
}

void FluidBase::onKeypress(int key, int action)
{
    printf("K: %i A: %i\n", key, action);
    switch (key)
    {
    case 'W':
        holdForward = action != GLFW_RELEASE;
        break;
    case 'S':
        holdBackward = action != GLFW_RELEASE;
        break;
    case 'A':
        holdLeft = action != GLFW_RELEASE;
        break;
    case 'D':
        holdRight = action != GLFW_RELEASE;
        break;
    case 32: // Space
        holdUp = action != GLFW_RELEASE;
        break;
    case 341: // Left control
        holdDown = action != GLFW_RELEASE;
        break;
    case 340: // Left shift
        holdShift = action != GLFW_RELEASE;
        break;
    }
    // Scene selection using F keys
    if (key >= 290 && key <= 295)
    {
        loadScene(key - 290);
    }
}
