#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <AntTweakBar.h>
#include <AntTweakBar_GLFW3.h>
#include "Program.h"

// Callback function called by GLFW when the cursor position changes.
void cursorPositionCallback(GLFWwindow* window, double xPos, double yPos)
{
    if (TwEventMousePosGLFW((int)xPos, (int)yPos))
        return;

    static double prevXPos;
    static double prevYPos;

    Program* program = static_cast<Program*>(glfwGetWindowUserPointer(window));

    program->onMouseMoved((float)(xPos - prevXPos), (float)(yPos - prevYPos));

    prevXPos = xPos;
    prevYPos = yPos;
}

// Callback function called by GLFW when the mouse wheel is scrolled.
void scrollCallback(GLFWwindow* window, double xOffset, double yOffset)
{
    if (TwEventMouseWheelGLFW((int)yOffset))
        return;

    Program* program = static_cast<Program*>(glfwGetWindowUserPointer(window));

    program->onMouseScrolled((float)yOffset);
}

// Callback function called by GLFW when a mouse button is pressed or released.
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    TwEventMouseButtonGLFW(button, action);
}

// Callback function called by GLFW when a keyboard key is pressed, repeated or released.
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    TwEventKeyGLFW(TwConvertKeyGLFW3to2(key), action);
}

// Callback function called by GLFW when a Unicode character is input.
void charCallback(GLFWwindow* window, unsigned int codepoint)
{
    TwEventCharGLFW(codepoint, GLFW_PRESS);
}

int main()
{
    glfwInit();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    const int windowSizeX = 1280;
    const int windowSizeY = 800;
    GLFWwindow* window = glfwCreateWindow(windowSizeX, windowSizeY, "Fluid simulation", nullptr, nullptr); // Windowed
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    glewInit();

    TwInit(TW_OPENGL_CORE, nullptr);
    TwWindowSize(windowSizeX, windowSizeY);

    Program fluidSimulationProgram(window);
    glfwSetWindowUserPointer(window, &fluidSimulationProgram);
    glfwSetCursorPosCallback(window, cursorPositionCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetKeyCallback(window, keyCallback);
    glfwSetCharCallback(window, charCallback);

    auto previousTime = glfwGetTime();
    while (!glfwWindowShouldClose(window))
    {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GL_TRUE);

        // Clear the screen to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        auto currentTime = glfwGetTime();
        fluidSimulationProgram.update((float)(currentTime - previousTime));
        previousTime = currentTime;

        // Draw AntTweakBar
        TwDraw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    TwTerminate();
    glfwTerminate();
    return 0;
}
