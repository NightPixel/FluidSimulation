#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "Program.h"

// Callback function called by GLFW when the cursor position changes.
void cursorPositionCallback(GLFWwindow* window, double xPos, double yPos)
{
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
    Program* program = static_cast<Program*>(glfwGetWindowUserPointer(window));

    program->onMouseScrolled((float)yOffset);
}

int main()
{
    glfwInit();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    GLFWwindow* window = glfwCreateWindow(1280, 800, "Fluid simulation", nullptr, nullptr); // Windowed
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    glewInit();

    Program fluidSimulationProgram(window);
    glfwSetWindowUserPointer(window, &fluidSimulationProgram);
    glfwSetCursorPosCallback(window, cursorPositionCallback);
    glfwSetScrollCallback(window, scrollCallback);

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

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
