#include "ProgramBase.h"
#include "OpenGLUtils.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <PolyVoxCore/VertexTypes.h>
#include <tuple>

ProgramBase::ProgramBase(GLFWwindow* window)
    : window(window)
{
    glfwGetWindowSize(window, &windowSizeX, &windowSizeY);
    antTweakBar = TwNewBar("Simulation settings");
    TwDefine("GLOBAL fontsize=3");

    // Initialize OpenGL
    glEnable(GL_DEPTH_TEST);
    glPointSize(5.0f);

    std::tie(simpleVertexShader, simpleFragmentShader, simpleShaderProgram) = createShaderProgram("Simple.vert", "Simple.frag", { { 0, "outColor" } });
    glUseProgram(simpleShaderProgram);

    // Set up model, view, projection matrices
    glUniformMatrix4fv(glGetUniformLocation(simpleShaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(glm::mat4{}));  // Identity matrix
    simpleViewUniform = glGetUniformLocation(simpleShaderProgram, "view");
    glUniformMatrix4fv(glGetUniformLocation(simpleShaderProgram, "proj"), 1, GL_FALSE,
        glm::value_ptr(glm::perspective(glm::radians(45.0f), (float)windowSizeX / windowSizeY, 1.0f, 25.0f))
    );

    // Create a Vertex Array Object for the particle points
    glGenVertexArrays(1, &pointsVAO);
    glBindVertexArray(pointsVAO);
    // Create a Vertex Buffer Object for the particle points
    glGenBuffers(1, &pointsVBO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);

    // Specify the layout of the particle points vertex data
    // glm::vec3 layout:
    //    (x, y, z)-position (3 floats)
    GLint pointsPosAttrib = glGetAttribLocation(simpleShaderProgram, "position");
    glEnableVertexAttribArray(pointsPosAttrib);
    glVertexAttribPointer(pointsPosAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), reinterpret_cast<void*>(0 * sizeof(float)));

    std::tie(waterVertexShader, waterFragmentShader, waterShaderProgram) = createShaderProgram("Water.vert", "Water.frag", { { 0, "outColor" } });
    glUseProgram(waterShaderProgram);

    // Set up model, view, projection matrices
    glUniformMatrix4fv(glGetUniformLocation(waterShaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(glm::mat4{}));  // Identity matrix
    waterViewUniform = glGetUniformLocation(waterShaderProgram, "view");
    glUniformMatrix4fv(glGetUniformLocation(waterShaderProgram, "proj"), 1, GL_FALSE,
        glm::value_ptr(glm::perspective(glm::radians(45.0f), (float)windowSizeX / windowSizeY, 1.0f, 25.0f))
    );

    // Set up fragment shader uniforms
    waterCamUniform = glGetUniformLocation(waterShaderProgram, "camPos");
    glUniform3fv(glGetUniformLocation(waterShaderProgram, "lightPos"), 1, glm::value_ptr(glm::vec3{ 3.0f, 3.0f, 3.0f }));

    glUniform3fv(glGetUniformLocation(waterShaderProgram, "ambientSceneColor"), 1, glm::value_ptr(glm::vec3{ 0.5f, 0.5f, 0.5f }));
    glUniform3fv(glGetUniformLocation(waterShaderProgram, "ambientLightColor"), 1, glm::value_ptr(glm::vec3{ 0.1f, 0.1f, 0.1f }));

    glUniform3fv(glGetUniformLocation(waterShaderProgram, "diffuseMaterialColor"), 1, glm::value_ptr(glm::vec3{ 0.5f, 0.5f, 0.95f }));
    glUniform3fv(glGetUniformLocation(waterShaderProgram, "diffuseLightColor"), 1, glm::value_ptr(glm::vec3{ 1.0f, 1.0f, 1.0f }));

    glUniform3fv(glGetUniformLocation(waterShaderProgram, "specularMaterialColor"), 1, glm::value_ptr(glm::vec3{ 1.0f, 1.0f, 1.0f }));
    glUniform3fv(glGetUniformLocation(waterShaderProgram, "specularLightColor"), 1, glm::value_ptr(glm::vec3{ 1.0f, 1.0f, 1.0f }));
    glUniform1i(glGetUniformLocation(waterShaderProgram, "shininess"), 32);

    // Create a Vertex Array Object for the surface mesh
    glGenVertexArrays(1, &meshVAO);
    glBindVertexArray(meshVAO);
    // Create a Vertex Buffer Object for the surface mesh
    glGenBuffers(1, &meshVBO);
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    // Create an Element Buffer Object for the surface mesh
    glGenBuffers(1, &meshEBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshEBO);

    // Specify the layout of the surface mesh vertex data
    // PolyVox::PositionMaterialNormal layout:
    //    (x, y, z)-position (3 floats)
    //    (x, y, z)-normal   (3 floats)
    //    material           (1 float)
    GLint meshPosAttrib = glGetAttribLocation(waterShaderProgram, "position");
    glEnableVertexAttribArray(meshPosAttrib);
    glVertexAttribPointer(meshPosAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(0 * sizeof(float)));
    GLint meshNormAttrib = glGetAttribLocation(waterShaderProgram, "normal");
    glEnableVertexAttribArray(meshNormAttrib);
    glVertexAttribPointer(meshNormAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(3 * sizeof(float)));
    GLint meshMatAttrib = glGetAttribLocation(waterShaderProgram, "material");
    glEnableVertexAttribArray(meshMatAttrib);
    glVertexAttribPointer(meshMatAttrib, 1, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(6 * sizeof(float)));
}

ProgramBase::~ProgramBase()
{
    glDeleteProgram(simpleShaderProgram);
    glDeleteShader(simpleFragmentShader);
    glDeleteShader(simpleVertexShader);
    glDeleteProgram(waterShaderProgram);
    glDeleteShader(waterFragmentShader);
    glDeleteShader(waterVertexShader);

    glDeleteBuffers(1, &meshEBO);
    glDeleteBuffers(1, &meshVBO);
    glDeleteVertexArrays(1, &meshVAO);
    glDeleteBuffers(1, &pointsVBO);
    glDeleteVertexArrays(1, &pointsVAO);
}

void ProgramBase::onMouseMoved(float dxPos, float dyPos)
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

void ProgramBase::onMouseScrolled(float yOffset)
{
    camera.zoom(yOffset / 30.0f);
}

void ProgramBase::onKeypress(int key, int action)
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
}
