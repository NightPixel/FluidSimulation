#include "Program.h"
#include "OpenGLUtils.h"
#include "Randomizer.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

Program::Program(GLFWwindow* window)
    : window(window)
{
    glfwGetWindowSize(window, &windowSizeX, &windowSizeY);

    // Initialize OpenGL
    glEnable(GL_DEPTH_TEST);
    glPointSize(5.0f);

    // Create Vertex Array Object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Create cube vertices
    for (int z = -cubeSize / 2; z <= cubeSize / 2; ++z)
        for (int y = -cubeSize / 2; y <= cubeSize / 2; ++y)
            for (int x = -cubeSize / 2; x <= cubeSize / 2; ++x)
            {
                printf("[%i, %i, %i] = (%f, %f, %f)\n",
                    3 * (cubeSize * cubeSize * (z + cubeSize / 2) + cubeSize * (y + cubeSize / 2) + (x + cubeSize / 2)) + 0,
                    3 * (cubeSize * cubeSize * (z + cubeSize / 2) + cubeSize * (y + cubeSize / 2) + (x + cubeSize / 2)) + 1,
                    3 * (cubeSize * cubeSize * (z + cubeSize / 2) + cubeSize * (y + cubeSize / 2) + (x + cubeSize / 2)) + 2,
                    0.3f * x, 0.3f * y, 0.3f * z);
                vertices[3 * (cubeSize * cubeSize * (z + cubeSize / 2) + cubeSize * (y + cubeSize / 2) + (x + cubeSize / 2)) + 0] = 0.3f * x;
                vertices[3 * (cubeSize * cubeSize * (z + cubeSize / 2) + cubeSize * (y + cubeSize / 2) + (x + cubeSize / 2)) + 1] = 0.3f * y;
                vertices[3 * (cubeSize * cubeSize * (z + cubeSize / 2) + cubeSize * (y + cubeSize / 2) + (x + cubeSize / 2)) + 2] = 0.3f * z;
            }

    // Create a Vertex Buffer Object and copy the vertex data to it
    glGenBuffers(1, &vbo);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STREAM_DRAW);

    // Create and compile the vertex shader
    vertexShader = createShaderFromSource("Simple.vert", GL_VERTEX_SHADER);
    auto shaderInfo = checkShaderCompilation(vertexShader);
    if (!shaderInfo.first)
        printf("Vertex shader failed to compile!\n%s\n", shaderInfo.second.c_str());
    else if (!shaderInfo.second.empty())
        printf("Vertex shader compiled with warnings.\n%s\n", shaderInfo.second.c_str());

    // Create and compile the fragment shader
    fragmentShader = createShaderFromSource("Simple.frag", GL_FRAGMENT_SHADER);
    shaderInfo = checkShaderCompilation(fragmentShader);
    if (!shaderInfo.first)
        printf("Vertex shader failed to compile!\n%s\n", shaderInfo.second.c_str());
    else if (!shaderInfo.second.empty())
        printf("Vertex shader compiled with warnings.\n%s\n", shaderInfo.second.c_str());

    // Link the vertex and fragment shader into a shader program
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glBindFragDataLocation(shaderProgram, 0, "outColor");
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    // Specify the layout of the vertex data
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 0, 0);

    // Set up model, view, projection matrices
    glm::mat4 model{}; // Identity matrix
    GLint uniModel = glGetUniformLocation(shaderProgram, "model");
    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

    uniView = glGetUniformLocation(shaderProgram, "view");

    glm::mat4 proj = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 1.0f, 10.0f);
    GLint uniProj = glGetUniformLocation(shaderProgram, "proj");
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));
}

Program::~Program()
{
    glDeleteProgram(shaderProgram);
    glDeleteShader(fragmentShader);
    glDeleteShader(vertexShader);

    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}

// Updates the state of the program, and draws a new frame.
// The dt parameter is the elapsed frame time since last frame, in ms.
void Program::update(float dt)
{
    glm::mat4 view = camera.getViewMatrix();
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

    // Randomize vertices
    for (auto& vert : vertices)
        vert += Randomizer::random(-0.001f, 0.001f);

    // "orphan" vertices array; we no longer need it
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), nullptr, GL_STREAM_DRAW);
    // Upload new vertex data
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);

    // Draw cube
    glDrawArrays(GL_POINTS, 0, cubeSize * cubeSize * cubeSize);
}

// Called when the mouse cursor is moved.
void Program::onMouseMoved(float dxPos, float dyPos)
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

// Called when the mouse wheel is scrolled.
void Program::onMouseScrolled(float yOffset)
{
    camera.zoom(yOffset / 30.0f);
}
