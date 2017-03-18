#include "Program.h"
#include "Kernels.h"
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

    // Create cube positions
    for (int z = 0; z != cubeSize; ++z)
        for (int y = 0; y != cubeSize; ++y)
            for (int x = 0; x != cubeSize; ++x)
            {
                r[z * cubeSize * cubeSize + y * cubeSize + x] =
                    glm::vec3(0.3f * (x - cubeSize / 2), 0.3f * (y - cubeSize / 2), 0.3f * (z - cubeSize / 2));
                printf("[%i] = (%f, %f, %f)\n",
                    z * cubeSize * cubeSize + y * cubeSize + x,
                    0.3f * (x - cubeSize / 2), 0.3f * (y - cubeSize / 2), 0.3f * (z - cubeSize / 2));
            }

    // Create a Vertex Buffer Object and copy the vertex data to it
    glGenBuffers(1, &vbo);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(r), r, GL_STREAM_DRAW);

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
    glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    // Set up model, view, projection matrices
    glm::mat4 model{}; // Identity matrix
    GLint uniModel = glGetUniformLocation(shaderProgram, "model");
    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

    uniView = glGetUniformLocation(shaderProgram, "view");

    glm::mat4 proj = glm::perspective(glm::radians(45.0f), (float)windowSizeX / windowSizeY, 1.0f, 10.0f);
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
    /* DEBUG: Randomize positions to test CPU-to-GPU vertex streaming */
    //for (auto& vert : r)
    //    vert += Randomizer::random(-0.001f, 0.001f);
    /* END DEBUG */

    float rho[particleCount] = {}; // Particle densities
    float p[particleCount] = {}; // Particle pressure
    // First, calculate density and pressure at each particle position
    for (size_t i = 0; i != particleCount; ++i)
    {
        for (size_t j = 0; j != particleCount; ++j)
            rho[i] += m * poly6(r[i] - r[j], h);
        p[i] = k * (rho[i] - rho0);
    }

    // Calculate the force acting on each particle
    glm::vec3 forces[particleCount];
    for (size_t i = 0; i != particleCount; ++i)
    {
        glm::vec3 pressureForce;
        for (size_t j = 0; j != particleCount; ++j)
            pressureForce += -m * ((p[i] + p[j]) / (2 * rho[j])) * spikyGradient(r[i] - r[j], h);

        glm::vec3 viscosityForce;
        for (size_t j = 0; j != particleCount; ++j)
            viscosityForce += mu * m * ((v[j] - v[i]) / rho[j]) * viscosityLaplacian(r[i] - r[j], h);

        glm::vec3 surfaceForce;
        glm::vec3 n; // Gradient field of the color field
        for (size_t j = 0; j != particleCount; ++j)
            n += m * (1 / rho[j]) * poly6Gradient(r[i] - r[j], h);
        if (glm::length(n) > csNormThreshold)
        {
            float csLaplacian = 0.0f; // Laplacian of the color field
            for (size_t j = 0; j != particleCount; ++j)
                csLaplacian += m * (1 / rho[j]) * poly6Laplacian(r[i] - r[j], h);
            surfaceForce = -sigma * csLaplacian * glm::normalize(n);
        }

        glm::vec3 gravityForce = gravity * rho[i];

        forces[i] = pressureForce + viscosityForce + surfaceForce + gravityForce;
    }

    // TODO: add external forces (for example, forces created by user input)

    // Move each particle
    for (size_t i = 0; i != particleCount; ++i)
    {
        glm::vec3 a = forces[i] / rho[i]; // Acceleration
        // Semi-implicit Euler integration (TODO: better integration?)
        v[i] += a * dt;
        r[i] += v[i] * dt;

        // Rudimentary collision; the particles reside inside a hard-coded AABB
        // Upon collision with bounds, push particles out of objects, and reflect their velocity vector
        for (int dim = 0; dim != 3; ++dim) // Loop over x, y and z components
        {
            if (r[i][dim] < minPos[dim])
            {
                r[i][dim] = minPos[dim];
                v[i][dim] = -v[i][dim];
            }
            else if (r[i][dim] > maxPos[dim])
            {
                r[i][dim] = maxPos[dim];
                v[i][dim] = -v[i][dim];
            }
        }
    }

    // "orphan" positions array; we no longer need it
    glBufferData(GL_ARRAY_BUFFER, sizeof(r), nullptr, GL_STREAM_DRAW);
    // Upload new vertex data
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(r), r);

    // Update view matrix
    glm::mat4 view = camera.getViewMatrix();
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

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
