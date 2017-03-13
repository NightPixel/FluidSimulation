#include "Program.h"
#include "OpenGLUtils.h"
#include "Randomizer.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

Program::Program()
{
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

    uniModel = glGetUniformLocation(shaderProgram, "model");

    // Set up projection
    glm::mat4 view = glm::lookAt(
        glm::vec3(3.5f, 3.5f, 2.0f),
        glm::vec3(0.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 1.0f)
    );
    GLint uniView = glGetUniformLocation(shaderProgram, "view");
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

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
    // Rotate cube
    rotationAngle += 0.01f * dt * glm::radians(15.0f);
    glm::mat4 model;
    model = glm::rotate(
        model,
        rotationAngle,
        glm::vec3(0.0f, 0.0f, 1.0f)
    );
    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

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
