#include "FluidVisualizer.h"
#include "OpenGLUtils.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

FluidVisualizer::FluidVisualizer(GLFWwindow* window)
    : FluidSimulator(window),
    voxelVolume({
        worldPosToVoxelIndex(minPos) - PolyVox::Vector3DInt32{
            (int)std::ceil(h * voxelVolumeResolutionScale),
            (int)std::ceil(h * voxelVolumeResolutionScale),
            (int)std::ceil(h * voxelVolumeResolutionScale)
        },
        worldPosToVoxelIndex(maxPos) + PolyVox::Vector3DInt32{
            (int)std::ceil(h * voxelVolumeResolutionScale),
            (int)std::ceil(h * voxelVolumeResolutionScale),
            (int)std::ceil(h * voxelVolumeResolutionScale)
        },
    }),
    surfaceExtractor(&voxelVolume, voxelVolume.getEnclosingRegion(), &surfaceMesh)
{
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

    // Set up fragment shader uniforms
    glUniform3fv(glGetUniformLocation(simpleShaderProgram, "minWorldPos"), 1, glm::value_ptr(minPos));
    glUniform3fv(glGetUniformLocation(simpleShaderProgram, "maxWorldPos"), 1, glm::value_ptr(maxPos));

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

FluidVisualizer::~FluidVisualizer()
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

// Draws a new frame.
void FluidVisualizer::draw()
{
    // Run marching cubes using PolyVox, and retrieve the vertex and index buffers
    fillVoxelVolume();
    surfaceExtractor.execute();
    const auto& lowerCorner = voxelVolume.getEnclosingRegion().getLowerCorner();
    surfaceMesh.translateVertices({ (float)lowerCorner.getX(), (float)lowerCorner.getY(), (float)lowerCorner.getZ() });
    surfaceMesh.scaleVertices(1.0f / voxelVolumeResolutionScale);
    const std::vector<uint32_t>& indices = surfaceMesh.getIndices();
    std::vector<PolyVox::PositionMaterialNormal>& vertices = surfaceMesh.getRawVertexData();
    for (auto& vert : vertices) // Clamp vertex locations to world boundaries
    {
        vert.position.setElements(
            clamp(minPos.x + gridOffset.x, maxPos.x + gridOffset.x, vert.position.getX() + gridOffset.x),
            clamp(minPos.y + gridOffset.y, maxPos.y + gridOffset.y, vert.position.getY() + gridOffset.y),
            clamp(minPos.z + gridOffset.z, maxPos.z + gridOffset.z, vert.position.getZ() + gridOffset.z)
        );
    }

    const glm::mat4 view = camera.getViewMatrix();
    glUseProgram(waterShaderProgram);
    glUniformMatrix4fv(waterViewUniform, 1, GL_FALSE, glm::value_ptr(view));
    glUniform3fv(waterCamUniform, 1, glm::value_ptr(camera.getPosition()));

    // Draw surface mesh
    glBindVertexArray(meshVAO);
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(PolyVox::PositionMaterialNormal), vertices.data(), GL_STREAM_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(uint32_t), indices.data(), GL_STREAM_DRAW);
    glDrawElements(GL_TRIANGLES, (GLsizei)indices.size(), GL_UNSIGNED_INT, nullptr);

    // Draw obstacle
    glBufferData(GL_ARRAY_BUFFER, objectVertices.size() * sizeof(PolyVox::PositionMaterialNormal), objectVertices.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)objectVertices.size());

    glUseProgram(simpleShaderProgram);
    glUniformMatrix4fv(simpleViewUniform, 1, GL_FALSE, glm::value_ptr(view));

    // Draw points
    glBindVertexArray(pointsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(r), r, GL_STREAM_DRAW);
    glDrawArrays(GL_POINTS, 0, (GLsizei)std::size(r));

    // Draw world bounds
    glBufferData(GL_ARRAY_BUFFER, sizeof(worldBoundsVertices), worldBoundsVertices, GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)std::size(worldBoundsVertices));
}

PolyVox::Vector3DInt32 FluidVisualizer::worldPosToVoxelIndex(const glm::vec3& worldPos) const
{
    return {
        (int)((worldPos.x - gridOffset.x) * voxelVolumeResolutionScale),
        (int)((worldPos.y - gridOffset.y) * voxelVolumeResolutionScale),
        (int)((worldPos.z - gridOffset.z) * voxelVolumeResolutionScale)
    };
}

glm::vec3 FluidVisualizer::voxelIndexToWorldPos(int voxelX, int voxelY, int voxelZ) const
{
    static const float invVoxelVolumeResolutionScale = 1.0f / voxelVolumeResolutionScale;

    return {
        voxelX * invVoxelVolumeResolutionScale + gridOffset.x,
        voxelY * invVoxelVolumeResolutionScale + gridOffset.y,
        voxelZ * invVoxelVolumeResolutionScale + gridOffset.z,
    };
}

void FluidVisualizer::fillVoxelVolume()
{
    const PolyVox::Region volumeRegion = voxelVolume.getEnclosingRegion();
    const PolyVox::Vector3DInt32& lowerCorner = volumeRegion.getLowerCorner();
    const PolyVox::Vector3DInt32& upperCorner = volumeRegion.getUpperCorner();

#pragma omp parallel for
    for (int z = lowerCorner.getZ(); z <= upperCorner.getZ(); z++)
        for (int y = lowerCorner.getY(); y <= upperCorner.getY(); y++)
            for (int x = lowerCorner.getX(); x <= upperCorner.getX(); x++)
                voxelVolume.setVoxelAt(x, y, z, calcDensity(voxelIndexToWorldPos(x, y, z)));
}
