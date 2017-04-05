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

    std::tie(simpleVertexShader, simpleFragmentShader, simpleShaderProgram) = createShaderProgram("Simple.vert", "Simple.frag", {{ 0, "outColor" }});
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

    std::tie(phongVertexShader, phongFragmentShader, phongShaderProgram) = createShaderProgram("Phong.vert", "Phong.frag", {{ 0, "outColor" }});
    glUseProgram(phongShaderProgram);

    // Set up model, view, projection matrices
    glUniformMatrix4fv(glGetUniformLocation(phongShaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(glm::mat4{}));  // Identity matrix
    phongViewUniform = glGetUniformLocation(phongShaderProgram, "view");
    glUniformMatrix4fv(glGetUniformLocation(phongShaderProgram, "proj"), 1, GL_FALSE,
        glm::value_ptr(glm::perspective(glm::radians(45.0f), (float)windowSizeX / windowSizeY, 1.0f, 25.0f))
    );

    // Set up fragment shader uniforms
    phongCamUniform = glGetUniformLocation(phongShaderProgram, "camPos");
    glUniform3fv(glGetUniformLocation(phongShaderProgram, "lightPos"), 1, glm::value_ptr(glm::vec3{3.0f, 3.0f, 3.0f}));

    phongAmbientUniform = glGetUniformLocation(phongShaderProgram, "ambientMaterialColor");
    glUniform3fv(glGetUniformLocation(phongShaderProgram, "ambientLightColor"), 1, glm::value_ptr(glm::vec3{0.1f, 0.1f, 0.1f}));

    phongDiffuseUniform = glGetUniformLocation(phongShaderProgram, "diffuseMaterialColor");
    glUniform3fv(glGetUniformLocation(phongShaderProgram, "diffuseLightColor"), 1, glm::value_ptr(glm::vec3{1.0f, 1.0f, 1.0f}));

    phongSpecularUniform = glGetUniformLocation(phongShaderProgram, "specularMaterialColor");
    glUniform3fv(glGetUniformLocation(phongShaderProgram, "specularLightColor"), 1, glm::value_ptr(glm::vec3{1.0f, 1.0f, 1.0f}));

    phongShininessUniform = glGetUniformLocation(phongShaderProgram, "shininess");

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
    GLint meshPosAttrib = glGetAttribLocation(phongShaderProgram, "position");
    glEnableVertexAttribArray(meshPosAttrib);
    glVertexAttribPointer(meshPosAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(0 * sizeof(float)));
    GLint meshNormAttrib = glGetAttribLocation(phongShaderProgram, "normal");
    glEnableVertexAttribArray(meshNormAttrib);
    glVertexAttribPointer(meshNormAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(3 * sizeof(float)));
    GLint meshMatAttrib = glGetAttribLocation(phongShaderProgram, "material");
    glEnableVertexAttribArray(meshMatAttrib);
    glVertexAttribPointer(meshMatAttrib, 1, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(6 * sizeof(float)));
}

FluidVisualizer::~FluidVisualizer()
{
    glDeleteProgram(simpleShaderProgram);
    glDeleteShader(simpleFragmentShader);
    glDeleteShader(simpleVertexShader);
    glDeleteProgram(phongShaderProgram);
    glDeleteShader(phongFragmentShader);
    glDeleteShader(phongVertexShader);

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
    surfaceMesh.translateVertices({(float)lowerCorner.getX(), (float)lowerCorner.getY(), (float)lowerCorner.getZ()});
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
    glUseProgram(phongShaderProgram);
    glUniformMatrix4fv(phongViewUniform, 1, GL_FALSE, glm::value_ptr(view));
    glUniform3fv(phongCamUniform, 1, glm::value_ptr(camera.getPosition()));


    glUniform3fv(phongAmbientUniform, 1, glm::value_ptr(glm::vec3{0.5f, 0.5f, 0.5f}));
    glUniform3fv(phongDiffuseUniform, 1, glm::value_ptr(glm::vec3{0.5f, 0.5f, 0.95f}));
    glUniform3fv(phongSpecularUniform, 1, glm::value_ptr(glm::vec3{1.0f, 1.0f, 1.0f}));
    glUniform1f(phongShininessUniform, 32.0f);

    // Draw surface mesh
    glBindVertexArray(meshVAO);
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(PolyVox::PositionMaterialNormal), vertices.data(), GL_STREAM_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(uint32_t), indices.data(), GL_STREAM_DRAW);
    glDrawElements(GL_TRIANGLES, (GLsizei)indices.size(), GL_UNSIGNED_INT, nullptr);

    // Draw models
    for (const auto& model : models)
    {
        glUniform3fv(phongAmbientUniform, 1, glm::value_ptr(model.ambientColor));
        glUniform3fv(phongDiffuseUniform, 1, glm::value_ptr(model.diffuseColor));
        glUniform3fv(phongSpecularUniform, 1, glm::value_ptr(model.specularColor));
        glUniform1f(phongShininessUniform, model.specularExponent);
        glBufferData(GL_ARRAY_BUFFER, model.vertexData.size() * sizeof(PolyVox::PositionMaterialNormal), model.vertexData.data(), GL_STREAM_DRAW);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)model.vertexData.size());
    }

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
    return{
        (int)((worldPos.x - gridOffset.x) * voxelVolumeResolutionScale),
        (int)((worldPos.y - gridOffset.y) * voxelVolumeResolutionScale),
        (int)((worldPos.z - gridOffset.z) * voxelVolumeResolutionScale)
    };
}

glm::vec3 FluidVisualizer::voxelIndexToWorldPos(int voxelX, int voxelY, int voxelZ) const
{
    static const float invVoxelVolumeResolutionScale = 1.0f / voxelVolumeResolutionScale;

    return{
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
