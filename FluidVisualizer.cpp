#include "FluidVisualizer.h"
#include "OpenGLUtils.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>

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
    // Allow for toggling mesh construction
    TwAddVarRW(antTweakBar, "Construct mesh", TW_TYPE_BOOLCPP, &meshConstruction, "");

    // Initialize OpenGL
    glEnable(GL_DEPTH_TEST);
    glPointSize(5.0f);

    // Compile vertex and pixel Simple shader
    std::tie(simpleVertexShader, simpleFragmentShader, simpleShaderProgram) = createShaderProgram("Simple.vert", "Simple.frag", {{ 0, "outColor" }});
    glUseProgram(simpleShaderProgram);

    // Set up model, view, projection matrices
    simpleModelUniform = glGetUniformLocation(simpleShaderProgram, "model");
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
    GLint simpleShaderPosAttrib = glGetAttribLocation(simpleShaderProgram, "position");
    glEnableVertexAttribArray(simpleShaderPosAttrib);
    glVertexAttribPointer(simpleShaderPosAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), reinterpret_cast<void*>(0 * sizeof(float)));

    // Create a Vertex Array Object for the world bounds
    glGenVertexArrays(1, &boundsVAO);
    glBindVertexArray(boundsVAO);
    // Create a Vertex Buffer Object for the world bounds
    glGenBuffers(1, &boundsVBO);
    glBindBuffer(GL_ARRAY_BUFFER, boundsVBO);
    // Upload vertex data
    glBufferData(GL_ARRAY_BUFFER, sizeof(worldBoundsVertices), worldBoundsVertices, GL_STATIC_DRAW);

    // Specify the layout of the world bounds vertex data
    // glm::vec3 layout:
    //    (x, y, z)-position (3 floats)
    glEnableVertexAttribArray(simpleShaderPosAttrib);
    glVertexAttribPointer(simpleShaderPosAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), reinterpret_cast<void*>(0 * sizeof(float)));

    // Compile vertex and pixel Phong shader
    std::tie(phongVertexShader, phongFragmentShader, phongShaderProgram) = createShaderProgram("Phong.vert", "Phong.frag", {{ 0, "outColor" }});
    glUseProgram(phongShaderProgram);

    // Set up model, view, projection, normal matrices
    phongModelUniform = glGetUniformLocation(phongShaderProgram, "model");
    phongViewUniform = glGetUniformLocation(phongShaderProgram, "view");
    glUniformMatrix4fv(glGetUniformLocation(phongShaderProgram, "proj"), 1, GL_FALSE,
        glm::value_ptr(glm::perspective(glm::radians(45.0f), (float)windowSizeX / windowSizeY, 1.0f, 25.0f))
    );
    phongNormalMatUniform = glGetUniformLocation(phongShaderProgram, "normalMat");

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
    phongAlphaUniform = glGetUniformLocation(phongShaderProgram, "alpha");

    // Setup buffers for the surface mesh and scene meshes
    setupMeshBuffers();
}

void FluidVisualizer::setupMeshBuffers()
{
    // Create a Vertex Array Object for the surface mesh
    glGenVertexArrays(1, &fluidVAO);
    glBindVertexArray(fluidVAO);
    // Create a Vertex Buffer Object for the surface mesh
    glGenBuffers(1, &fluidVBO);
    glBindBuffer(GL_ARRAY_BUFFER, fluidVBO);
    // Create an Element Buffer Object for the surface mesh
    glGenBuffers(1, &fluidEBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, fluidEBO);

    // Specify the layout of the surface mesh vertex data
    // PolyVox::PositionMaterialNormal layout:
    //    (x, y, z)-position (3 floats)
    //    (x, y, z)-normal   (3 floats)
    //    material           (1 float)
    GLint phongShaderPosAttrib = glGetAttribLocation(phongShaderProgram, "position");
    glEnableVertexAttribArray(phongShaderPosAttrib);
    glVertexAttribPointer(phongShaderPosAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(0 * sizeof(float)));
    GLint phongShaderNormAttrib = glGetAttribLocation(phongShaderProgram, "normal");
    glEnableVertexAttribArray(phongShaderNormAttrib);
    glVertexAttribPointer(phongShaderNormAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(3 * sizeof(float)));
    GLint phongShaderMatAttrib = glGetAttribLocation(phongShaderProgram, "material");
    glEnableVertexAttribArray(phongShaderMatAttrib);
    glVertexAttribPointer(phongShaderMatAttrib, 1, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(6 * sizeof(float)));

    // Create Vertex Array Objects for the models
    modelVAOs.resize(models.size());
    glGenVertexArrays((GLsizei)models.size(), modelVAOs.data());
    // Create Vertex Buffer Objects for the models
    modelVBOs.resize(models.size());
    glGenBuffers((GLsizei)models.size(), modelVBOs.data());
    for (size_t i = 0; i != models.size(); ++i)
    {
        glBindVertexArray(modelVAOs[i]);
        glBindBuffer(GL_ARRAY_BUFFER, modelVBOs[i]);
        // Upload vertex data
        glBufferData(GL_ARRAY_BUFFER, models[i].vertexData.size() * sizeof(PolyVox::PositionMaterialNormal), models[i].vertexData.data(), GL_STATIC_DRAW);

        // Specify the layout of the model vertex data
        // PolyVox::PositionMaterialNormal layout:
        //    (x, y, z)-position (3 floats)
        //    (x, y, z)-normal   (3 floats)
        //    material           (1 float)
        glEnableVertexAttribArray(phongShaderPosAttrib);
        glVertexAttribPointer(phongShaderPosAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(0 * sizeof(float)));
        glEnableVertexAttribArray(phongShaderNormAttrib);
        glVertexAttribPointer(phongShaderNormAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(3 * sizeof(float)));
        glEnableVertexAttribArray(phongShaderMatAttrib);
        glVertexAttribPointer(phongShaderMatAttrib, 1, GL_FLOAT, GL_FALSE, sizeof(PolyVox::PositionMaterialNormal), reinterpret_cast<void*>(6 * sizeof(float)));
    }
}

FluidVisualizer::~FluidVisualizer()
{
    glDeleteProgram(simpleShaderProgram);
    glDeleteShader(simpleFragmentShader);
    glDeleteShader(simpleVertexShader);
    glDeleteProgram(phongShaderProgram);
    glDeleteShader(phongFragmentShader);
    glDeleteShader(phongVertexShader);

    glDeleteBuffers(1, &fluidEBO);
    glDeleteBuffers(1, &fluidVBO);
    glDeleteVertexArrays(1, &fluidVAO);
    glDeleteBuffers(1, &pointsVBO);
    glDeleteVertexArrays(1, &pointsVAO);
    glDeleteBuffers(1, &boundsVBO);
    glDeleteVertexArrays(1, &boundsVAO);

    glDeleteBuffers((GLsizei)modelVBOs.size(), modelVBOs.data());
    glDeleteVertexArrays((GLsizei)modelVAOs.size(), modelVAOs.data());
}

// Draws a new frame.
void FluidVisualizer::draw()
{

    const glm::mat4 view = camera.getViewMatrix();
    glUseProgram(phongShaderProgram);
    glUniformMatrix4fv(phongViewUniform, 1, GL_FALSE, glm::value_ptr(view));
    glUniform3fv(phongCamUniform, 1, glm::value_ptr(camera.getPosition()));

    // Draw models
    glm::mat4 modelMatrix = glm::translate(sceneOffset);
    glm::mat3 normalMatrix = glm::transpose(glm::inverse(modelMatrix));
    glUniformMatrix4fv(phongModelUniform, 1, GL_FALSE, glm::value_ptr(modelMatrix));
    glUniformMatrix3fv(phongNormalMatUniform, 1, GL_FALSE, glm::value_ptr(normalMatrix));
    for (size_t i = 0; i != models.size(); ++i)
    {
        glBindVertexArray(modelVAOs[i]);
        glUniform3fv(phongAmbientUniform, 1, glm::value_ptr(models[i].ambientColor));
        glUniform3fv(phongDiffuseUniform, 1, glm::value_ptr(models[i].diffuseColor));
        glUniform3fv(phongSpecularUniform, 1, glm::value_ptr(models[i].specularColor));
        glUniform1f(phongShininessUniform, models[i].specularExponent);
        glUniform1f(phongAlphaUniform, 1.0f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)models[i].vertexData.size());
    }

    glUniform3fv(phongAmbientUniform, 1, glm::value_ptr(glm::vec3{0.5f, 0.5f, 0.5f}));
    glUniform3fv(phongDiffuseUniform, 1, glm::value_ptr(glm::vec3{0.5f, 0.5f, 0.95f}));
    glUniform3fv(phongSpecularUniform, 1, glm::value_ptr(glm::vec3{1.0f, 1.0f, 1.0f}));
    glUniform1f(phongShininessUniform, 32.0f);

    // Draw surface mesh
    if (meshConstruction)
    {
        // Run marching cubes using PolyVox, and retrieve the vertex and index buffers
        fillVoxelVolume();
        surfaceExtractor.execute();
        const auto& lowerCorner = voxelVolume.getEnclosingRegion().getLowerCorner();
        surfaceMesh.translateVertices({ (float)lowerCorner.getX(), (float)lowerCorner.getY(), (float)lowerCorner.getZ() });
        surfaceMesh.scaleVertices(1.0f / voxelVolumeResolutionScale);
        const std::vector<uint32_t>& waterMeshIndices = surfaceMesh.getIndices();
        std::vector<PolyVox::PositionMaterialNormal>& waterMeshVertices = surfaceMesh.getRawVertexData();
        for (auto& vert : waterMeshVertices) // Clamp vertex locations to world boundaries
        {
            vert.position.setElements(
                clamp(minPos.x + sceneOffset.x, maxPos.x + sceneOffset.x, vert.position.getX() + sceneOffset.x),
                clamp(minPos.y + sceneOffset.y, maxPos.y + sceneOffset.y, vert.position.getY() + sceneOffset.y),
                clamp(minPos.z + sceneOffset.z, maxPos.z + sceneOffset.z, vert.position.getZ() + sceneOffset.z)
            );
        }

        glUniformMatrix4fv(phongModelUniform, 1, GL_FALSE, glm::value_ptr(glm::mat4{}));  // Identity matrix
        glUniformMatrix4fv(phongNormalMatUniform, 1, GL_FALSE, glm::value_ptr(glm::mat4{}));  // Identity matrix
        glUniform1f(phongAlphaUniform, 0.75f);
        glBindVertexArray(fluidVAO);
        glBindBuffer(GL_ARRAY_BUFFER, fluidVBO);
        glBufferData(GL_ARRAY_BUFFER, waterMeshVertices.size() * sizeof(PolyVox::PositionMaterialNormal), waterMeshVertices.data(), GL_STREAM_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, waterMeshIndices.size() * sizeof(uint32_t), waterMeshIndices.data(), GL_STREAM_DRAW);
        glDrawElements(GL_TRIANGLES, (GLsizei)waterMeshIndices.size(), GL_UNSIGNED_INT, nullptr);
    }
    

    glUseProgram(simpleShaderProgram);
    glUniformMatrix4fv(simpleViewUniform, 1, GL_FALSE, glm::value_ptr(view));

    // Draw points
    glUniformMatrix4fv(simpleModelUniform, 1, GL_FALSE, glm::value_ptr(glm::mat4{}));  // Identity matrix
    glBindVertexArray(pointsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(r), r, GL_STREAM_DRAW);
    glDrawArrays(GL_POINTS, 0, (GLsizei)std::size(r));

    // Draw world bounds
    glUniformMatrix4fv(simpleModelUniform, 1, GL_FALSE, glm::value_ptr(glm::translate(sceneOffset)));
    glBindVertexArray(boundsVAO);
    glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)std::size(worldBoundsVertices));
}

// Returns the voxel index of a world position
PolyVox::Vector3DInt32 FluidVisualizer::worldPosToVoxelIndex(const glm::vec3& worldPos) const
{
    return {
        (int)((worldPos.x - sceneOffset.x) * voxelVolumeResolutionScale),
        (int)((worldPos.y - sceneOffset.y) * voxelVolumeResolutionScale),
        (int)((worldPos.z - sceneOffset.z) * voxelVolumeResolutionScale)
    };
}

// Returns the world position of a voxel index
glm::vec3 FluidVisualizer::voxelIndexToWorldPos(int voxelX, int voxelY, int voxelZ) const
{
    static const float invVoxelVolumeResolutionScale = 1.0f / voxelVolumeResolutionScale;

    return {
        voxelX * invVoxelVolumeResolutionScale + sceneOffset.x,
        voxelY * invVoxelVolumeResolutionScale + sceneOffset.y,
        voxelZ * invVoxelVolumeResolutionScale + sceneOffset.z,
    };
}

// Fills the voxel volume: for each voxel, a density is calculated. If the density is above a certain threshold, the voxel will be considered solid
// For the threshold, see PolyVox::DefaultMarchingCubesController<float>::getThreshold() at the start of FluidVisualizer.h
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

// Load a scene with a given scene number. The fluid is not reset when this is done.
void FluidVisualizer::loadScene(int sceneNumber)
{
    FluidBase::loadScene(sceneNumber);
    setupMeshBuffers();

    // fillTriangleGrid() assumes that the sceneOffset is 0
    glm::vec3 sceneOffsetBackup = sceneOffset;
    sceneOffset = glm::vec3();
    fillTriangleGrid();
    sceneOffset = sceneOffsetBackup;
}
