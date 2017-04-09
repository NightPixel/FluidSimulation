#pragma once

#include "FluidSimulator.h"
#include <PolyVoxCore/MarchingCubesSurfaceExtractor.h>
#include <PolyVoxCore/SurfaceMesh.h>
#include <PolyVoxCore/SimpleVolume.h>

// Returns the density value above which a voxel is considered solid.
// Template specialization of PolyVox::DefaultMarchingCubesController<float>::getThreshold()
template<> inline PolyVox::DefaultMarchingCubesController<float>::DensityType
PolyVox::DefaultMarchingCubesController<float>::getThreshold()
{
    return 5.0f;
}

class FluidVisualizer : public FluidSimulator
{
public:
    explicit FluidVisualizer(GLFWwindow* window);
    ~FluidVisualizer();

    // Whether or not to execute mesh construction using marching cubes
    bool meshConstruction = true;

    // Draws a new frame.
    void draw();

private:
    GLuint fluidVAO, pointsVAO, boundsVAO;
    GLuint fluidVBO, pointsVBO, boundsVBO;
    GLuint fluidEBO;
    GLuint simpleVertexShader, simpleFragmentShader, simpleShaderProgram;
    GLuint phongVertexShader, phongFragmentShader, phongShaderProgram;
    GLint simpleModelUniform, simpleViewUniform;
    GLint phongModelUniform, phongViewUniform, phongNormalMatUniform, phongCamUniform, phongAmbientUniform,
        phongDiffuseUniform, phongSpecularUniform, phongShininessUniform, phongAlphaUniform;
    std::vector<GLuint> modelVAOs;
    std::vector<GLuint> modelVBOs;

    PolyVox::Vector3DInt32 worldPosToVoxelIndex(const glm::vec3& worldPos) const;
    glm::vec3 voxelIndexToWorldPos(int voxelX, int voxelY, int voxelZ) const;
    void fillVoxelVolume();

    // World positions will be multiplied by this scale for the purposes of voxel indexing.
    // Example: if this scale is 10, a world pos of (-1, 0, 2.5) will map to the voxel at (-10, 0, 25).
    // Higher values result in a voxel grid of a higher resolution; individual voxels would be smaller.
    float voxelVolumeResolutionScale = 10.0f;
    PolyVox::SimpleVolume<float> voxelVolume;
    PolyVox::SurfaceMesh<PolyVox::PositionMaterialNormal> surfaceMesh;
    PolyVox::MarchingCubesSurfaceExtractor<PolyVox::SimpleVolume<float>> surfaceExtractor;
};
