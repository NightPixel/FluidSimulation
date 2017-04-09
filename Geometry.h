#pragma once

#include <glm/glm.hpp>
#include <PolyVoxCore/VertexTypes.h>
#include <utility>
#include <algorithm>
#include <vector>

struct Triangle
{
    glm::vec3 positions[3];
    glm::vec3 normal;

    Triangle(const glm::vec3& pos0, const glm::vec3& pos1, const glm::vec3& pos2, const glm::vec3& norm)
        : normal(norm)
    {
        positions[0] = pos0;
        positions[1] = pos1;
        positions[2] = pos2;
    }

    std::pair<glm::vec3, glm::vec3> getBoundingBox() const
    {
        return std::make_pair(
            glm::vec3{
                std::min({positions[0].x, positions[1].x, positions[2].x}),
                std::min({positions[0].y, positions[1].y, positions[2].y}),
                std::min({positions[0].z, positions[1].z, positions[2].z}),
            },
            glm::vec3{
                std::max({positions[0].x, positions[1].x, positions[2].x}),
                std::max({positions[0].y, positions[1].y, positions[2].y}),
                std::max({positions[0].z, positions[1].z, positions[2].z}),
            }
        );
    }
};

struct Model
{
    std::vector<Triangle> triangles;
    glm::vec3 ambientColor;
    glm::vec3 diffuseColor;
    glm::vec3 specularColor;
    float specularExponent;
    std::vector<PolyVox::PositionMaterialNormal> vertexData;

    Model(
        const std::vector<Triangle>& tris,
        const glm::vec3& ambient,
        const glm::vec3& diffuse,
        const glm::vec3& specular,
        float exponent,
        const std::vector<PolyVox::PositionMaterialNormal> vertices
    ) :
        triangles(tris),
        ambientColor(ambient),
        diffuseColor(diffuse),
        specularColor(specular),
        specularExponent(exponent),
        vertexData(vertices)
    { }
};

std::vector<Model> loadOBJFile(const std::string& fileName,
    const glm::vec3& offset = {}, const glm::vec3& rotation = {}, const glm::vec3& scale = glm::vec3{1.0f, 1.0f, 1.0f});

bool triangleBoxIntersection(const Triangle& triangle, const glm::vec3& boxCenter, const glm::vec3& boxHalfSize);

bool planeBoxIntersection(const glm::vec3& normal, const glm::vec3& vertex, const glm::vec3& boxHalfSize);

std::pair<bool, float> triangleLineSegmentIntersection(const Triangle& triangle, const glm::vec3& segmentStart, const glm::vec3& segmentEnd, const glm::vec3& sceneOffset);
