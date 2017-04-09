#include "Geometry.h"
#include "Definitions.h"
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include <tiny_obj_loader.h>
#include <glm/gtx/transform.hpp>

std::vector<Model> loadOBJFile(const std::string& fileName, const glm::vec3& offset, const glm::vec3& rotation, const glm::vec3& scale)
{
    std::vector<Model> models;
    glm::vec3 radiansRotation = glm::radians(rotation);
    glm::mat4 T = glm::translate(offset);
    glm::mat4 R =
        glm::rotate(radiansRotation.z, glm::vec3{0.0f, 0.0f, 1.0f}) *
        glm::rotate(radiansRotation.y, glm::vec3{0.0f, 1.0f, 0.0f}) *
        glm::rotate(radiansRotation.x, glm::vec3{1.0f, 0.0f, 0.0f});
    glm::mat4 S = glm::scale(scale);
    glm::mat4 modelMatrix = T * R * S;
    glm::mat3 normalMatrix = glm::transpose(glm::inverse(modelMatrix));

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    tinyobj::LoadObj(&attrib, &shapes, &materials, &err, ("data/" + fileName).c_str(), "data/", true);

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++)
    {
        std::vector<Triangle> triangles;
        std::vector<PolyVox::PositionMaterialNormal> vertexData;

        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
        {
            // tinyobj::LoadObj()s 'triangulate' parameter was set to true,
            // so we should only get faces with three vertices
            assert(shapes[s].mesh.num_face_vertices[f] == 3);
            glm::vec3 positions[3];
            glm::vec3 normals[3];

            int fv = shapes[s].mesh.num_face_vertices[f];
            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++)
            {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                float vx = attrib.vertices[3 * idx.vertex_index + 0];
                float vy = attrib.vertices[3 * idx.vertex_index + 1];
                float vz = attrib.vertices[3 * idx.vertex_index + 2];
                float nx = attrib.normals[3 * idx.normal_index + 0];
                float ny = attrib.normals[3 * idx.normal_index + 1];
                float nz = attrib.normals[3 * idx.normal_index + 2];
                printf("Vertex [%i](%f, %f, %f), Normal [%i](%f, %f, %f)\n", idx.vertex_index, vx, vy, vz, idx.normal_index, nx, ny, nz);

                positions[v] = modelMatrix * glm::vec4{vx, vy, vz, 1.0f};
                normals[v] = glm::normalize(normalMatrix * glm::vec3{nx, ny, nz});

                vertexData.emplace_back(
                    PolyVox::Vector3DFloat{positions[v].x, positions[v].y, positions[v].z},
                    PolyVox::Vector3DFloat{normals[v].x, normals[v].y, normals[v].z},
                    (float)shapes[s].mesh.material_ids[f]);
            }
            triangles.emplace_back(positions[0], positions[1], positions[2], glm::normalize(normals[0] + normals[1] + normals[2]));

            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];
        }

        // We only support one material per shape (not one material per triangle).
        // We simply pick the first material, and use its properties.
        glm::vec3 ambientColor{materials[0].ambient[0], materials[0].ambient[1], materials[0].ambient[2]};
        glm::vec3 diffuseColor{materials[0].diffuse[0], materials[0].diffuse[1], materials[0].diffuse[2]};
        glm::vec3 specularColor{materials[0].specular[0], materials[0].specular[1], materials[0].specular[2]};
        float specularExponent = materials[0].shininess;
        models.emplace_back(triangles, ambientColor, diffuseColor, specularColor, specularExponent, vertexData);
    }

    return models;
}

bool triangleBoxIntersection(const Triangle& triangle, const glm::vec3& boxCenter, const glm::vec3& boxHalfSize)
{
    // First, we create a new, translated triangle such that the
    // box we're testing it against is centered around the origin.
    const Triangle translatedTriangle{
        triangle.positions[0] - boxCenter,
        triangle.positions[1] - boxCenter,
        triangle.positions[2] - boxCenter,
        triangle.normal
    };

    // Each triangle's bounding box is described by a (lower corner, upper corner) vector pair.
    glm::vec3 minCoords, maxCoords;
    std::tie(minCoords, maxCoords) = translatedTriangle.getBoundingBox();

    // 1. Three tests: triangle AABB against box
    if (minCoords.x > boxHalfSize.x || maxCoords.x < -boxHalfSize.x ||
        minCoords.y > boxHalfSize.y || maxCoords.y < -boxHalfSize.y ||
        minCoords.z > boxHalfSize.z || maxCoords.z < -boxHalfSize.z)
        return false;

    // 2. One test: triangle plane against box
    if (!planeBoxIntersection(translatedTriangle.normal, translatedTriangle.positions[0], boxHalfSize))
        return false;

    // 3. Nine tests: each combination of some global axis crossed with a triangle edge
    const glm::vec3 axes[3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f}
    };
    const glm::vec3 edges[3] = {
        translatedTriangle.positions[1] - translatedTriangle.positions[0],
        translatedTriangle.positions[2] - translatedTriangle.positions[1],
        translatedTriangle.positions[0] - translatedTriangle.positions[2]
    };

    for (size_t axis = 0; axis != 3; ++axis)
        for (size_t edge = 0; edge != 3; ++edge)
        {
            const glm::vec3 a = glm::cross(axes[axis], edges[edge]);
            // Project the triangle vertices onto a
            const float p0 = glm::dot(a, translatedTriangle.positions[0]);
            const float p1 = glm::dot(a, translatedTriangle.positions[1]);
            const float p2 = glm::dot(a, translatedTriangle.positions[2]);
            // Calculate the "radius" of the box projected on a
            const float r = boxHalfSize.x * abs(a.x) + boxHalfSize.y * abs(a.y) + boxHalfSize.z * abs(a.z);
            if (std::min({p0, p1, p2}) > r || std::max({p0, p1, p2}) < -r)
                return false;
        }

    return true;
}

bool planeBoxIntersection(const glm::vec3& normal, const glm::vec3& vertex, const glm::vec3& maxBox)
{
    glm::vec3 vMin, vMax;
    for (int q = 0; q != 2; ++q)
    {
        const float v = vertex[q];
        if (normal[q] > 0.0f)
        {
            vMin[q] = -maxBox[q] - v;
            vMax[q] =  maxBox[q] - v;
        }
        else
        {
            vMin[q] =  maxBox[q] - v;
            vMax[q] = -maxBox[q] - v;
        }
    }
    if (glm::dot(normal, vMin) > 0.0f)
        return false;
    if (glm::dot(normal, vMax) >= 0.0f)
        return true;

    return false;
}

std::pair<bool, float> triangleLineSegmentIntersection(const Triangle& triangle, const glm::vec3& segmentStart, const glm::vec3& segmentEnd)
{
    // First test whether the line segment intersects with the plane.
    // If it does, we test whether the intersection point also lies inside the triangle.
    // Only if that is true we report an intersection.
    const float planeDenominator = glm::dot(triangle.normal, segmentEnd - segmentStart);
    if (std::abs(planeDenominator) < EPSILON)
        return {false, {}}; // The segment is parallel to the plane.

    const float planeNumerator = glm::dot(triangle.normal, triangle.positions[0] - segmentStart);

    // The line segment intersects with the plane at
    // segmentStart + r * (segmentEnd - segmentStart)
    const float r = planeNumerator / planeDenominator;
    if (r < 0.0f)
        return {false, {}}; // The segment points away from the triangle
    if (r > 1.0f)
        return {false, {}}; // The segment ends before hitting the plane
    // 0 <= r <= 1: the segment hits the plane. Does the intersection point also lie inside the triangle?
    const glm::vec3 intersection = segmentStart + r * (segmentEnd - segmentStart);

    // Points in the plane can be described by the equation
    // v(s, t) = v0 + su + tv
    // where u is one triangle edge pointing away from v0, and v is the other edge.
    const glm::vec3 u = triangle.positions[1] - triangle.positions[0];
    const glm::vec3 v = triangle.positions[2] - triangle.positions[0];
    const glm::vec3 w = intersection - triangle.positions[0]; // w is a vector in the plane

    // Pre-calculate some dot products
    const float uu = glm::dot(u, u);
    const float uv = glm::dot(u, v);
    const float vv = glm::dot(v, v);
    const float wu = glm::dot(w, u);
    const float wv = glm::dot(w, v);
    const float triangleDenominator = uv * uv - uu * vv;

    // The intersection point lies inside the triangle if
    // s >= 0, t >= 0, and s + t <= 1
    const float s = (uv * wv - vv * wu) / triangleDenominator;
    if (s < 0.0f || s > 1.0f)
        return {false, {}};
    const float t = (uv * wu - uu * wv) / triangleDenominator;
    if (t < 0.0f || s + t > 1.0f)
        return {false, {}};

    return {true, r};
}
