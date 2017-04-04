#include "Geometry.h"
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include <tiny_obj_loader.h>
#include <glm/gtx/rotate_vector.hpp>

std::vector<Model> loadOBJFile(const std::string& fileName, const glm::vec3& offset, const glm::vec3& rotation, const glm::vec3& scale)
{
    std::vector<Model> models;
    glm::vec3 radiansRotation = glm::radians(rotation);

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

                positions[v] = glm::vec3{ vx, vy, vz };
                normals[v] = glm::vec3{ nx, ny, nz };

                positions[v] *= scale;
                positions[v] = glm::rotateX(positions[v], radiansRotation.x);
                positions[v] = glm::rotateY(positions[v], radiansRotation.y);
                positions[v] = glm::rotateZ(positions[v], radiansRotation.z);
                positions[v] += offset;
                normals[v] = glm::rotateX(normals[v], radiansRotation.x);
                normals[v] = glm::rotateY(normals[v], radiansRotation.y);
                normals[v] = glm::rotateZ(normals[v], radiansRotation.z);

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
        glm::vec3 ambientColor{ materials[0].ambient[0], materials[0].ambient[1], materials[0].ambient[2] };
        glm::vec3 diffuseColor{ materials[0].diffuse[0], materials[0].diffuse[1], materials[0].diffuse[2] };
        glm::vec3 specularColor{ materials[0].specular[0], materials[0].specular[1], materials[0].specular[2] };
        float specularExponent = materials[0].shininess;
        models.emplace_back(triangles, ambientColor, diffuseColor, specularColor, specularExponent, vertexData);
    }

    return models;
}
