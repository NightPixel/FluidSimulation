#include "FluidBase.h"
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include <tiny_obj_loader.h>

FluidBase::FluidBase(GLFWwindow* window)
    : window(window)
{
    glfwGetWindowSize(window, &windowSizeX, &windowSizeY);
    antTweakBar = TwNewBar("Simulation settings");
    TwDefine("GLOBAL fontsize=3");


    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, "data/cube.obj", "data/");

    for (size_t v = 0; v < attrib.vertices.size(); v += 3)
        objectPositions.emplace_back(attrib.vertices[v], attrib.vertices[v + 1], attrib.vertices[v + 2]);

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++)
    {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
        {
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
                objectVertices.emplace_back(PolyVox::Vector3DFloat{ vx, vy, vz }, PolyVox::Vector3DFloat{ nx, ny, nz }, shapes[s].mesh.material_ids[f]);
                //float tx = attrib.texcoords[2 * idx.texcoord_index + 0];
                //float ty = attrib.texcoords[2 * idx.texcoord_index + 1];
            }
            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];
        }
    }
}

void FluidBase::onMouseMoved(float dxPos, float dyPos)
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

void FluidBase::onMouseScrolled(float yOffset)
{
    camera.zoom(yOffset / 30.0f);
}

void FluidBase::onKeypress(int key, int action)
{
    printf("K: %i A: %i\n", key, action);
    switch (key)
    {
    case 'W':
        holdForward = action != GLFW_RELEASE;
        break;
    case 'S':
        holdBackward = action != GLFW_RELEASE;
        break;
    case 'A':
        holdLeft = action != GLFW_RELEASE;
        break;
    case 'D':
        holdRight = action != GLFW_RELEASE;
        break;
    case 32: // Space
        holdUp = action != GLFW_RELEASE;
        break;
    case 341: // Left control
        holdDown = action != GLFW_RELEASE;
        break;
    case 340: // Left shift
        holdShift = action != GLFW_RELEASE;
        break;
    }
}
