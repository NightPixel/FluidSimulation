#pragma once

#include "Camera.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <AntTweakBar.h>
#include <PolyVoxCore/VertexTypes.h>
#include <vector>
#include <algorithm>

class FluidBase
{
public:
    explicit FluidBase(GLFWwindow* window);

    // Called when the mouse cursor is moved.
    void onMouseMoved(float dxPos, float dyPos);
    // Called when the mouse wheel is scrolled.
    void onMouseScrolled(float yOffset);
    // Called when a key is pressed or released.
    void onKeypress(int key, int action);

protected:
    Camera camera;
    GLFWwindow* window;
    TwBar* antTweakBar;
    int windowSizeX;
    int windowSizeY;

    bool holdForward = false;
    bool holdBackward = false;
    bool holdRight = false;
    bool holdLeft = false;
    bool holdUp = false;
    bool holdDown = false;
    bool holdShift = false;

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
    std::vector<Triangle> objectTriangles;
    std::vector<PolyVox::PositionMaterialNormal> objectVertices;
};
