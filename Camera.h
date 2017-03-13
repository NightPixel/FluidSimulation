#pragma once

#include <glm/glm.hpp>

// A trackball camera, based on https://computergraphics.stackexchange.com/a/168
class Camera
{
public:
    // Rotates the camera around the target position it's focusing on.
    void rotate(float dTheta, float dPhi);

    // Zooms the camera towards the target.
    void zoom(float distance);

    // Pans the camera, moving the target position parallel to the camera plane.
    void pan(float dx, float dy);

    // Gets the camera's view matrix.
    glm::mat4 getViewMatrix() const;

private:
    // Gets the camera position in spherical coordinates; the target position is the origin.
    glm::vec3 toCartesian() const;
    
    // Gets the camera position in world coordinates.
    glm::vec3 getPosition() const;

    // The camera position is stored in spherical coordinates.
    // Note that a different convention is used from the source implementation on StackOverflow;
    // the y-axis is the upwards axis in this program, not the z-axis.

    // Theta is the angle on the XZ plane; 0 theta is on the x-axis, pi/2 theta is on the z-axis.
    float theta = glm::radians(45.0f);
    // Phi is the angle away from the y-axis; 0 phi is on the y-axis, pi/2 phi is on the XZ plane.
    float phi = glm::radians(45.0f);
    // Radius is the radius of the sphere around the target on which the camera is located.
    float radius = 5.0f;

    // Target is a position in world coordinates; the camera looks at this position.
    glm::vec3 target;
};
