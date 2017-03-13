#include "Camera.h"
#include <glm/gtc/matrix_transform.hpp>

// Rotates the camera around the target position it's focusing on.
void Camera::rotate(float dTheta, float dPhi)
{
    theta += dTheta;
    phi += dPhi;
}

// Zooms the camera towards the target.
void Camera::zoom(float distance)
{
    radius -= distance * radius;
}

// Pans the camera, moving the target position parallel to the camera plane.
void Camera::pan(float dx, float dy)
{
    glm::vec3 look = glm::normalize(toCartesian());
    glm::vec3 worldUp(0.0f, 1.0f, 0.0f);
    glm::vec3 right = glm::cross(look, worldUp);
    glm::vec3 up = glm::cross(look, right);
    target += right * dx + up * dy;
}

// Gets the camera's view matrix.
glm::mat4 Camera::getViewMatrix() const
{
    return glm::lookAt(getPosition(), target, glm::vec3(0.0f, 1.0f, 0.0f));
}

// Gets the camera position in spherical coordinates; the target position is the origin.
glm::vec3 Camera::toCartesian() const
{
    float x = radius * sin(phi) * sin(theta);
    float y = radius * cos(phi);
    float z = radius * sin(phi) * cos(theta);
        
    return glm::vec3(x, y, z);
}

// Gets the camera position in world coordinates.
glm::vec3 Camera::getPosition() const
{
    return target + toCartesian();
}

