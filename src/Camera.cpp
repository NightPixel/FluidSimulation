#include "Camera.h"
#include <glm/gtc/matrix_transform.hpp>

// Rotates the camera around the target position it's focusing on.
void Camera::rotate(float dTheta, float dPhi)
{
    static const float pi = 3.14159265358979323846f;
    static const float twoPi = 2.0f * pi;

    if (up > 0.0f)
        theta += dTheta;
    else
        theta -= dTheta;

    phi += dPhi;

    // Phi should be in the range [-2pi, 2pi]
    if (phi > twoPi)
        phi -= twoPi;
    else if (phi < -twoPi)
        phi += twoPi;

    // up is positive Y if phi is in [0, pi] or [-2pi, -pi].
    // Otherwise, it is negative Y.
    if (0 < phi && phi < pi || -twoPi < phi && phi < -pi)
        up = 1.0f;
    else
        up = -1.0f;
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
    glm::vec3 worldUp(0.0f, up, 0.0f);
    glm::vec3 right = glm::cross(look, worldUp);
    glm::vec3 up = glm::cross(look, right);
    target += right * dx + up * dy;
}

// Gets the camera's view matrix.
glm::mat4 Camera::getViewMatrix() const
{
    return glm::lookAt(getPosition(), target, glm::vec3(0.0f, up, 0.0f));
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
