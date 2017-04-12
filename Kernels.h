#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <algorithm>

static const float pi = 3.14159265358979323846f;

inline float poly6(const glm::vec3& rVec, float h)
{
    return 315 / (64 * pi * pow(h, 9)) * pow(std::max(0.0f, h*h - glm::length2(rVec)), 3);
}

// Overload of poly6 that takes the squared length of r directly, allowing for the glm::length2 call to be 
// skipped when the squared length of r is already known
inline float poly6(const float rSqLen, float h)
{
    return 315 / (64 * pi * pow(h, 9)) * pow(std::max(0.0f, h*h - rSqLen), 3);
}

inline glm::vec3 poly6Gradient(const glm::vec3& rVec, float h)
{
    return -945 / (32 * pi * pow(h, 9)) * rVec * pow(std::max(0.0f, h*h - glm::length2(rVec)), 2);
}

inline float poly6Laplacian(const glm::vec3& rVec, float h)
{
    const float rNorm2 = glm::length2(rVec);
    return -945 / (32 * pi * pow(h, 9)) * std::max(0.0f, h*h - rNorm2) * (3*h*h - 7*rNorm2);
}

inline glm::vec3 spikyGradient(const glm::vec3& rVec, float h)
{
    return -45 / (pi * pow(h, 6)) * glm::normalize(rVec) * pow(std::max(0.0f, h - glm::length(rVec)), 2);
}

inline float viscosityLaplacian(const glm::vec3& rVec, float h)
{
    return 45 / (pi * pow(h, 6)) * std::max(0.0f, h - glm::length(rVec));
}