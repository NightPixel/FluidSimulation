#version 150

in vec3 Color;
in vec3 Normal;
in float Material;

out vec4 outColor;

uniform vec3 minWorldPos;
uniform vec3 maxWorldPos;

float invLerp(float lower, float upper, float value)
{
    return (value - lower) / (upper - lower);
}

vec3 invLerp(vec3 lower, vec3 upper, vec3 value)
{
    return vec3(
        invLerp(lower.x, upper.x, value.x),
        invLerp(lower.y, upper.y, value.y),
        invLerp(lower.z, upper.z, value.z)
    );
}

void main()
{
    outColor = vec4(invLerp(minWorldPos, maxWorldPos, Color), 1.0);
}