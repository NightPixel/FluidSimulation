#version 150

in vec3 position;
in vec3 normal;
in float material;

out vec3 WorldPos;
out vec3 Normal;
out float Material;

uniform mat4 model;
uniform mat4 view;
uniform mat4 proj;

void main()
{
    WorldPos = vec3(model * vec4(position, 1.0));
    Normal = normal;
    Material = material;
    gl_Position = proj * view * model * vec4(position, 1.0);
}