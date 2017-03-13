#version 150

in vec3 Color;

out vec4 outColor;

void main()
{
    vec3 pixelColor = Color + vec3(0.9, 0.9, 0.9);
    pixelColor /= 2;
    outColor = vec4(pixelColor, 1.0);
}