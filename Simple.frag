#version 150

in vec3 Color;

out vec4 outColor;

void main()
{
    vec3 pixelColor = Color + vec3(2.0, 2.0, 2.0);
    pixelColor /= 4;
    outColor = vec4(pixelColor, 1.0);
}