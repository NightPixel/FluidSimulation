#version 150

in vec3 WorldPos;
in vec3 Normal;
in float Material;

out vec4 outColor;

uniform vec3 camPos;
uniform vec3 lightPos;

uniform vec3 ambientMaterialColor;
uniform vec3 ambientLightColor;

uniform vec3 diffuseMaterialColor;
uniform vec3 diffuseLightColor;

uniform vec3 specularMaterialColor;
uniform vec3 specularLightColor;
uniform float shininess;
uniform float alpha;

void main()
{
    vec3 ambient = ambientMaterialColor * ambientLightColor;

    vec3 N = normalize(Normal);
    vec3 L = normalize(lightPos - WorldPos);
    float diffuseFactor = max(dot(N, L), 0.0);
    vec3 diffuse = diffuseMaterialColor * diffuseFactor * diffuseLightColor;

    vec3 V = normalize(camPos - WorldPos);
    vec3 R = reflect(-V, N);
    float specularFactor = pow(max(dot(R, V), 0.0), shininess) * ceil(diffuseFactor);
    vec3 specular = specularMaterialColor * specularFactor * specularLightColor;

    outColor = vec4(ambient + diffuse + specular, alpha);
}