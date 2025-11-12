#version 330

// Input vertex attributes (from vertex shader)
in vec2 fragTexCoord;
in vec4 fragColor;

uniform float radius;
uniform sampler2D texture0;

// Output fragment color
out vec4 finalColor;

void main() {
    float gamma = 2.2;
    vec3 hdrColor = texture(texture0, fragTexCoord).rgb;
    vec3 mapped = hdrColor / (hdrColor + vec3(1.0));
    mapped = pow(mapped, vec3(1.0 / gamma));
    finalColor = vec4(mapped, fragColor.a);
}
