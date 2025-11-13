#version 330

// Input vertex attributes (from vertex shader)
in vec2 fragTexCoord;
in vec4 fragColor;

uniform float radius;
uniform sampler2D texture0;

// Output fragment color
out vec4 finalColor;

vec3 simple_tm(vec3 hdrColor) {
    return hdrColor / (hdrColor + vec3(1.0));
}

// adapted from https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
vec3 aces_tm(vec3 x) {
    vec3 a = vec3(2.51f);
    vec3 b = vec3(0.03f);
    vec3 c = vec3(2.43f);
    vec3 d = vec3(0.59f);
    vec3 e = vec3(0.14f);
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

vec3 exposure_tm(vec3 hdrColor) {
    vec3 exposure = vec3(0.5);
    return vec3(1.0) - exp(-hdrColor * exposure);
}

void main() {
    float gamma = 2.2;
    vec3 hdrColor = texture(texture0, fragTexCoord).rgb;
    vec3 mapped = simple_tm(hdrColor);
    mapped = pow(mapped, vec3(1.0 / gamma));
    finalColor = vec4(mapped, fragColor.a);
}
