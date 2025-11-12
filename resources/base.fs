#version 330

// Input vertex attributes (from vertex shader)
in vec2 fragTexCoord;
in vec4 fragColor;

uniform float radius;

// Output fragment color
out vec4 finalColor;

void main() {
    vec2 p = (fragTexCoord - vec2(0.5)) * 2;
    float l = sqrt(p.x*p.x + p.y*p.y);
    float alpha = l < radius ? 1 : 0;
    float glow = 0.6*exp((-l/radius)*1.0);
    // float intensity = smoothstep(0.05, 0, l);
    // glow += intensity;
    // float glow = 2-exp(l);
    // float glow = 1-sqrt(l);
    // finalColor = vec4(fragColor.rgb, fragColor.a*glow);
    alpha += glow;
    finalColor = vec4(fragColor.rgb, fragColor.a*alpha);
}
