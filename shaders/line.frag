#version 450 core

in vec3 vColor;

uniform int uRound;  // 1 : les points sont dessinés en disques

out vec4 fragColor;

void main() {
    if (uRound == 1) {
        vec2 c = gl_PointCoord * 2.0 - 1.0;
        if (dot(c, c) > 1.0) discard;
    }
    fragColor = vec4(vColor, 1.0);
}
