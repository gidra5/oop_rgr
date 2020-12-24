#version 150 core

out vec3 normal;

in ivec3 a_pos;
in ivec3 a_normal;

uniform mat4 u_proj;
uniform mat4 u_view;

void main() {
    gl_Position = u_proj * u_view * vec4(a_pos, 1.0);
    normal = a_normal;
}