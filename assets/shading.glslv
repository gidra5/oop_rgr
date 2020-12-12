#version 150 core

out vec3 frag_pos;
out vec3 normal;

in ivec3 a_pos;
in ivec3 a_normal;

uniform float scaling;
uniform mat4 u_proj;
uniform mat4 u_view;

void main() {
    gl_Position = u_proj * (
        u_view * vec4(a_pos, 1.0) - scaling * vec4(0., 0., 10., 0.)
    );

    frag_pos = (u_view * vec4(a_pos, 1.0)).xyz;
    normal = a_normal;
}
