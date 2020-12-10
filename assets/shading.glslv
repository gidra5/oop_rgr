#version 150 core

out float diffusion;
out vec3 frag_pos;
out vec3 normal;

in ivec3 a_pos;
in ivec3 a_normal;

uniform mat4 u_proj;
uniform mat4 u_view;
uniform mat4 u_model;
// uniform vec3 view_pos;
uniform vec3 light_pos;
// uniform vec3 light_color;

void main() {
    vec4 pos = u_model * vec4(a_pos, 1.0);
    gl_Position = u_view * pos;
    gl_Position.z -= 8.;
    gl_Position = u_proj * gl_Position;

    frag_pos = pos.xyz;
    normal = a_normal;

    vec3 obj_color = vec3(0.5, 0., 0.);

    vec3 norm = normalize(normal);
    vec3 light_dir = normalize(light_pos - frag_pos);

    diffusion = max(dot(norm, light_dir), 0.0);
}
