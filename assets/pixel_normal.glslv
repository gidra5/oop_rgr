#version 150 core
precision highp float;

out vec3 normal;

in ivec3 a_pos;

uniform mat4 u_proj;
uniform mat4 u_view;

void main() {
    mat4 rot = mat4(
        vec4(u_view[0].xyz, 0.),
        vec4(u_view[1].xyz, 0.),
        vec4(u_view[2].xyz, 0.),
        vec4(0., 0., 0., 1.)
    );

    mat4 translate = mat4(
        vec4(1., 0., 0., 0.), 
        vec4(0., 1., 0., 0.),
        vec4(0., 0., 1., 0.),
        u_view[3]
    );

    gl_Position = u_proj * rot * translate * vec4(a_pos, 1.);
    // vec4 p = rot * tm;

    // gl_Position = u_proj * p;

    // gl_Position = u_proj * u_view * vec4(a_pos, 1.0);
    // gl_Position = u_proj * inverse(transpose(rot) * inverse(translate)) * vec4(a_pos, 1.0);
    // gl_Position = u_proj * (u_view[0] * a_pos.x + u_view[1] * a_pos.y + u_view[2] * a_pos.z + u_view[3]);
    // gl_Position = u_proj * (vec4(rot * a_pos, 0.) + u_view[3]);
    // gl_Position = vec4(a_pos, 1.0);
    // normal = vec3(1.);
}