#version 150 core
precision highp float;

in ivec2 a_pos;
out vec3 ray_dir;
uniform mat4 u_view;
uniform mat4 u_proj;

void main() {
  gl_Position = vec4(a_pos, 0., 1.);
  // mat4 rot = mat4(
  //     vec4(u_view[0].xyz, 0.),
  //     vec4(u_view[1].xyz, 0.),
  //     vec4(u_view[2].xyz, 0.),
  //     vec4(0., 0., 0., 1.)
  // );
  mat4 rot = mat4(mat3(u_view[0].xyz, u_view[1].xyz, u_view[2].xyz));
  vec4 r1 = transpose(rot) * inverse(u_proj) * vec4(a_pos, 0., u_view[3].w);
  
  ray_dir = vec3(-1, 1, -1) * r1.xyz / r1.w;
}