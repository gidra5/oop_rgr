#version 150 core

in ivec2 a_pos;
out vec3 ray_dir;
uniform mat4 u_view;
uniform mat4 u_proj;

void main() {
  gl_Position = vec4(a_pos, 0., 1.);
  
  mat4 inv = inverse(u_proj * u_view);
  vec4 r1 = inv * vec4(a_pos, 1., 1.);
  vec4 r2 = inv * vec4(a_pos, 0., 1.);
  ray_dir = r1.xyz / r1.w - r2.xyz / r2.w;
}