#version 150 core
precision highp float;

in ivec2 a_pos;

void main() {
  gl_Position = vec4(a_pos, 0., 1.);
}