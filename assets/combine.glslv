#version 150 core
precision highp float;

#define PI 3.1415925359
#define TWO_PI 6.2831852

in ivec2 a_pos;

void main() {
  gl_Position = vec4(a_pos, 0., 1.);
}