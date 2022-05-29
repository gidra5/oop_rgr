#version 150 core

#define PHI 1.61803398874989484820459 // Golden Ratio   
#define SRT 1.41421356237309504880169 // Square Root of Two
#define PI 3.14159265358979323846264
#define TWO_PI 6.28318530717958647692528

precision highp float;

uniform sampler2D prev;
uniform sampler2D current;
uniform float dt;
uniform float t;
// uniform uint samples;

out vec4 frag_color;

float random_0t1(in vec2 coordinate, in float seed) {
  int base = 2<<8;
  int modulo = 2<<16;
  return fract(sin(dot(coordinate * (fract(seed / modulo) * modulo + base), vec2(PHI, PI)) * .1) * SRT * 10000.0);
}
vec2 random_0t1_2(in vec2 coordinate, in float seed) {
  return vec2(random_0t1(coordinate, seed), random_0t1(coordinate, seed * 0.5 + 3.));
}
vec3 random_0t1_3(in vec2 coordinate, in float seed) {
  return vec3(random_0t1_2(coordinate, seed), random_0t1(coordinate, seed * 0.75 + 2.));
}
vec4 random_0t1_4(in vec2 coordinate, in float seed) {
  return vec4(random_0t1_3(coordinate, seed), random_0t1(coordinate, seed * 0.85 + 1.));
}

void main() {
  // frag_color = frag_color * vec4(0.9, 1., 1., 1.);
  // frag_color = texture(current, gl_FragCoord.xy / vec2(1280, 720));
  // frag_color = frag_color;
  // frag_color = gl_FragCoord;
  // frag_color = vec4(gl_FragCoord.xy / vec2(1280, 720), 0., 1.);
}