#version 150 core
precision highp float;

out vec4 frag_color;

in vec3 normal;

void main() {
  // frag_color = vec4(normalize(normal), 1.);
  frag_color = vec4(vec3(gl_FragCoord.z), 1.);
}