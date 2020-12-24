#version 150 core

out vec4 frag_color;

in vec3 normal;

void main() {
  frag_color = vec4(normalize(normal), gl_FragCoord.z);
}