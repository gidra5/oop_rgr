#version 150 core
precision highp float;

#define PI 3.1415925359
#define TWO_PI 6.2831852

in ivec2 a_pos;
// out vec3 ray_dir;
// uniform mat4 u_view;
// uniform mat4 u_proj;
// uniform vec2 u_resolution;

// float cameraFovAngle = PI * 2. / 3.;
// float paniniDistance = 0.75;
// float verticalCompression = 0.1;
// float halfFOV = cameraFovAngle / 2.f;
// vec2 p = vec2(sin(halfFOV), cos(halfFOV) + paniniDistance);
// float M = sqrt(dot(p, p));
// float halfPaniniFOV = atan(p.x, p.y);

// vec3 paniniRay(vec2 pixel) {
//   vec2 hvPan = pixel * vec2(halfPaniniFOV, halfFOV);
//   float x = sin(hvPan.x) * M;
//   float z = cos(hvPan.x) * M - paniniDistance;
//   float y = tan(hvPan.y) * (z + verticalCompression);

//   return vec3(x, y, z);
// }

void main() {
  gl_Position = vec4(a_pos, 0., 1.);
  
  // ray_dir = paniniRay(a_pos);
}