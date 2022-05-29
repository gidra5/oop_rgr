#version 150 core

#define PHI 1.61803398874989484820459 // Golden Ratio   
#define SRT 1.41421356237309504880169 // Square Root of Two
#define PI 3.14159265358979323846264
#define TWO_PI 6.28318530717958647692528

precision highp float;

uniform mat4 u_proj;
uniform mat4 u_view;
uniform vec2 u_resolution;
uniform samplerCube skybox;
uniform float t;
uniform uint samples;

const   float max_dist  = 1000.;
const   float min_dist  = 0.0005;
const   int   max_steps = 1000;

uniform vec3 light_pos;
uniform vec3 light_color;
uniform vec3 sphere_center;
uniform vec3 plane_center;
uniform vec3 cylinder_center;
const mat3 triangle_pts = mat3(
  vec3(0.),
  vec3(0., 1., 0.),
  vec3(1., 1., 0.)
);

struct Ray {
  vec3 pos; // Origin
  vec3 dir; // Direction (normalized)
};

in vec3 ray_dir;

out vec4 frag_color;

float random_0t1(in vec2 coordinate, in float seed) {
  int base = 2<<8;
  int modulo = 2<<16;
  return fract(sin(dot(coordinate * (fract(seed / modulo) * modulo + base), vec2(PHI * .1, PI * .1))) * SRT * 10000.0);
  // return fract(sin(dot(coordinate * seed, vec2(PHI * .1, PI * .1))) * SRT * 10000.0);
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

//primitives are "centered" at (0, 0, 0)
float box(vec3 p, vec3 half_sides) {
  vec3 q = abs(p) - half_sides;
  return length(max(q, 0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sphere(vec3 p, float r) {
  return length(p) - r;
}
float sphere(Ray ray, float r) {
  float t = dot(ray.dir, ray.pos);
  float r_min = r + min_dist;
  float d = t + sqrt(r_min * r_min - dot(ray.pos, ray.pos) + t * t);
  // if (d < 0.) return sphere(p, r);
  // else        return d;
  if (t < 0.) return max_dist; 
  else return d;
}

float plane(vec3 p, vec3 norm) {
  return abs(dot(p, norm)) - min_dist;
}

float plane(Ray ray, vec3 norm) {
  float d1 = -dot(ray.pos, norm);
  float d2 = d1 + min_dist;
  float d3 = dot(ray.dir, norm);
  if (d1 * d3 < 0.) return max_dist;
  else              return d2 / d3;
}

float cylinder(vec3 p, vec3 dir_c, float r) {
  return length(cross(p, dir_c)) - r;
}

float triangle(vec3 p, vec3 a, vec3 b, vec3 c) {
  vec3 ba = normalize(b - a); vec3 pa = p - a;
  vec3 cb = normalize(c - b); vec3 pb = p - b;
  vec3 ac = normalize(a - c); vec3 pc = p - c;
  mat3 t = mat3(
    cross(pc, ac),
    cross(pb, cb),
    cross(pa, ba)
  );

  vec3 d = transpose(t) * cross(ac, ba);

  ivec3 n = ivec3(lessThan(d, vec3(0.)));
  ivec3 not_n = 1 - n;
  int f = int(dot(n, n));
  if (f == 3) return abs(dot(t[2], ac)) - min_dist; // directly above triangle

  mat3 vecs = mat3(pa, pb, pc);
  int index = int(dot(f == 2 ? not_n : n, ivec3(0, 1, 2)));
  return length(f == 2 ? t[index] : vecs[index]); 
}

//can_* show if primitive could be intersected
bool can_plane(Ray ray, vec3 norm) {
  return dot(ray.dir, norm) * dot(ray.pos, norm) < 0.;
}

bool can_cylinder(Ray ray, vec3 dir_c, float r) {
  return abs(dot(ray.dir, cross(ray.pos, dir_c))) < r;
}

bool can_sphere(Ray ray, float r) {
  return length(cross(ray.pos, ray.dir)) < r;
}

bool can_triangle(Ray ray, vec3 a, vec3 b, vec3 c) {
  mat3 A = mat3(
    b - a,
    c - b,
    a - c
  );
  vec3 d = vec3(
    dot(b, b) - dot(a, a),
    dot(c, c) - dot(b, b),
    dot(a, a) - dot(c, c)
  );
  vec3 center = 2 * inverse(transpose(A)) * d;

  // return can_sphere(dir, p - center, length(a - center)) && can_plane(dir, p, normalize(cross(A[0], A[2])));
  // return can_sphere(dir, p - center, length(a - center));
  return can_plane(ray, normalize(cross(A[0], A[2])));
  // return can_cylinder(dir, p - point_ba, , max(dot(A[0], A[2]), max(dot(A[0], A[1]), dot(A[0], A[0]))) / (2. * length(A[0])));
  // return true;
}

float raw_ds(vec3 p) {
  float d = max_dist;

  d = min(d, sphere(p - sphere_center, 1.));
  d = min(d, plane(p - plane_center, -normalize(vec3(0., 1., 0.))));
  d = min(d, cylinder(p - cylinder_center, vec3(0., 1., 0.), 1.));

  d = min(d, triangle(p - vec3(2., 2., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]));

  return d;
}

float dist_scene(Ray ray) {
  vec3 p = ray.pos;
  vec3 dir = ray.dir;
  float d = max_dist;

  // d = min(d, sphere(p - sphere_center, dir, 1.));
  if (can_sphere(Ray(p - sphere_center, dir), 1.))
    d = min(d, sphere(p - sphere_center, 1.));

  d = min(d, plane(Ray(p - plane_center, dir), vec3(0., 1., 0.)));
  // if (can_plane(p - plane_center, dir, vec3(0., 1., 0.)))
  // d = min(d, plane(p - plane_center, vec3(0., 1., 0.)));

  if (can_cylinder(Ray(p - cylinder_center, dir), vec3(0., 1., 0.), 1.))
    d = min(d, cylinder(p - cylinder_center, vec3(0., 1., 0.), 1.));

  // if (can_triangle(dir, p - vec3(2., 0., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]))
  //   d = min(d, triangle(p - vec3(2., 0., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]));

  return d;
}

vec3 dist_scene_gradient(vec3 p) {
    float d = raw_ds(p);
    vec2 e = vec2(min_dist,0);

    return (vec3(raw_ds(p + e.xyy), raw_ds(p + e.yxy), raw_ds(p + e.yyx)) - d) / min_dist;
}

vec3 raymarch(Ray ray, out bool hit) {
  vec3 origin = ray.pos;
  vec3 dir = ray.dir;
  vec3 p = origin;
  hit = false;

  float d = 0.;
  for (int steps = 0; steps < max_steps && d < max_dist && !hit; ++steps) {
    float ds = dist_scene(Ray(p, dir));
    // float ds = raw_ds(p);

    d += ds;
    p += dir * ds;

    hit = ds < min_dist;
  }

  return p;
}

vec3 sample_light_shape(float t) {
  return light_pos + vec3(cos(t * TWO_PI), 0., sin(t * TWO_PI));
}

// float in_shadow_simple(vec3 pos, vec3 lp) {
//   vec3 d = lp - pos;
//   float mag_sq = dot(d, d);
//   vec3 light_dir = d * inversesqrt(mag_sq);

//   vec3 p = pos;
//   bool hit = false;
//   float d = 0.;

//   for (int i = 0; i < max_steps && !hit; ++i) {
//     float ds = dist_scene(p, light_dir);

//     d += ds;
//     p += light_dir * ds;

//     hit = ds < min_dist;
//     if ( d * d >= mag_sq ) {
//       return 1.;
//     }
//   }

//   return 0.;
// }

vec3 in_shadow(vec3 pos, vec3 light_dir, float mag_sq) {
  if (1. > min_dist * min_dist * mag_sq) {
    // float sin_a = 1.;
    vec3 p = pos;
    bool hit = false;
    float d = 0.;

    for (int i = 0; i < max_steps && !hit; ++i) {
      float ds = dist_scene(Ray(p, light_dir));
      // float ds = raw_ds(p);

      d += ds;
      p += light_dir * ds;

      hit = ds < min_dist;
      if ( d * d >= mag_sq ) {
        // return light_dir * sin_a;
        // return light_dir * sin_a / mag_sq;
        return light_dir;
      }

      // sin_a = min(sin_a, ds / d);
    }
  }

  return vec3(0.);
}

vec3 sun_color = vec3(0x92, 0x97, 0xC4) / 0xff * 0.9;
vec3 light(vec3 pos, vec3 norm, vec3 light_pos) {
    vec3 d_color = vec3(0.);
    vec2 e = vec2(min_dist, 0);

    float ambience = 0.02;
    // float pos_in_shadow = in_shadow_simple(pos, light_pos);

    vec3 p = pos + 1.1 * norm * min_dist;
    vec3 d = light_pos - p;
    float mag_sq = dot(d, d);
    d_color += in_shadow(p, normalize(d), mag_sq) / mag_sq;

    // d_color.x += (in_shadow_simple(pos, light_pos + e.xyy) - pos_in_shadow) * sample_light_shape(t)

    return max(dot(d_color, norm) / samples, 0.) * light_color + ambience * sun_color;
}
vec3 sun(vec3 pos, vec3 norm) {
    vec3 d_color = vec3(0.);
    vec2 e = vec2(min_dist, 0);

    float ambience = 0.02;
    vec3 ambiance_color = sun_color;
    // float pos_in_shadow = in_shadow_simple(pos, light_pos);

    // Shadows
    vec3 p = pos + 1.1 * norm * min_dist;
    vec3 d = vec3(1.);
    d_color += in_shadow(p, normalize(d), 0.99 / (min_dist * min_dist));

    // d_color.x += (in_shadow_simple(pos, light_pos + e.xyy) - pos_in_shadow) * sample_light_shape(t)

    return (max(dot(d_color, norm) / samples, 0.) + ambience) * sun_color;
}

float cameraFovAngle = PI * 2. / 3.;
float paniniDistance = 0.75;
float verticalCompression = 0.1;
float halfFOV = cameraFovAngle / 2.f;
vec2 p = vec2(sin(halfFOV), cos(halfFOV) + paniniDistance);
float M = sqrt(dot(p, p));
float halfPaniniFOV = atan(p.x, p.y);

vec3 pinholeRay(vec2 pixel) { 
  return vec3(pixel, 1/tan(halfFOV));
}

vec3 paniniRay(vec2 pixel) {
  vec2 hvPan = pixel * vec2(halfPaniniFOV, halfFOV);
  float x = sin(hvPan.x) * M;
  float z = cos(hvPan.x) * M - paniniDistance;
  float y = tan(hvPan.y) * (z + verticalCompression);

  return vec3(x, y, z);
}

float imagePlaneDistance = 1.;
float lensFocalLength = 2.1;
float circleOfConfusionRadius = 0.04;

Ray thinLensRay(vec3 ray, vec2 lensOffset) {
  float theta = lensOffset.x * TWO_PI;
  float radius = sqrt(lensOffset.y);
  vec2 uv = vec2(cos(theta), sin(theta)) * radius;

  float focusPlane = (imagePlaneDistance * lensFocalLength) / (imagePlaneDistance - lensFocalLength);
  vec3 focusPoint = ray * (focusPlane / ray.z);

  vec3 origin = vec3(uv * circleOfConfusionRadius, 0.f);
  vec3 direction = -normalize(focusPoint + origin);
  return Ray(origin, direction);
}

void main() {
  for (int x = 0; x < int(samples); ++x) {
    vec2 subpixel = random_0t1_2(gl_FragCoord.xy, t * x + 0.5);
    vec2 uv = (2. * (gl_FragCoord.xy + subpixel) - u_resolution) / u_resolution.x;
    // vec3 rayDirection = normalize(paniniRay(uv));
    vec3 rayDirection = normalize(pinholeRay(uv));
    vec3 subpixel_color = vec3(0.);

    for (int x = 0; x < int(samples) * 4; ++x) {
      Ray ray = thinLensRay(rayDirection, normalize(random_0t1_2(uv, t * x)));
      ray.dir = (u_view * vec4(normalize(ray.dir), 0.)).xyz;
      vec4 ray_pos = u_view * vec4(ray.pos, 1.);
      ray.pos = ray_pos.xyz / ray_pos.w;

      vec3 obj_color = vec3(1.);

      bool hit = false;

      vec3 p = raymarch(ray, hit);

      vec3 normal = normalize(dist_scene_gradient(p));

      if (hit) {
        vec3 color = sun(p, normal) * obj_color;

        for (int x = 0; x < int(samples); ++x) {
          vec3 light_pos = sample_light_shape(random_0t1(uv, t * x));
          color += light(p, normal, light_pos) * obj_color / samples; // Diffuse lighting
        }

        subpixel_color += color / (int(samples) * 4); 
      } else {
        // subpixel_color = texture(skybox, ray.dir);
        subpixel_color = vec3(0.);
      }
    }
    frag_color += vec4(subpixel_color, 1);
  }
  frag_color /= samples;
}