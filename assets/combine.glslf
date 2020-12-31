#version 150 core

#define PI 3.1415925359
#define TWO_PI 6.2831852

precision highp float;

uniform mat4 u_proj;
uniform mat4 u_view;
uniform vec2 u_resolution;

const   float max_dist  = 1000.;
const   float min_dist  = 0.005;
const   int   samples   = 1;
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

in vec3 ray_dir;

out vec4 frag_color;

const float dx = 1. / samples;

float sphere(vec3 p, float r);
float plane(vec3 p, vec3 dir, vec3 norm);
float cylinder(vec3 p, vec3 dir_c, float r)

//can_* show if primitive could be intersected
bool can_cylinder(vec3 dir, vec3 p, vec3 dir_c, float r);
bool can_sphere(vec3 dir, vec3 p, float r);

float raw_ds(vec3 p);
float dist_scene(vec3 p, vec3 dir);
vec3 dist_scene_gradient(vec3 p);
vec3 raymarch(vec3 origin, vec3 dir, out bool hit);
vec3 sample_light_shape(float t);
vec3 in_shadow(vec3 pos, vec3 lp);
vec4 light(vec3 pos, vec3 norm, vec3 light_pos, vec3 color);

void main() {
  vec2 uv = gl_FragCoord.xy / u_resolution;

  vec3 obj_color = vec3(1.);

  bool hit = false;
  vec3 ray_origin = vec3(1, -1, 1) * u_view[3].xyz;
  vec3 ray_dir_n = normalize(ray_dir);

  vec3 p = raymarch(ray_origin, ray_dir_n, hit);

  vec3 normal = normalize(dist_scene_gradient(p));

  if (hit) {
    frag_color += light(p, normal, light_pos, obj_color);
  } else {
    frag_color = vec4(0.);
  }
}

float raw_ds(vec3 p) {
  float d = max_dist;

  d = min(d, sphere(p - sphere_center, 1.));
  d = min(d, plane(p - plane_center, -normalize(vec3(0., 1., 0.))));
  d = min(d, cylinder(p - cylinder_center, vec3(0., 1., 0.), 1.));

  return d;
}

//primitives are "centered" at (0, 0, 0)
float sphere(vec3 p, float r) {
  return length(p) - r;
}

float plane(vec3 p, vec3 dir, vec3 norm) {
  float d1 = -dot(p, norm);
  float d2 = d1 + min_dist;
  float d3 = dot(dir, norm);
  if (d1 * d3 < 0.) return max_dist;
  else              return d2 / d3;
}

float cylinder(vec3 p, vec3 dir_c, float r) {
  return length(cross(p, dir_c)) - r;
}

//can_* show if primitive could be intersected
bool can_cylinder(vec3 dir, vec3 p, vec3 dir_c, float r) {
  return abs(dot(dir, cross(p, dir_c))) < r;
}

bool can_sphere(vec3 dir, vec3 p, float r) {
  return length(cross(p, dir)) < r;
}


float dist_scene(vec3 p, vec3 dir) {
  float d = max_dist;

  if (can_sphere(dir, p - sphere_center, 1.))
    d = min(d, sphere(p - sphere_center, 1.));

  d = min(d, plane(p - plane_center, dir, vec3(0., 1., 0.)));

  if (can_cylinder(dir, p - cylinder_center, vec3(0., 1., 0.), 1.))
    d = min(d, cylinder(p - cylinder_center, vec3(0., 1., 0.), 1.));

  return d;
}

vec3 dist_scene_gradient(vec3 p) {
    float d = raw_ds(p);
    vec2 e = vec2(min_dist,0);

    return (vec3(raw_ds(p + e.xyy), raw_ds(p + e.yxy), raw_ds(p + e.yyx)) - d) / min_dist;
}

vec3 raymarch(vec3 origin, vec3 dir, out bool hit) {
  vec3 p = origin;
  hit = false;

  float d = 0.;
  for (int steps = 0; steps < max_steps && d < max_dist && !hit; ++steps) {
    float ds = dist_scene(p, dir);

    d += ds;
    p += dir * ds;

    hit = ds < min_dist;
  }

  return p;
}

vec3 sample_light_shape(float t) {
  return light_pos + vec3(cos(t * TWO_PI), 0., sin(t * TWO_PI));
}

vec3 in_shadow(vec3 pos, vec3 lp) {
  vec3 d = lp - pos;
  float mag_sq = dot(d, d);
  vec3 light_dir = d * inversesqrt(mag_sq);

  if (1. > min_dist * min_dist * mag_sq) {
    vec3 p = pos;
    bool hit = false;
    float d = 0.;

    for (int i = 0; i < max_steps && !hit; ++i) {
      float ds = dist_scene(p, light_dir);

      d += ds;
      p += light_dir * ds;

      hit = ds < min_dist;
      if ( d * d >= mag_sq ) {
        return light_dir / mag_sq;
      }
    }
  }

  return vec3(0.);
}

vec4 light(vec3 pos, vec3 norm, vec3 light_pos, vec3 color) {
    vec3 d_color = vec3(0.);

    float ambience = 0.05;

    // Shadows
    for (float t = 0.; t < 1.; t += dx) {
      d_color += in_shadow(pos + 1.1 * norm * min_dist, sample_light_shape(t));
    }

    return vec4((max(dot(d_color, norm) * dx, 0.) + ambience) * light_color * color, 1);
}