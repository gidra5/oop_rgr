#version 150 core
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

struct Ray {
  vec3 pos; // Origin
  vec3 dir; // Direction (normalized)
};

in vec3 ray_dir;

out vec4 frag_color;

#define PI 3.1415925359
#define TWO_PI 6.2831852

const float dx = 1. / samples;

//primitives are "centered" at (0, 0, 0)
float box(vec3 p, vec3 half_sides) {
  vec3 q = abs(p) - half_sides;
  return length(max(q, 0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sphere(vec3 p, float r) {
  return length(p) - r;
}
float sphere(vec3 p, vec3 dir, float r) {
  float t = dot(dir, p);
  float d = t + sqrt((r + min_dist) * (r + min_dist) - dot(p, p) + t * t);
  // if (d < 0.) return sphere(p, r);
  // else        return d;
  if (t < 0.) return max_dist; 
  else return d;
}

float plane(vec3 p, vec3 norm) {
  return abs(dot(p, norm)) - min_dist;
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

float triangle( vec3 p, vec3 a, vec3 b, vec3 c ) {
  vec3 ba = normalize(b - a); vec3 pa = p - a;
  vec3 cb = normalize(c - b); vec3 pb = p - b;
  vec3 ac = normalize(a - c); vec3 pc = p - c;
  mat3 t = mat3(
    cross(pc, ac),
    cross(pb, cb),
    cross(pa, ba)
  );

  vec3 d = transpose(t) * cross( ac, ba );

  bvec3 n = lessThan(d, vec3(0.));

  if (n.x) {
    if (n.z) {
      if (n.y)  return abs(dot(t[2], ac)) - min_dist;  // directly above triangle
      else      return length(t[1]);        // near side cb?
    } else {
      if (n.y)  return length(t[2]);        // near side ba?
      else      return length(pb);          // near b
    }
  } else {
    if (n.z) {
      if (n.y)  return length(t[0]);        // near side ac?
      else      return length(pc);          // near c
    } else      return length(pa);          // near a
  }
}

//can_* show if primitive could be intersected
bool can_plane(vec3 dir, vec3 p, vec3 norm) {
  return dot(dir, norm) * dot(p, norm) < 0.;
}

bool can_cylinder(vec3 dir, vec3 p, vec3 dir_c, float r) {
  return abs(dot(dir, cross(p, dir_c))) < r;
}

bool can_sphere(vec3 dir, vec3 p, float r) {
  return length(cross(p, dir)) < r;
}

bool can_triangle(vec3 dir, vec3 p, vec3 a, vec3 b, vec3 c) {
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
  return can_plane(dir, p, normalize(cross(A[0], A[2])));
  // return can_cylinder(dir, p - point_ba, , max(dot(A[0], A[2]), max(dot(A[0], A[1]), dot(A[0], A[0]))) / (2. * length(A[0])));
  // return true;
}

float raw_ds(vec3 p) {
  float d = max_dist;

  d = min(d, sphere(p - sphere_center, 1.));
  d = min(d, plane(p - plane_center, -normalize(vec3(0., 1., 0.))));
  d = min(d, cylinder(p - cylinder_center, vec3(0., 1., 0.), 1.));

  // d = min(d, triangle(p - vec3(2., 0., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]));

  return d;
}

float dist_scene(vec3 p, vec3 dir) {
  float d = max_dist;

  // d = min(d, sphere(p - sphere_center, dir, 1.));
  if (can_sphere(dir, p - sphere_center, 1.))
    d = min(d, sphere(p - sphere_center, 1.));

  d = min(d, plane(p - plane_center, dir, vec3(0., 1., 0.)));
  // if (can_plane(p - plane_center, dir, vec3(0., 1., 0.)))
  // d = min(d, plane(p - plane_center, vec3(0., 1., 0.)));

  if (can_cylinder(dir, p - cylinder_center, vec3(0., 1., 0.), 1.))
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

vec3 raymarch(vec3 origin, vec3 dir, out bool hit) {
  vec3 p = origin;
  hit = false;

  float d = 0.;
  for (int steps = 0; steps < max_steps && d < max_dist && !hit; ++steps) {
    float ds = dist_scene(p, dir);
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

vec3 in_shadow(vec3 pos, vec3 lp) {
  vec3 d = lp - pos;
  float mag_sq = dot(d, d);
  vec3 light_dir = d * inversesqrt(mag_sq);

  if (1. > min_dist * min_dist * mag_sq) {
    // float sin_a = 1.;
    vec3 p = pos;
    bool hit = false;
    float d = 0.;

    for (int i = 0; i < max_steps && !hit; ++i) {
      float ds = dist_scene(p, light_dir);
      // float ds = raw_ds(p);

      d += ds;
      p += light_dir * ds;

      hit = ds < min_dist;
      if ( d * d >= mag_sq ) {
        // return light_dir * sin_a;
        // return light_dir * sin_a / mag_sq;
        return light_dir / mag_sq;
      }

      // sin_a = min(sin_a, ds / d);
    }
  }

  return vec3(0.);
}

vec4 light(vec3 pos, vec3 norm, vec3 light_pos, vec3 color) {
    vec3 d_color = vec3(0.);
    vec2 e = vec2(min_dist,0);

    float ambience = 0.02;
    // float pos_in_shadow = in_shadow_simple(pos, light_pos);

    // Shadows
    for (float t = 0.; t < 1.; t += dx) {
      d_color += in_shadow(pos + 1.1 * norm * min_dist, sample_light_shape(t));
    }

    // d_color.x += (in_shadow_simple(pos, light_pos + e.xyy) - pos_in_shadow) * sample_light_shape(t)

    return vec4((max(dot(d_color, norm) * dx, 0.) + ambience) * light_color * color, 1);
}

float cameraFovAngle = PI * 2. / 3.;
float paniniDistance = 0.75;
float verticalCompression = 0.1;
float halfFOV = cameraFovAngle / 2.f;
vec2 p = vec2(sin(halfFOV), cos(halfFOV) + paniniDistance);
float M = sqrt(dot(p, p));
float halfPaniniFOV = atan(p.x, p.y);

vec3 paniniRay(vec2 pixel) {
  vec2 hvPan = pixel * vec2(halfPaniniFOV, halfFOV);
  float x = sin(hvPan.x) * M;
  float z = cos(hvPan.x) * M - paniniDistance;
  float y = tan(hvPan.y) * (z + verticalCompression);

  return vec3(x, y, z);
}

float imagePlaneDistance = 2.;
float lensFocalLength = 5.;
float focalLength = 1.;
float fStop = 50.;

Ray thinLensRay(vec3 ray, vec2 lensOffset) {
  float theta = lensOffset.x * TWO_PI;
  float radius = sqrt(lensOffset.y);
  vec2 uv = vec2(cos(theta), sin(theta)) * radius;

  float focusPlane = (imagePlaneDistance * lensFocalLength) / (imagePlaneDistance - lensFocalLength);
  vec3 focusPoint = ray * (focusPlane / ray.z);
  float circleOfConfusionRadius = focalLength / (2.f * fStop);

  vec3 origin = vec3(uv * circleOfConfusionRadius, 0.f);
  vec3 direction = -normalize(focusPoint + origin);
  return Ray(origin, direction);
}

float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

void main() {
  vec2 uv = (2. * gl_FragCoord.xy - u_resolution) / u_resolution.x;

  mat3 rot = mat3(u_view[0].xyz, u_view[1].xyz, u_view[2].xyz);
  vec3 paniniRayDirection = paniniRay(uv);
  // vec3 paniniRayDirection = ray_dir;

  for (int x = 0; x < 64; ++x) {
    Ray ray = thinLensRay(paniniRayDirection, vec2(rand(vec2(1. / x)), rand(vec2(1. / x + 1.))));
    ray.dir = (u_view * vec4(normalize(ray.dir), 0.)).xyz;
    vec4 ray_pos = u_view * vec4(ray.pos, 1.);
    ray.pos = ray_pos.xyz / ray_pos.w;

    vec3 obj_color = vec3(1.);

    bool hit = false;

    vec3 p = raymarch(ray.pos, ray.dir, hit);

    vec3 normal = normalize(dist_scene_gradient(p));

    if (hit) {
      frag_color += light(p, normal, light_pos, obj_color); // Diffuse lighting
    } else {
      frag_color = vec4(0.);
    }
  }
  frag_color /= 64;
}