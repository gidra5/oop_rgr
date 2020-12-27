#version 150 core
precision highp float;

uniform vec3 light_pos;
uniform vec3 light_color;

uniform mat4 u_proj;
uniform mat4 u_view;
uniform vec2 u_resolution;

const vec3 sphere_center = vec3(0.);
const vec3 plane_center = vec3(0., -1., 0.);
const vec3 cylinder_center = vec3(1., 0., 4.);
const mat3 triangle_pts = mat3(
  vec3(0.),
  vec3(0., 1., 0.),
  vec3(1., 1., 0.)
);

in vec3 ray_end;

out vec4 frag_color;

#define PI 3.1415925359
#define TWO_PI 6.2831852
#define MAX_STEPS 200
#define MAX_DIST 1000.
#define MIN_DIST .0001
#define SAMPLES 1

const float dx = 1. / SAMPLES;

//primitives are "centered" at (0, 0, 0)
float box( vec3 p, vec3 half_sides )
{
  vec3 q = abs(p) - half_sides;
  return length(max(q, 0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sphere(vec3 p, float r) {
  return length(p) - r;
}

float plane(vec3 p, vec3 norm) {
  return abs(dot(p, norm)) - MIN_DIST;
}

float cylinder(vec3 p, vec3 dir, float r) {
  return length(cross(p, dir)) - r;
}

float triangle( vec3 p, vec3 a, vec3 b, vec3 c )
{
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
      if (n.y)  return abs(dot(t[2], ac)) - MIN_DIST;  // directly above triangle
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

bool can_cylinder(vec3 dir, vec3 p, vec3 c_dir, float r) {
  return abs(dot(dir, cross(p, c_dir))) < r;
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
  vec3 center = inverse(transpose(A)) * d;

  // return can_sphere(dir, p - center, length(a - center)) && can_plane(dir, p, normalize(cross(A[0], A[2])));
  return can_sphere(dir, p - center, 2.);
  // return can_plane(dir, p, normalize(cross(A[0], A[2])));
  // return true;
}

float raw_ds(vec3 p) {
  float d = MAX_DIST;

  d = min(d, sphere(p - sphere_center, 1.));
  d = min(d, plane(p - plane_center, -normalize(vec3(0., 1., 0.))));
  d = min(d, cylinder(p - cylinder_center, vec3(0., 1., 0.), 1.));

  d = min(d, triangle(p - vec3(2., 0., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]));

  return d;
}

float dist_scene(vec3 p, vec3 dir) {
  float d = MAX_DIST;

  if (can_sphere(dir, p - sphere_center, 1.))
    d = min(d, sphere(p - sphere_center, 1.));

  if (can_plane(dir, p - plane_center, vec3(0., 1., 0.)))
    d = min(d, plane(p - plane_center, vec3(0., 1., 0.)));

  if (can_cylinder(dir, p - cylinder_center, vec3(0., 1., 0.), 1.))
    d = min(d, cylinder(p - cylinder_center, vec3(0., 1., 0.), 1.));

  if (can_triangle(dir, p - vec3(2., 0., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]))
    d = min(d, triangle(p - vec3(2., 0., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]));

  return d;
}

vec3 dist_scene_gradient(vec3 p)
{
    float d = raw_ds(p);
    vec2 e = vec2(MIN_DIST,0);
    vec3 n = d - vec3(
      raw_ds(p + e.xyy),
      raw_ds(p + e.yxy),
      raw_ds(p + e.yyx)
    );

    return -normalize(n);
}

vec3 raymarch(vec3 origin, vec3 dir, out bool hit) {
  vec3 p = origin;
  hit = false;

  float d = 0.;
  for (int i = 0; i < MAX_STEPS && d < MAX_DIST && !hit; ++i) {
    float ds = dist_scene(p, dir);
    // float ds = raw_ds(p);

    d += ds;
    p += dir * ds;

    hit = ds < MIN_DIST;
  }

  return p;
}

vec4 light(vec3 pos, vec3 norm, vec3 light_pos, vec3 color)
{
    float d_color = 0.;

    float ambience = 0.02;

    // Shadows
    for (int i = 0; i < SAMPLES; i += 1) {
      vec3 lp = light_pos + vec3(cos(i * dx * TWO_PI), 0., sin(i * dx * TWO_PI));
      vec3 d = lp - pos;
      float b = dot(d, d);
      vec3 light_dir = d * inversesqrt(b);

      if (1. > MIN_DIST * b) {
        vec3 p = pos + norm * MIN_DIST * 1.5;
        bool hit;
        p = raymarch(p, light_dir, hit);
        vec3 r = p - pos;

        if ( dot(r, r) >= b ) {
          float diffusion = dot(norm, light_dir) / b;
          d_color += max(diffusion, 0.);
        }
      }
    }
    return vec4((d_color * dx + ambience) * light_color * color, 1);
}

void main() {
  vec2 uv = gl_FragCoord.xy / u_resolution;

  vec3 obj_color = vec3(1.);

  vec3 p;
  vec3 normal;
  bool hit = false;
  vec3 ray_origin = vec3(1, -1, 1) * u_view[3].xyz;
  vec3 ray_dir = normalize(ray_end);

  p = raymarch(ray_origin, ray_dir, hit);

  normal = dist_scene_gradient(p);

  if (hit) {
    frag_color += light(p, normal, light_pos, obj_color); // Diffuse lighting
  } else {
    frag_color = vec4(0.);
  }
}