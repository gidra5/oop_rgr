#version 150 core
precision highp float;

uniform vec3 view_pos;
uniform vec3 light_pos;
uniform vec3 light_color;

uniform sampler2D normal_texture;
uniform sampler2D depth_texture;
uniform mat4 u_proj;
uniform mat4 u_view;
uniform vec2 u_resolution;
uniform sampler1D vertex_buf;
uniform sampler1D index_buf;

in vec3 ray_end;

out vec4 frag_color;

#define PI 3.1415925359
#define TWO_PI 6.2831852
#define MAX_STEPS 100
#define MAX_DIST 1000.
#define MIN_DIST .01

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
  return dot(p, norm);
}

float cylinder(vec3 p, vec3 dir, float r) {
  return length(cross(p, dir)) - r;
}

float triangle(vec3 p, vec3 a, vec3 b, vec3 c, vec3 norm) {
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 ac = a - c; vec3 pc = p - c;

  if (dot(pa, ba) < 0.) {
    if (dot(pa, ac) > 0.) return length(pa);
    else return cylinder(p - c, ac, 0.);
  }

  if (dot(pb, cb) < 0.) {
    if (dot(pb, ba) > 0.) return length(pb);
    else return cylinder(p - a, ba, 0.);
  }

  if (dot(pc, ac) < 0.) {
    if (dot(pc, cb) > 0.) return length(pc);
    else return cylinder(p - b, cb, 0.);
  }

  return abs(plane(p - a, norm));
}

float triangle( vec3 p, vec3 a, vec3 b, vec3 c )
{
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 ac = a - c; vec3 pc = p - c;
  vec3 nor = normalize(cross( ba, ac ));

  if (dot(pa, ba) < 0.) {
    return dot(pa, ac) > 0. ?
      length(pa) :
      cylinder(p - c, ac, 0.);
  }

  if (dot(pb, cb) < 0.) {
    return dot(pb, ba) > 0. ?
      length(pb) :
      cylinder(p - a, ba, 0.);
  }

  if (dot(pc, ac) < 0.) {
    return dot(pc, cb) > 0. ?
      length(pc) :
      cylinder(p - b, cb, 0.);
  }

  return abs(plane(p - a, nor));
}

//can_* show if primitive could be intersected
bool can_plane(vec3 dir, vec3 p, vec3 norm) {
  return dot(dir, norm) * plane(p, norm) < 0.;
}

bool can_cylinder(vec3 dir, vec3 p, vec3 c_dir, float r) {
  return abs(dot(dir, cross(p, c_dir))) < r;
}

bool can_sphere(vec3 dir, vec3 p, float r) {
  return length(cross(p, dir)) < r;
}

float raw_ds(vec3 p) {
  float d = MAX_DIST;

  d = min(d, sphere(p - vec3(-1., 1., -1.), 1.));
  d = min(d, plane(p - vec3(0., -1., 0.), vec3(0., 1., 0.)));
  d = min(d, cylinder(p - vec3(1., 0., 4.), vec3(0., 1., 0.), 1.));

  for (int i = 0; i < 12; ++i) {
    vec4 index = texture(index_buf, float(i) / 12.);
    vec4 a = texture(vertex_buf, index.x / 8.);
    vec4 b = texture(vertex_buf, index.y / 8.);
    vec4 c = texture(vertex_buf, index.z / 8.);
    d = min(d, triangle(p, a.xyz, b.xyz, c.xyz));
  }

  return d;
}

float dist_scene(vec3 p, vec3 dir) {
  float d = MAX_DIST;

  if(can_sphere(dir, p - vec3(-1., 1., -1.), 1.))
    d = min(d, sphere(p - vec3(-1., 1., -1.), 1.));

  if(can_plane(dir, p - vec3(0., -1., 0.), vec3(0., 1., 0.)))
    d = min(d, plane(p - vec3(0., -1., 0.), vec3(0., 1., 0.)));

  if(can_cylinder(dir, p - vec3(1., 0., 4.), vec3(0., 1., 0.), 1.))
    d = min(d, cylinder(p - vec3(1., 0., 4.), vec3(0., 1., 0.), 1.));

  for (int i = 0; i < 12; ++i) {
    vec4 index = texture(index_buf, i);
    vec4 a = texture(vertex_buf, index.x);
    vec4 b = texture(vertex_buf, index.y);
    vec4 c = texture(vertex_buf, index.z);
    d = min(d, triangle(p, a.xyz, b.xyz, c.xyz));
  }

  return d;
}

vec3 dist_scene_gradient(vec3 p, vec3 dir)
{
    float d = raw_ds(p);
    vec2 e = vec2(MIN_DIST,0);
    vec3 n = d - vec3(
      raw_ds(p - e.xyy),
      raw_ds(p - e.yxy),
      raw_ds(p - e.yyx)
    );

    return normalize(n);
}

vec3 raymarch(vec3 origin, vec3 dir, out bool hit) {
  vec3 p = origin;
  hit = false;

  float d = 0.;
  for (int i = 0; i < MAX_STEPS && d < MAX_DIST && !hit; ++i) {
    float ds = dist_scene(p, dir);

    d += ds;
    p += dir * ds;

    hit = ds < MIN_DIST;
  }

  return p;
}

vec4 light(vec3 pos, vec3 norm, vec3 light_pos, vec3 color)
{
    vec3 d_color = vec3(0.);
    float dx = 1.;

    float ambience = 0.02;

    // Shadows
    for (float i = 0.; i < 1.; i += dx) {
      vec3 lp = light_pos + vec3(cos(i * TWO_PI), 0., sin(i * TWO_PI));
      vec3 d = lp - pos;
      float b = dot(d, d);
      vec3 light_dir = d / sqrt(b); // Light Vector

      if (1. > MIN_DIST * b) {
        vec3 p = pos + norm * MIN_DIST * 1.5;
        bool hit;
        p = raymarch(p, light_dir, hit);
        vec3 r = p - pos;

        if( dot(r, r) >= b ) {
          float diffusion = dot(norm, light_dir) / b; // Diffuse light
          d_color += clamp(diffusion, 0., 1.); // Clamp so it doesnt go below 0
        }
      }
    }
    return vec4((d_color * dx + ambience) * light_color * color, 1);
}

void main() {
  vec2 u_resolution = vec2(1280., 720.);
  vec2 uv = gl_FragCoord.xy / u_resolution;

  vec3 obj_color = vec3(1.);
  float depth = texture(depth_texture, uv).x;

  vec3 p;
  vec3 normal;
  bool hit = false;

  // if (false) {
  if (depth != 1.) {
    vec4 pos_texel = inverse(u_proj) * vec4(uv, depth, 1.);
    p = pos_texel.xyz / pos_texel.w;
    normal = texture(normal_texture, uv).xyz;

    hit = true;
  } else {
    vec3 ray_origin = vec3(1, -1, 1) * u_view[3].xyz;
    vec3 ray_dir = normalize(ray_end);

    p = raymarch(ray_origin, ray_dir, hit);

    normal = dist_scene_gradient(p, ray_dir);
  }

  if (hit) {
    frag_color += light(p, normal, light_pos, obj_color); // Diffuse lighting
  } else {
    frag_color = vec4(0.);
  }
}