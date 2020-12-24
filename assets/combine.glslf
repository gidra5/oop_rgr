#version 150 core

uniform vec3 view_pos;
uniform vec3 light_pos;
uniform vec3 light_color;

uniform sampler2D normal_texture;
uniform sampler2D depth_texture;
uniform mat4 u_proj;
uniform mat4 u_view;
uniform sampler1D vertex_buf[];
uniform sampler1D index_buf[];

in vec3 ray_dir;

out vec4 frag_color;

#define PI 3.1415925359
#define TWO_PI 6.2831852
#define MAX_STEPS 500
#define MAX_DIST 1000.
#define MIN_DIST .0001

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  // if (all(greaterThan(q, vec3(0.)))) {
  //   return length(q);
  // } else {
  //   return max(q.x,max(q.y,q.z));
  // }
  return length(max(q, 0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float dist_scene(vec3 p) 
{
  vec3 c = vec3(16);
  p = mod(p + 0.5 * c, c) - 0.5 * c;

  vec4 s = vec4(0, 0, 0, 1); //Sphere xyz is position w is radius
  // vec4 s2 = vec4(4, 4, 4, 1); //Sphere xyz is position w is radius
  float sphereDist = length(p - s.xyz) - s.w;
  // float planeDist  = p.y + 2;
  // float boxDist = sdBox(p - s.xyz, vec3(s.w));
  // return min(min(sphereDist, planeDist), boxDist);
  // return min(sphereDist, boxDist);
  // return  boxDist;
  return sphereDist;
  // return planeDist;
  // return sdBox(p - s.xyz, vec3(s.w));
}
 
vec3 normal(vec3 p)
{
    float d = dist_scene(p);
    vec2 e = vec2(MIN_DIST,0);
    vec3 n = d - vec3(
      dist_scene(p - e.xyy),
      dist_scene(p - e.yxy),
      dist_scene(p - e.yyx)
    );
 
    return normalize(n);
}
 
vec4 light(vec3 pos, vec3 norm, vec3 color)
{ 
    // Light (directional diffuse)
    vec3 light_dir = normalize(light_pos - pos); // Light Vector
    // vec3 view_dir = normalize(view_pos - frag_pos);
    // vec3 reflect_dir = -reflect(light_dir, norm);

    float ambience = 0.02;
   
    float diffusion = dot(norm, light_dir); // Diffuse light
    diffusion = clamp(diffusion, 0., 1.); // Clamp so it doesnt go below 0

    // float specular_strength = 0.1;
    // float specularity = specular_strength * pow(max(dot(view_dir, reflect_dir), 0.0), 3);
 
    return vec4(
        (diffusion + ambience) * light_color * color,
        1.
    );
}
const vec2 u_resolution = vec2(1280., 720.);

void main() {
  vec2 u_resolution = vec2(1280., 720.);
  vec2 uv = gl_FragCoord.xy / u_resolution;

  vec3 obj_color = vec3(0.5, 0.5, 0.);
  vec4 norm_texel = texture(normal_texture, uv);
  // vec4 pos_texel = inverse(u_proj) * vec4(uv, norm_texel.w, 1.);
  // vec3 pos = pos_texel.xyz / pos_texel.w;
  // vec4 pos_texel = inverse(u_proj) * vec4(uv, 1., 1.);
  vec4 pos_texel = vec4(1,1, norm_texel.w, 1.);

  if (norm_texel.w != 1.) {
    // frag_color = light(pos_texel.xyz / pos_texel.w, norm_texel.xyz, obj_color);
    // frag_color = light(vec3(0), normalize(norm_texel.xyz), obj_color);
    // frag_color = light(vec3(uv,1), norm_texel.xyz, obj_color);
  } 
  uv = (u_resolution - 2. * gl_FragCoord.xy) / u_resolution.y;

  vec3 ray_origin = -view_pos;
  vec3 n_ray_dir = normalize(ray_dir);
  vec4 proj_ro = u_proj * u_view * vec4(ray_origin, 1);
  vec4 proj_rd = u_proj * u_view * vec4(n_ray_dir, 0);

  vec3 p = ray_origin;
  bool hit = false;

  int i = 0;
  float d = 0.;
  // vec4 tmp = proj_ro;
  for (; i < MAX_STEPS && d < MAX_DIST && !hit; ++i) {
  // for (; i < MAX_STEPS && tmp.z < norm_texel.w * tmp.w; ++i) {
  // for (; i < MAX_STEPS && p.z * pos_texel.w < length(pos_texel.xyz - view_pos * pos_texel.w); ++i) {
  // for (; i < MAX_STEPS && d < (1 - norm_texel.w) * MIN_DIST + norm_texel.w * MAX_DIST; ++i) {
    float ds = dist_scene(p); 
    // vec3 t = (p + sign(n_ray_dir) * 8) / n_ray_dir;
    // ds = min(ds, min(min(t.x, t.y), t.z));
    d += ds;
    p += n_ray_dir * ds;
    hit = ds < MIN_DIST;
    // tmp += proj_rd * ds;
  }

  if (hit) {
    // if (norm_texel.w != 1.) {
    //   if (d < -pos_texel.z)
    //   // if (pos_texel.z / pos_texel.w < 0)
    //   // if (d < (1 - norm_texel.w) * MIN_DIST + norm_texel.w * MAX_DIST)
    //     frag_color = vec4(vec3(4. * d / MAX_DIST), 1.);
    //     // frag_color = (vec4(vec3(4. * d / MAX_DIST), 1.) + frag_color) / 2.; // Diffuse lighting
    // //   // frag_color = 16. * light(p, normal(p), vec3(1.)) / dot(ro - p, ro - p); // Diffuse lighting
      frag_color = light(p, normal(p), vec3(1.)); // Diffuse lighting
    // } else {
    // // frag_color = 16. * light(p, normal(p), vec3(1.)) / dot(ro - p, ro - p); // Diffuse lighting
        // frag_color = vec4(vec3(4. * d / MAX_DIST), 1.);
    // frag_color = (light(p, normal(p), vec3(1.)) + frag_color) / 2.; // Diffuse lighting
    // frag_color = vec4(1.);
    // }
  }

  // frag_color = norm_texel;
  // frag_color = texture(depth_texture, uv);
}