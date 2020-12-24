#version 150 core
precision highp float;

uniform vec3 view_pos;
uniform vec3 light_pos;
uniform vec3 light_color;

uniform sampler2D normal_texture;
uniform sampler2D depth_texture;
uniform mat4 u_proj;
uniform mat4 u_view;
uniform sampler1D vertex_buf[];
uniform sampler1D index_buf[];

in vec3 ray_end;

out vec4 frag_color;

#define PI 3.1415925359
#define TWO_PI 6.2831852
#define MAX_STEPS 200
#define MAX_DIST 1000.
#define MIN_DIST .0001

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q, 0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float dist_scene(vec3 p) 
{
  // vec3 c = vec3(16);
  // p = mod(p + 0.5 * c, c) - 0.5 * c;

  vec4 s = vec4(0, 0, 0, 1.2); //Sphere xyz is position w is radius
  vec4 s2 = vec4(4, 4, 4, 1); //Sphere xyz is position w is radius
  float sphereDist = length(p - s.xyz) - s.w;
  float planeDist  = p.y + 2;
  float boxDist = sdBox(p - s2.xyz, vec3(s2.w));
  return min(min(sphereDist, planeDist), boxDist);
  // return min(sphereDist, boxDist);
  // return min(planeDist, boxDist);
  // return  boxDist;
  // return sphereDist;
  // return planeDist;
  // return sdBox(p - s.xyz, vec3(s.w));
}
 
vec3 dist_scene_gradient(vec3 p)
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
    
    // Shadows
    vec3 p = pos + norm * MIN_DIST * 1.5;
    bool hit = false;

    float d = 0.;
    for (int i = 0; i < MAX_STEPS && d < MAX_DIST && !hit; ++i) {
      float ds = dist_scene(p); 
      d += ds;
      p += light_dir * ds;
      hit = ds < MIN_DIST;
    }
     
    if(d < length(light_pos - p)) diffusion = 0.;
 
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
  float depth = texture(depth_texture, uv).x;

  vec3 p;
  vec3 normal;
  bool hit = false;

  // if (depth != 1.) {
    // vec4 pos_texel = inverse(u_proj) * vec4(uv, depth, 1.);
    // p = pos_texel.xyz / pos_texel.w;
    // normal = texture(normal_texture, uv).xyz;

    // hit = true;
  // } 
  // else {
    vec3 ray_origin = vec3(1, -1, 1) * u_view[3].xyz;
    vec3 ray_dir = normalize(ray_end);

    p = ray_origin;

    float d = 0.;
    for (int i = 0; i < MAX_STEPS && d < MAX_DIST && !hit; ++i) {
      float ds = dist_scene(p); 
      d += ds;
      p += ray_dir * ds;
      hit = ds < MIN_DIST;
    }

    normal = dist_scene_gradient(p);
  // }

  if (hit) {
      frag_color = light(p, normal, vec3(1.)); // Diffuse lighting
  }

  frag_color = vec4(vec3(depth), 1);
}