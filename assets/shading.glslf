#version 150 core

vec2 u_resolution;

in float diffusion;
in vec3 frag_pos;
in vec3 normal;

out vec4 frag_color;
uniform vec3 view_pos;
uniform vec3 light_pos;
uniform vec3 light_color;

// void main() {
//     vec3 obj_color = vec3(0.0, 0.5, 0.0);

//     vec3 norm = normalize(normal);
//     vec3 light_dir = normalize(light_pos - frag_pos);
//     vec3 view_dir = normalize(view_pos - frag_pos);
//     vec3 reflect_dir = reflect(-light_dir, norm);

//     // float diffusion = max(dot(norm, light_dir), 0.0);
//     float ambience = 0.1;

//     float specular_strength = 0.1;
//     float specularity = specular_strength * pow(max(dot(view_dir, reflect_dir), 0.0), 3);

//     frag_color = vec4(
//         (ambience + diffusion + specularity) * light_color * obj_color,
//         1.
//     );
// }


// Constants
#define PI 3.1415925359
#define TWO_PI 6.2831852
#define MAX_STEPS 100
#define MAX_DIST 100.
#define SURFACE_DIST .01
 
float GetDist(vec3 p) 
{
    vec4 s = vec4(0,1,6,1); //Sphere xyz is position w is radius
    float sphereDist = length(p-s.xyz) - s.w;
    float planeDist  = p.y;
    float d = min(sphereDist,planeDist);
    return d;
}
 
float RayMarch(vec3 ro, vec3 rd) 
{
    float dO = 0.; //Distane Origin
    for(int i=0;i<MAX_STEPS;i++)
    {
        vec3 p = ro + rd * dO;
        float ds = GetDist(p); // ds is Distance Scene
        dO += ds;
        if(dO > MAX_DIST || ds < SURFACE_DIST) break;
    }
    return dO;
}
 
vec3 GetNormal(vec3 p)
{
    float d = GetDist(p);
    vec2 e = vec2(.01,0);
    vec3 n = d - vec3(
    GetDist(p-e.xyy),
    GetDist(p-e.yxy),
    GetDist(p-e.yyx));
 
    return normalize(n);
}
 
float GetLight(vec3 p)
{ 
    // Light (directional diffuse)
    vec3 l = normalize(light_pos - p); // Light Vector
    vec3 n = GetNormal(p); // Normal Vector
   
    float dif = dot(n,l); // Diffuse light
    dif = clamp(dif,0.,1.); // Clamp so it doesnt go below 0
 
    return dif;
}

void main()
{
    u_resolution = vec2(1280., 720.);
    vec2 uv = (gl_FragCoord.xy - .5 * u_resolution.xy) / u_resolution.y;
    vec3 ro = vec3(0,1,0); // Ray Origin/Camera
    vec3 rd = normalize(vec3(uv.x,uv.y,1)); // Ray Direction
 
    float d = RayMarch(ro,rd); // Distance
   
    vec3 p = ro + rd * d;
    float dif = GetLight(p); // Diffuse lighting
    d*= .2;
    vec3 color = vec3(dif);
    //color += GetNormal(p);
    //float color = GetLight(p);
 
    // Set the output color
    frag_color = vec4(color,1.0);
}