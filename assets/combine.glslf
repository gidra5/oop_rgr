#version 150 core
precision highp float;

#define PHI 1.61803398874989484820459 // Golden Ratio   
#define SRT 1.41421356237309504880169 // Square Root of Two
#define PI 3.14159265358979323846264
#define TWO_PI 6.28318530717958647692528

uniform mat4 u_proj;
uniform mat4 u_view;
uniform vec2 u_resolution;
uniform samplerCube skybox;
uniform float t;
uniform uint aliasing_samples;
uniform uint lens_samples;
uniform uint light_samples;
uniform uint reflection_samples;
uniform uint gi_reflection_depth;

const   float max_dist  = 1000000.;
const   float min_dist  = 0.0005;
// const   int   max_steps = 1000;

uniform float cameraFovAngle;
uniform float paniniDistance;
// uniform float verticalCompression;
// const   float cameraFovAngle = PI * 2. / 3.;
// const   float paniniDistance = 0.75;
// const   float verticalCompression = 0.1;

// const   float imagePlaneDistance = 1.;
// const   float lensFocalLength = 1.1;
// const   float circleOfConfusionRadius = 0.01;
uniform float imagePlaneDistance;
uniform float lensFocalLength;
uniform float circleOfConfusionRadius;

const mat3 triangle_pts = mat3(
  vec3(0.),
  vec3(0., 1., 0.),
  vec3(1., 1., 0.)
);
const vec3 sun_color = vec3(0x92, 0x97, 0xC4) / 0xff * 0.9;
const float ambience = 0.05;

uniform vec3 light_pos;
uniform vec3 light_color;
uniform vec3 sphere_center;
uniform vec3 plane_center;
uniform vec3 cylinder_center;

struct Ray {
  vec3 pos; // Origin
  vec3 dir; // Direction (normalized)
};
struct Hit {
  vec3 normal;
  vec3 color;
};

in vec3 ray_dir;

out vec4 frag_color;






float random_0t1(in vec2 coordinate, in float seed) {
  int base = 1<<9;
  int modulo = 1<<10;
  float seed_mod = fract(seed / modulo) * modulo + base;
  return fract(sin(dot(coordinate * seed_mod, vec2(PHI, PI)) * .1) * SRT * 10000.0);
  // return fract(sin(dot(coordinate * seed, vec2(PHI, PI)) * .1) * SRT * 10000.0);
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






float dot2( in vec3 v ) { return dot(v,v); }

// Plane 
float iPlane( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
              in vec3 planeNormal, in float planeDist) {
    float a = dot(rd, planeNormal);
    float d = -(dot(ro, planeNormal)+planeDist)/a;
    if (a > 0. || d < distBound.x || d > distBound.y) {
        return max_dist;
    } else {
        normal = planeNormal;
    	return d;
    }
}

// Sphere:          https://www.shadertoy.com/view/4d2XWV
float iSphere( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
               float sphereRadius ) {
    float b = dot(ro, rd);
    float c = dot(ro, ro) - sphereRadius*sphereRadius;
    float h = b*b - c;
    if (h < 0.) {
        return max_dist;
    } else {
	    h = sqrt(h);
        float d1 = -b-h;
        float d2 = -b+h;
        if (d1 >= distBound.x && d1 <= distBound.y) {
            normal = normalize(ro + rd*d1);
            return d1;
        } else if (d2 >= distBound.x && d2 <= distBound.y) { 
            normal = normalize(ro + rd*d2);            
            return d2;
        } else {
            return max_dist;
        }
    }
}

// Box:             https://www.shadertoy.com/view/ld23DV
float iBox( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal, 
            in vec3 boxSize ) {
    vec3 m = sign(rd)/max(abs(rd), 1e-8);
    vec3 n = m*ro;
    vec3 k = abs(m)*boxSize;
	
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;

	float tN = max( max( t1.x, t1.y ), t1.z );
	float tF = min( min( t2.x, t2.y ), t2.z );
	
    if (tN > tF || tF <= 0.) {
        return max_dist;
    } else {
        if (tN >= distBound.x && tN <= distBound.y) {
        	normal = -sign(rd)*step(t1.yzx,t1.xyz)*step(t1.zxy,t1.xyz);
            return tN;
        } else if (tF >= distBound.x && tF <= distBound.y) { 
        	normal = -sign(rd)*step(t1.yzx,t1.xyz)*step(t1.zxy,t1.xyz);
            return tF;
        } else {
            return max_dist;
        }
    }
}

// Capped Cylinder: https://www.shadertoy.com/view/4lcSRn
float iCylinder( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
                 in vec3 pa, in vec3 pb, float ra ) {
    vec3 ca = pb-pa;
    vec3 oc = ro-pa;

    float caca = dot(ca,ca);
    float card = dot(ca,rd);
    float caoc = dot(ca,oc);
    
    float a = caca - card*card;
    float b = caca*dot( oc, rd) - caoc*card;
    float c = caca*dot( oc, oc) - caoc*caoc - ra*ra*caca;
    float h = b*b - a*c;
    
    if (h < 0.) return max_dist;
    
    h = sqrt(h);
    float d = (-b-h)/a;

    float y = caoc + d*card;
    if (y > 0. && y < caca && d >= distBound.x && d <= distBound.y) {
        normal = (oc+d*rd-ca*y/caca)/ra;
        return d;
    }

    d = ((y < 0. ? 0. : caca) - caoc)/card;
    
    if( abs(b+a*d) < h && d >= distBound.x && d <= distBound.y) {
        normal = normalize(ca*sign(y)/caca);
        return d;
    } else {
        return max_dist;
    }
}

// Torus:           https://www.shadertoy.com/view/4sBGDy
float iTorus( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
              in vec2 torus ) {
    // bounding sphere
    vec3 tmpnormal;
    if (iSphere(ro, rd, distBound, tmpnormal, torus.y+torus.x) > distBound.y) {
        return max_dist;
    }
    
    float po = 1.0;
    
	float Ra2 = torus.x*torus.x;
	float ra2 = torus.y*torus.y;
	
	float m = dot(ro,ro);
	float n = dot(ro,rd);

#if 1
	float k = (m + Ra2 - ra2)/2.0;
  float k3 = n;
	float k2 = n*n - Ra2*dot(rd.xy,rd.xy) + k;
  float k1 = n*k - Ra2*dot(rd.xy,ro.xy);
  float k0 = k*k - Ra2*dot(ro.xy,ro.xy);
#else
	float k = (m - Ra2 - ra2)/2.0;
	float k3 = n;
	float k2 = n*n + Ra2*rd.z*rd.z + k;
	float k1 = k*n + Ra2*ro.z*rd.z;
	float k0 = k*k + Ra2*ro.z*ro.z - Ra2*ra2;
#endif
    
#if 1
  // prevent |c1| from being too close to zero
  if (abs(k3*(k3*k3-k2)+k1) < 0.01) {
      po = -1.0;
      float tmp=k1; k1=k3; k3=tmp;
      k0 = 1.0/k0;
      k1 = k1*k0;
      k2 = k2*k0;
      k3 = k3*k0;
  }
#endif
    
    // reduced cubic
    float c2 = k2*2.0 - 3.0*k3*k3;
    float c1 = k3*(k3*k3-k2)+k1;
    float c0 = k3*(k3*(c2+2.0*k2)-8.0*k1)+4.0*k0;
    
    c2 /= 3.0;
    c1 *= 2.0;
    c0 /= 3.0;

    float Q = c2*c2 + c0;
    float R = c2*c2*c2 - 3.0*c2*c0 + c1*c1;
    
    float h = R*R - Q*Q*Q;
    float t = max_dist;
    
    if (h>=0.0) {
        // 2 intersections
        h = sqrt(h);
        
        float v = sign(R+h)*pow(abs(R+h),1.0/3.0); // cube root
        float u = sign(R-h)*pow(abs(R-h),1.0/3.0); // cube root

        vec2 s = vec2( (v+u)+4.0*c2, (v-u)*sqrt(3.0));
    
        float y = sqrt(0.5*(length(s)+s.x));
        float x = 0.5*s.y/y;
        float r = 2.0*c1/(x*x+y*y);

        float t1 =  x - r - k3; t1 = (po<0.0)?2.0/t1:t1;
        float t2 = -x - r - k3; t2 = (po<0.0)?2.0/t2:t2;

        if (t1 >= distBound.x) t=t1;
        if (t2 >= distBound.x) t=min(t,t2);
	} else {
        // 4 intersections
        float sQ = sqrt(Q);
        float w = sQ*cos( acos(-R/(sQ*Q)) / 3.0 );

        float d2 = -(w+c2); if( d2<0.0 ) return max_dist;
        float d1 = sqrt(d2);

        float h1 = sqrt(w - 2.0*c2 + c1/d1);
        float h2 = sqrt(w - 2.0*c2 - c1/d1);
        float t1 = -d1 - h1 - k3; t1 = (po<0.0)?2.0/t1:t1;
        float t2 = -d1 + h1 - k3; t2 = (po<0.0)?2.0/t2:t2;
        float t3 =  d1 - h2 - k3; t3 = (po<0.0)?2.0/t3:t3;
        float t4 =  d1 + h2 - k3; t4 = (po<0.0)?2.0/t4:t4;

        if (t1 >= distBound.x) t=t1;
        if (t2 >= distBound.x) t=min(t,t2);
        if (t3 >= distBound.x) t=min(t,t3);
        if (t4 >= distBound.x) t=min(t,t4);
    }
    
	if (t >= distBound.x && t <= distBound.y) {
        vec3 pos = ro + rd*t;
        normal = normalize( pos*(dot(pos,pos) - torus.y*torus.y - torus.x*torus.x*vec3(1,1,-1)));
        return t;
    } else {
        return max_dist;
    }
}

// Capsule:         https://www.shadertoy.com/view/Xt3SzX
float iCapsule( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
                in vec3 pa, in vec3 pb, in float r ) {
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;

    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoa = dot(ba,oa);
    float rdoa = dot(rd,oa);
    float oaoa = dot(oa,oa);

    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - r*r*baba;
    float h = b*b - a*c;
    if (h >= 0.) {
        float t = (-b-sqrt(h))/a;
        float d = max_dist;
        
        float y = baoa + t*bard;
        
        // body
        if (y > 0. && y < baba) {
            d = t;
        } else {
            // caps
            vec3 oc = (y <= 0.) ? oa : ro - pb;
            b = dot(rd,oc);
            c = dot(oc,oc) - r*r;
            h = b*b - c;
            if( h>0.0 ) {
                d = -b - sqrt(h);
            }
        }
        if (d >= distBound.x && d <= distBound.y) {
            vec3  pa = ro + rd * d - pa;
            float h = clamp(dot(pa,ba)/dot(ba,ba),0.0,1.0);
            normal = (pa - h*ba)/r;
            return d;
        }
    }
    return max_dist;
}

// Capped Cone:     https://www.shadertoy.com/view/llcfRf
float iCone( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
             in vec3  pa, in vec3  pb, in float ra, in float rb ) {
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    vec3  ob = ro - pb;
    
    float m0 = dot(ba,ba);
    float m1 = dot(oa,ba);
    float m2 = dot(ob,ba); 
    float m3 = dot(rd,ba);

    //caps
    if (m1 < 0.) { 
        if( dot2(oa*m3-rd*m1)<(ra*ra*m3*m3) ) {
            float d = -m1/m3;
            if (d >= distBound.x && d <= distBound.y) {
                normal = -ba*inversesqrt(m0);
                return d;
            }
        }
    }
    else if (m2 > 0.) { 
        if( dot2(ob*m3-rd*m2)<(rb*rb*m3*m3) ) {
            float d = -m2/m3;
            if (d >= distBound.x && d <= distBound.y) {
                normal = ba*inversesqrt(m0);
                return d;
            }
        }
    }
                       
    // body
    float m4 = dot(rd,oa);
    float m5 = dot(oa,oa);
    float rr = ra - rb;
    float hy = m0 + rr*rr;
    
    float k2 = m0*m0    - m3*m3*hy;
    float k1 = m0*m0*m4 - m1*m3*hy + m0*ra*(rr*m3*1.0        );
    float k0 = m0*m0*m5 - m1*m1*hy + m0*ra*(rr*m1*2.0 - m0*ra);
    
    float h = k1*k1 - k2*k0;
    if( h < 0. ) return max_dist;

    float t = (-k1-sqrt(h))/k2;

    float y = m1 + t*m3;
    if (y > 0. && y < m0 && t >= distBound.x && t <= distBound.y) {
        normal = normalize(m0*(m0*(oa+t*rd)+rr*ba*ra)-ba*hy*y);
        return t;
    } else {   
	    return max_dist;
    }
}

// Ellipsoid:       https://www.shadertoy.com/view/MlsSzn
float iEllipsoid( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
                  in vec3 rad ) {
    vec3 ocn = ro / rad;
    vec3 rdn = rd / rad;
    
    float a = dot( rdn, rdn );
	float b = dot( ocn, rdn );
	float c = dot( ocn, ocn );
	float h = b*b - a*(c-1.);
    
    if (h < 0.) {
        return max_dist;
    }
    
	float d = (-b - sqrt(h))/a;
    
    if (d < distBound.x || d > distBound.y) {
        return max_dist;
    } else {
        normal = normalize((ro + d*rd)/rad);
    	return d;
    }
}

// Rounded Cone:    https://www.shadertoy.com/view/MlKfzm
float iRoundedCone( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
                    in vec3  pa, in vec3  pb, in float ra, in float rb ) {
    vec3  ba = pb - pa;
	vec3  oa = ro - pa;
	vec3  ob = ro - pb;
    float rr = ra - rb;
    float m0 = dot(ba,ba);
    float m1 = dot(ba,oa);
    float m2 = dot(ba,rd);
    float m3 = dot(rd,oa);
    float m5 = dot(oa,oa);
	float m6 = dot(ob,rd);
    float m7 = dot(ob,ob);
    
    float d2 = m0-rr*rr;
    
	float k2 = d2    - m2*m2;
    float k1 = d2*m3 - m1*m2 + m2*rr*ra;
    float k0 = d2*m5 - m1*m1 + m1*rr*ra*2. - m0*ra*ra;
    
	float h = k1*k1 - k0*k2;
    if (h < 0.0) {
        return max_dist;
    }
    
    float t = (-sqrt(h)-k1)/k2;
    
    float y = m1 - ra*rr + t*m2;
    if (y>0.0 && y<d2) {
        if (t >= distBound.x && t <= distBound.y) {
        	normal = normalize( d2*(oa + t*rd)-ba*y );
            return t;
        } else {
            return max_dist;
        }
    } else {
        float h1 = m3*m3 - m5 + ra*ra;
        float h2 = m6*m6 - m7 + rb*rb;

        if (max(h1,h2)<0.0) {
            return max_dist;
        }

        vec3 n = vec3(0);
        float r = max_dist;

        if (h1 > 0.) {        
            r = -m3 - sqrt( h1 );
            n = (oa+r*rd)/ra;
        }
        if (h2 > 0.) {
            t = -m6 - sqrt( h2 );
            if( t<r ) {
                n = (ob+t*rd)/rb;
                r = t;
            }
        }
        if (r >= distBound.x && r <= distBound.y) {
            normal = n;
            return r;
        } else {
            return max_dist;
        }
    }
}

// Triangle:        https://www.shadertoy.com/view/MlGcDz
float iTriangle( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
                 in vec3 v0, in vec3 v1, in vec3 v2 ) {
    vec3 v1v0 = v1 - v0;
    vec3 v2v0 = v2 - v0;
    vec3 rov0 = ro - v0;

    vec3  n = cross( v1v0, v2v0 );
    vec3  q = cross( rov0, rd );
    float d = 1.0/dot( rd, n );
    float u = d*dot( -q, v2v0 );
    float v = d*dot(  q, v1v0 );
    float t = d*dot( -n, rov0 );

    if( u<0. || v<0. || (u+v)>1. || t<distBound.x || t>distBound.y) {
        return max_dist;
    } else {
        normal = normalize(-n);
        return t;
    }
}

// Sphere4:         https://www.shadertoy.com/view/3tj3DW
float iSphere4( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
                in float ra ) {
    // -----------------------------
    // solve quartic equation
    // -----------------------------
    
    float r2 = ra*ra;
    
    vec3 d2 = rd*rd; vec3 d3 = d2*rd;
    vec3 o2 = ro*ro; vec3 o3 = o2*ro;

    float ka = 1.0/dot(d2,d2);

    float k0 = ka* dot(ro,d3);
    float k1 = ka* dot(o2,d2);
    float k2 = ka* dot(o3,rd);
    float k3 = ka*(dot(o2,o2) - r2*r2);

    // -----------------------------
    // solve cubic
    // -----------------------------

    float c0 = k1 - k0*k0;
    float c1 = k2 + 2.0*k0*(k0*k0 - (3.0/2.0)*k1);
    float c2 = k3 - 3.0*k0*(k0*(k0*k0 - 2.0*k1) + (4.0/3.0)*k2);

    float p = c0*c0*3.0 + c2;
    float q = c0*c0*c0 - c0*c2 + c1*c1;
    float h = q*q - p*p*p*(1.0/27.0);

    // -----------------------------
    // skip the case of 3 real solutions for the cubic, which involves 
    // 4 complex solutions for the quartic, since we know this objcet is 
    // convex
    // -----------------------------
    if (h<0.0) {
        return max_dist;
    }
    
    // one real solution, two complex (conjugated)
    h = sqrt(h);

    float s = sign(q+h)*pow(abs(q+h),1.0/3.0); // cuberoot
    float t = sign(q-h)*pow(abs(q-h),1.0/3.0); // cuberoot

    vec2 v = vec2( (s+t)+c0*4.0, (s-t)*sqrt(3.0) )*0.5;
    
    // -----------------------------
    // the quartic will have two real solutions and two complex solutions.
    // we only want the real ones
    // -----------------------------
    
    float r = length(v);
	float d = -abs(v.y)/sqrt(r+v.x) - c1/r - k0;

    if (d >= distBound.x && d <= distBound.y) {
	    vec3 pos = ro + rd * d;
	    normal = normalize( pos*pos*pos );
	    return d;
    } else {
        return max_dist;
    }
}

// Goursat:         https://www.shadertoy.com/view/3lj3DW
float cuberoot( float x ) { return sign(x)*pow(abs(x),1.0/3.0); }

float iGoursat( in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
                in float ra, float rb ) {
// hole: x4 + y4 + z4 - (r2^2)·(x2 + y2 + z2) + r1^4 = 0;
    float ra2 = ra*ra;
    float rb2 = rb*rb;
    
    vec3 rd2 = rd*rd; vec3 rd3 = rd2*rd;
    vec3 ro2 = ro*ro; vec3 ro3 = ro2*ro;

    float ka = 1.0/dot(rd2,rd2);

    float k3 = ka*(dot(ro ,rd3));
    float k2 = ka*(dot(ro2,rd2) - rb2/6.0);
    float k1 = ka*(dot(ro3,rd ) - rb2*dot(rd,ro)/2.0  );
    float k0 = ka*(dot(ro2,ro2) + ra2*ra2 - rb2*dot(ro,ro) );

    float c2 = k2 - k3*(k3);
    float c1 = k1 + k3*(2.0*k3*k3-3.0*k2);
    float c0 = k0 + k3*(k3*(c2+k2)*3.0-4.0*k1);

    c0 /= 3.0;

    float Q = c2*c2 + c0;
    float R = c2*c2*c2 - 3.0*c0*c2 + c1*c1;
    float h = R*R - Q*Q*Q;
    
    
    // 2 intersections
    if (h>0.0) {
        h = sqrt(h);

        float s = cuberoot( R + h );
        float u = cuberoot( R - h );
        
        float x = s+u+4.0*c2;
        float y = s-u;
        
        float k2 = x*x + y*y*3.0;
  
        float k = sqrt(k2);

		float d = -0.5*abs(y)*sqrt(6.0/(k+x)) 
                  -2.0*c1*(k+x)/(k2+x*k) 
                  -k3;
        
        if (d >= distBound.x && d <= distBound.y) {
            vec3 pos = ro + rd * d;
            normal = normalize( 4.0*pos*pos*pos - 2.0*pos*rb*rb );
            return d;
        } else {
            return max_dist;
        }
    } else {	
        // 4 intersections
        float sQ = sqrt(Q);
        float z = c2 - 2.0*sQ*cos( acos(-R/(sQ*Q)) / 3.0 );

        float d1 = z   - 3.0*c2;
        float d2 = z*z - 3.0*c0;

        if (abs(d1)<1.0e-4) {  
            if( d2<0.0) return max_dist;
            d2 = sqrt(d2);
        } else {
            if (d1<0.0) return max_dist;
            d1 = sqrt( d1/2.0 );
            d2 = c1/d1;
        }

        //----------------------------------

        float h1 = sqrt(d1*d1 - z + d2);
        float h2 = sqrt(d1*d1 - z - d2);
        float t1 = -d1 - h1 - k3;
        float t2 = -d1 + h1 - k3;
        float t3 =  d1 - h2 - k3;
        float t4 =  d1 + h2 - k3;

        if (t2<0.0 && t4<0.0) return max_dist;

        float result = 1e20;
             if (t1>0.0) result=t1;
        else if (t2>0.0) result=t2;
             if (t3>0.0) result=min(result,t3);
        else if (t4>0.0) result=min(result,t4);

        if (result >= distBound.x && result <= distBound.y) {
            vec3 pos = ro + rd * result;
            normal = normalize( 4.0*pos*pos*pos - 2.0*pos*rb*rb );
            return result;
        } else {
            return max_dist;
        }
    }
}

// Rounded Box:     https://www.shadertoy.com/view/WlSXRW
float iRoundedBox(in vec3 ro, in vec3 rd, in vec2 distBound, inout vec3 normal,
   				  in vec3 size, in float rad ) {
	// bounding box
    vec3 m = 1.0/rd;
    vec3 n = m*ro;
    vec3 k = abs(m)*(size+rad);
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;
	float tN = max( max( t1.x, t1.y ), t1.z );
	float tF = min( min( t2.x, t2.y ), t2.z );
    if (tN > tF || tF < 0.0) {
    	return max_dist;
    }
    float t = (tN>=distBound.x&&tN<=distBound.y)?tN:
    		  (tF>=distBound.x&&tF<=distBound.y)?tF:max_dist;

    // convert to first octant
    vec3 pos = ro+t*rd;
    vec3 s = sign(pos);
    vec3 ros = ro*s;
    vec3 rds = rd*s;
    pos *= s;
        
    // faces
    pos -= size;
    pos = max( pos.xyz, pos.yzx );
    if (min(min(pos.x,pos.y),pos.z)<0.0) {
        if (t >= distBound.x && t <= distBound.y) {
            vec3 p = ro + rd * t;
            normal = sign(p)*normalize(max(abs(p)-size,0.0));
            return t;
        }
    }
    
    // some precomputation
    vec3 oc = ros - size;
    vec3 dd = rds*rds;
	vec3 oo = oc*oc;
    vec3 od = oc*rds;
    float ra2 = rad*rad;

    t = max_dist;        

    // corner
    {
    float b = od.x + od.y + od.z;
	float c = oo.x + oo.y + oo.z - ra2;
	float h = b*b - c;
	if (h > 0.0) t = -b-sqrt(h);
    }

    // edge X
    {
	float a = dd.y + dd.z;
	float b = od.y + od.z;
	float c = oo.y + oo.z - ra2;
	float h = b*b - a*c;
	if (h>0.0) {
	  h = (-b-sqrt(h))/a;
      if (h>=distBound.x && h<t && abs(ros.x+rds.x*h)<size.x ) t = h;
    }
	}
    // edge Y
    {
	float a = dd.z + dd.x;
	float b = od.z + od.x;
	float c = oo.z + oo.x - ra2;
	float h = b*b - a*c;
	if (h>0.0) {
	  h = (-b-sqrt(h))/a;
      if (h>=distBound.x && h<t && abs(ros.y+rds.y*h)<size.y) t = h;
    }
	}
    // edge Z
    {
	float a = dd.x + dd.y;
	float b = od.x + od.y;
	float c = oo.x + oo.y - ra2;
	float h = b*b - a*c;
	if (h>0.0) {
	  h = (-b-sqrt(h))/a;
      if (h>=distBound.x && h<t && abs(ros.z+rds.z*h)<size.z) t = h;
    }
	}
    
	if (t >= distBound.x && t <= distBound.y) {
        vec3 p = ro + rd * t;
        normal = sign(p)*normalize(max(abs(p)-size,1e-16));
        return t;
    } else {
        return max_dist;
    };
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
    // d = min(d, triangle(p - vec3(2., 0., 2.), triangle_pts[0], triangle_pts[1], triangle_pts[2]));

  return d;
}

vec3 dist_scene_gradient(vec3 p) {
  float d = raw_ds(p);
  vec2 e = vec2(min_dist, 0);

  return (vec3(raw_ds(p + e.xyy), raw_ds(p + e.yxy), raw_ds(p + e.yyx)) - d) / min_dist;
}

// vec3 raymarch(Ray ray, out bool hit) {
//   vec3 origin = ray.pos;
//   vec3 dir = ray.dir;
//   vec3 p = origin;
//   hit = false;

//   float d = 0.;
//   for (int steps = 0; steps < max_steps && d < max_dist && !hit; ++steps) {
//     float ds = dist_scene(Ray(p, dir));
//     // float ds = raw_ds(p);

//     d += ds;
//     p += dir * ds;

//     hit = ds < min_dist;
//   }

//   return p;
// }

float scene(Ray ray, out bool hit, out Hit hitObj) {
  float d = max_dist;
  float d2;

  d2 = iPlane(ray.pos - plane_center, ray.dir, vec2(0., d), hitObj.normal, vec3(0., 1., 0.), 0.);
  hit = hit || d2 < d;
  if (d2 < d) hitObj.color = vec3(1., 1., 1.);
  d = min(d, d2);

  d2 = iSphere(ray.pos - sphere_center, ray.dir, vec2(0., d), hitObj.normal, 1.);
  hit = hit || d2 < d;
  if (d2 < d) hitObj.color = vec3(0., 1., 0.);
  d = min(d, d2);

  d2 = iCylinder(ray.pos - sphere_center, ray.dir, vec2(0., d), hitObj.normal, vec3(0., 1., 1.), vec3(0., 2., 2.), 1.);
  hit = hit || d2 < d;
  if (d2 < d) hitObj.color = vec3(1., 1., 1.);
  d = min(d, d2);

  d2 = iPlane(ray.pos - plane_center - vec3(-3., 0., 0.), ray.dir, vec2(0., d), hitObj.normal, vec3(1., 0., 0.), 0.);
  hit = hit || d2 < d;
  if (d2 < d) hitObj.color = vec3(1., 1., 1.);
  d = min(d, d2);

  d2 = iBox(ray.pos - plane_center - vec3(0., 0., 5.), ray.dir, vec2(0., d), hitObj.normal, vec3(20., 9., 0.1));
  hit = hit || d2 < d;
  if (d2 < d) hitObj.color = vec3(1., 1., 0.);
  d = min(d, d2);

  d2 = iBox(ray.pos - plane_center - vec3(0., 0., -7.), ray.dir, vec2(0., d), hitObj.normal, vec3(20., 9., 0.1));
  hit = hit || d2 < d;
  if (d2 < d) hitObj.color = vec3(0., 1., 1.);
  d = min(d, d2);

  return d;
}






vec2 sample_circle(float t) {
  return vec2(cos(t * TWO_PI), sin(t * TWO_PI));
}
vec2 sample_incircle(vec2 t) {
  float theta = t.x * TWO_PI;
  float radius = sqrt(t.y);
  return vec2(cos(theta), sin(theta)) * radius;
}
vec3 sample_sphere(vec2 uv) {
  float sinTheta = sqrt(1 - uv.x * uv.x); 
  float phi = TWO_PI * uv.y; 
  float x = sinTheta * cos(phi); 
  float z = sinTheta * sin(phi); 
  return vec3(x, uv.x, z); 
}

float in_shadow(Ray ray, float mag_sq) {
  bool hit;
  Hit hitObj;
  float ds = scene(ray, hit, hitObj);

  return !hit || ds * ds >= mag_sq ? 1. : 0.;
}

const float k = 1./(2.*TWO_PI);

vec3 light(vec3 pos, vec3 norm, vec3 light_pos, vec3 light_color) {
    vec3 d = light_pos - pos;
    float mag_sq = dot(d, d);
    float mag = sqrt(mag_sq);
    vec3 dir = d / mag;
    float vis = in_shadow(Ray(pos + min_dist * norm, dir), mag_sq);

    return vis * max(dot(dir, norm), 0.) * light_color / mag_sq;
}
vec3 sun(vec3 pos, vec3 norm) {
    vec3 dir = normalize(vec3(1., 1., -1.));
    float vis = in_shadow(Ray(pos + min_dist * norm, dir), 0.99 / (min_dist * min_dist));

    return (vis * max(dot(dir, norm), 0.) + ambience) * sun_color;
}






vec3 pinholeRay(vec2 pixel) { 
  return vec3(pixel, 1/tan(cameraFovAngle / 2.f));
}

vec3 paniniRay(vec2 pixel) {
  float halfFOV = cameraFovAngle / 2.f;
  vec2 p = vec2(sin(halfFOV), cos(halfFOV) + paniniDistance);
  float M = sqrt(dot(p, p));
  float halfPaniniFOV = atan(p.x, p.y);
  vec2 hvPan = pixel * vec2(halfPaniniFOV, halfFOV);
  float x = sin(hvPan.x) * M;
  float z = cos(hvPan.x) * M - paniniDistance;
  // float y = tan(hvPan.y) * (z + verticalCompression);
  float y = tan(hvPan.y) * (z + pow(max(0., (3. * cameraFovAngle/PI - 1.) / 8.), 0.92));

  return vec3(x, y, z);
}

Ray thinLensRay(vec3 ray, vec2 uv) {
  float focusPlane = (imagePlaneDistance * lensFocalLength) / (imagePlaneDistance - lensFocalLength);
  vec3 focusPoint = ray * (focusPlane / ray.z);

  vec3 origin = vec3(uv * circleOfConfusionRadius, 0.f);
  vec3 direction = -normalize(focusPoint + origin);
  return Ray(origin, direction);
}






void main() {
  vec3 _frag_color = vec3(0.);
  for (int x = 0; x < int(aliasing_samples); ++x) {
    vec2 subpixel = random_0t1_2(gl_FragCoord.xy * t, t * x + 0.5);
    vec2 uv = (2. * (gl_FragCoord.xy + subpixel) - u_resolution) / u_resolution.x;
    vec3 rayDirection = normalize(paniniRay(uv));
    // vec3 rayDirection = normalize(pinholeRay(uv));
    vec3 subpixel_color = vec3(0.);

    for (int x = 0; x < int(lens_samples); ++x) {
      Ray ray = thinLensRay(rayDirection, sample_incircle(random_0t1_2(uv * t, t * x)));
      ray.dir = (u_view * vec4(normalize(ray.dir), 0.)).xyz;
      ray.dir = (u_view * vec4(rayDirection, 0.)).xyz;
      vec4 ray_pos = u_view * vec4(ray.pos, 1.);
      ray.pos = ray_pos.xyz / ray_pos.w;

      vec3 color = vec3(0.);

      for (int x = 0; x < int(light_samples); ++x) {
        vec2 p = sample_circle(random_0t1(uv * t, t * x));
        vec3 _light_pos = light_pos + vec3(p.x, 0., p.y);

        vec3 refl_color = vec3(0.);
        for (int x = 0; x < int(reflection_samples); ++x) {
          Ray rays[10];
          Hit hitObjs[9];
          vec3 indirect_color = vec3(0.);

          // // generate ray path
          // bool hit = true;
          // int i = 0;
          // rays[0] = ray;

          // for (; i < int(gi_reflection_depth) + 2 && hit; ++i) {
          //   hit = false;
          //   float _t = scene(rays[i], hit, hitObjs[i]);
          //   vec3 d = sample_sphere(random_0t1_2(uv * t, t * x));
          //   vec3 pos = rays[i].pos + _t * rays[i].dir + min_dist * hitObjs[i].normal;
          //   if (dot(d, hitObjs[i].normal) < 0.) d = -d;
          //   rays[i + 1] = Ray(pos, d);
          // }

          // // brackpropagate color info
          // for (int j = i - 1; j >= 1; --j) {
          //   vec3 light = 
          //     dot(hitObjs[j - 1].normal, rays[j].dir) * indirect_color +
          //     light(rays[j].pos, hitObjs[j - 1].normal, _light_pos, light_color) +
          //     sun(rays[j].pos, hitObjs[j - 1].normal);
          //   indirect_color = hitObjs[j - 1].color * light;  
          // }

          // generate ray path
          bool hit = false;
          float t = scene(ray, hit, hitObjs[0]);
          rays[0] = Ray(ray.pos + t * ray.dir + min_dist * hitObjs[0].normal, vec3(0.));
          int i = 1;

          for (; i < int(gi_reflection_depth) + 1 && hit; ++i) {
            vec3 d = sample_sphere(random_0t1_2(uv * t, t * x));
            if (dot(d, hitObjs[i - 1].normal) < 0.) d = -d;
            rays[i - 1].dir = d;

            float t = scene(rays[i - 1], hit, hitObjs[i]);
            rays[i] = Ray(rays[i - 1].pos + t * d + min_dist * hitObjs[i].normal, vec3(0.));
          }

          // brackpropagate color info
          for (int j = i - 1; j >= 0; --j) {
            vec3 light = 
              dot(hitObjs[j].normal, rays[j].dir) * indirect_color +
              light(rays[j].pos, hitObjs[j].normal, _light_pos, light_color) +
              sun(rays[j].pos, hitObjs[j].normal);
            indirect_color = hitObjs[j].color * light;  
          }
          refl_color += indirect_color;
        }

        color += refl_color / reflection_samples;
      }
      subpixel_color += color / light_samples; 
    }
    _frag_color += subpixel_color / lens_samples;
  }
  frag_color = vec4(_frag_color / aliasing_samples, 1.);
}