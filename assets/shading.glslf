#version 150 core

in float diffusion;
in vec3 frag_pos;
in vec3 normal;

out vec4 frag_color;
uniform vec3 view_pos;
uniform vec3 light_pos;
uniform vec3 light_color;

void main() {
    vec3 obj_color = vec3(0.5, 0., 0.);

    vec3 norm = normalize(normal);
    vec3 light_dir = normalize(light_pos - frag_pos);
    vec3 view_dir = normalize(view_pos - frag_pos);
    vec3 reflect_dir = reflect(-light_dir, norm);

    // float diffusion = max(dot(norm, light_dir), 0.0);
    float ambience = 0.1;

    float specular_strength = 0.1;
    float specularity = specular_strength * pow(max(dot(view_dir, reflect_dir), 0.0), 3);

    frag_color = vec4(
        (ambience + diffusion + specularity) * light_color * obj_color,
        1.
    );
}
