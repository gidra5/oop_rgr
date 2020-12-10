extern crate piston_window;
extern crate vecmath;
extern crate camera_controllers;

#[macro_use]
extern crate gfx;
extern crate shader_version;
extern crate quaternion;
#[macro_use] extern crate conrod;

mod load_mesh;
// use load_mesh::*;

//----------------------------------------
// Cube associated data

gfx_defines!{
    vertex Vertex {
        pos: [i8; 4] = "a_pos",
        normal: [i8; 4] = "a_normal",
    }

    pipeline pipe {
        vbuf: gfx::VertexBuffer<Vertex> = (),
        u_resolution: gfx::Global<[f32; 2]> = "u_resolution",
        u_proj: gfx::Global<[[f32; 4]; 4]> = "u_proj",
        u_view: gfx::Global<[[f32; 4]; 4]> = "u_view",
        u_model: gfx::Global<[[f32; 4]; 4]> = "u_model",
        view_pos: gfx::Global<[f32; 3]> = "view_pos",
        light_pos: gfx::Global<[f32; 3]> = "light_pos",
        light_color: gfx::Global<[f32; 3]> = "light_color",
        out_color: gfx::RenderTarget<::gfx::format::Srgba8> = "frag_color",
        out_depth: gfx::DepthTarget<gfx::format::DepthStencil> =
            gfx::preset::depth::LESS_EQUAL_WRITE,
    }
}

impl Vertex {
    fn new(pos: [i8; 3], normal: [i8; 3]) -> Vertex {
        Vertex {
            pos: [pos[0], pos[1], pos[2], 1],
            normal: [normal[0], normal[1], normal[2], 1]
        }
    }
}

//----------------------------------------

fn main() {
    use piston_window::*;
    use gfx::traits::*;
    use camera_controllers::CameraPerspective;

    let opengl = OpenGL::V3_2;

    let mut window: PistonWindow =
        WindowSettings::new("piston: cube", [1280, 720])
        .exit_on_esc(true)
        .samples(4)
        .graphics_api(opengl)
        .build()
        .unwrap();
    // window.set_capture_cursor(true);

    let ref mut factory = window.factory.clone();

    let vertex_data = vec![
        Vertex::new([ 1,  1,  1], [ 1,  1,  1]),
        Vertex::new([ 1,  1, -1], [ 1,  1, -1]),
        Vertex::new([ 1, -1, -1], [ 1, -1, -1]),
        Vertex::new([ 1, -1,  1], [ 1, -1,  1]),
        Vertex::new([-1,  1,  1], [-1,  1,  1]),
        Vertex::new([-1,  1, -1], [-1,  1, -1]),
        Vertex::new([-1, -1, -1], [-1, -1, -1]),
        Vertex::new([-1, -1,  1], [-1, -1,  1]),
    ];

    let index_data: &[u16] = &[
        0, 1, 2, 0, 2, 3, //right
        6, 5, 4, 7, 6, 4, //left
        0, 1, 4, 1, 4, 5, //top
        0, 3, 7, 0, 7, 4, //front
        5, 2, 1, 6, 2, 5, //back
        2, 3, 7, 6, 7, 2, //bottom
    ];

    let (vbuf, slice) = factory.create_vertex_buffer_with_slice
        (&vertex_data, index_data);

    let pso = factory.create_pipeline_simple(
            include_bytes!("../assets/shading.glslv"),
            include_bytes!("../assets/shading.glslf"),
            pipe::new()
        ).unwrap();

    let get_projection = |w: &PistonWindow| {
        let draw_size = w.window.draw_size();
        CameraPerspective {
            fov: 60.0, near_clip: 0.1, far_clip: 1000.0,
            aspect_ratio: (draw_size.width as f32) / (draw_size.height as f32)
        }.projection()
    };

    let mut view_orientation = quaternion::id::<f32>();

    let mut data = pipe::Data {
        vbuf:           vbuf.clone(),
        u_resolution:   [1280., 720.],
        u_proj:         get_projection(&window),
        u_view:         vecmath::mat4_id(),
        u_model:        vecmath::mat4_id(),
        view_pos:       [0., 0., -8.],
        light_pos:      [50.; 3],
        light_color:    [1.; 3],
        out_color:      window.output_color.clone(),
        out_depth:      window.output_stencil.clone(),
    };

    let mut holding_mouse_button = None;

    while let Some(e) = window.next() {
        if holding_mouse_button != None {
            e.mouse_relative(|d| {
                let front = quaternion::rotate_vector(view_orientation, [0., 0., 1.]);
                let right = quaternion::rotate_vector(view_orientation, [1., 0., 0.]);
                let s = [front[2], 0., -front[0]]; //cross-product of front and y-axis

                let q_x = quaternion::axis_angle::<f32>([0., 1., 0.], -d[0] as f32 * 0.01);
                let q_y = quaternion::axis_angle::<f32>(vecmath::vec3_normalized(s), -d[1] as f32 * 0.01);
                let q_z = quaternion::rotation_from_to(right, [right[0], 0., right[2]]);

                view_orientation = quaternion::mul(q_x, view_orientation);
                view_orientation = quaternion::mul(q_y, view_orientation);
                view_orientation = quaternion::mul(q_z, view_orientation);
            });

            let (w, r) = view_orientation;

            let mut view = [
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 1.]
            ];

            let del = |i, j| if i == j { 1. } else { 0. };
            let eps = |i, j, k| -> f32 {
                ((i as i32 - j as i32) * (j as i32 - k as i32) * (k as i32 - i as i32)) as f32 / 2.
            };

            let mut cross_mat = [[0.; 3]; 3];

            (0..3).for_each(|m| (0..3).for_each(|k| {
                    cross_mat[m][k] = (0..3).map(|i|
                        (0..3).map(|j| del(m, i) * eps(i, j, k) * r[j]).sum::<f32>()
                    ).sum::<f32>();
                }
            ));

            (0..3).for_each(|i| (0..3).for_each(|j| {
                    view[j][i] = del(i, j) - 2. *
                        (w * cross_mat[i][j] - (0..3).map(|k| cross_mat[i][k] * cross_mat[k][j]).sum::<f32>());
                }
            ));

            data.u_view = view;
        }

        if let Some(Button::Mouse(button)) = e.press_args() {
            holding_mouse_button = Some(button);
        }

        if let Some(Button::Mouse(button)) = e.release_args() {
            if button == holding_mouse_button.unwrap() {
                holding_mouse_button = None;
            }
        };

        window.draw_3d(&e, |window| {
            let args = e.render_args().unwrap();

            window.encoder.clear(&window.output_color, [0.3, 0.3, 0.3, 1.0]);
            window.encoder.clear_depth(&window.output_stencil, 1.0);

            data.view_pos = quaternion::rotate_vector(view_orientation, [0., 0., -8.]);

            window.encoder.draw(&slice, &pso, &data);
        });

        if let Some(_) = e.resize_args() {
            data.u_proj = get_projection(&window);
            data.out_color = window.output_color.clone();
            data.out_depth = window.output_stencil.clone();
        }
    }
}
