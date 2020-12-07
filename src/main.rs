extern crate piston_window;
extern crate vecmath;
extern crate camera_controllers;

#[macro_use]
extern crate gfx;
extern crate shader_version;
extern crate dual_quaternion;
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
    use camera_controllers::{
        FirstPersonSettings,
        FirstPerson,
        CameraPerspective
    };

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

    let mut first_person = FirstPerson::new(
        [0.5, 0.5, 4.0],
        FirstPersonSettings{
            mouse_sensitivity_horizontal: 0.5,
            mouse_sensitivity_vertical: 0.5,
            ..FirstPersonSettings::keyboard_wasd()
        }
    );

    let view_orientation =
        dual_quaternion::from_rotation_and_translation((1., [0.;3]), [0.5, 0.5, 4.]);
    let view_mat = [[1., 0., 0.,  0. ],
                    [0., 1., 0.,  0. ],
                    [0., 0., 1.,  0. ],
                    [0., 0., 0.,  1. ]];

    let mut data = pipe::Data {
            vbuf:           vbuf.clone(),
            u_proj:         get_projection(&window),
            u_view:         view_mat,
            u_model:        vecmath::mat4_id(),
            view_pos:       [0.5, 0.5, 4.],
            light_pos:      [50.; 3],
            light_color:    [1.; 3],
            out_color:      window.output_color.clone(),
            out_depth:      window.output_stencil.clone(),
        };
    let mut holding_mouse_button = None;
    while let Some(e) = window.next() {
        let ((r_0, [r_1, r_2, r_3]), _) = view_orientation;
        let [x, y, z] = dual_quaternion::get_translation(view_orientation);
        let r = [r_1, r_2, r_3];

        // let mat1 = [
        //     [     r_0 * r_0 + r_1 * r_1, 2. * r_1 * r_2            , 2. * r_1 * r_3            ],
        //     [2. * r_2 * r_1            ,      r_0 * r_0 + r_2 * r_2, 2. * r_2 * r_3            ],
        //     [2. * r_3 * r_1            , 2. * r_3 * r_2            ,      r_0 * r_0 + r_3 * r_3]
        // ];
        // let mat2 = [
        //     [      r_2 * r_2 + r_3 * r_3, 2. * r_0 * r_3            , -2. * r_0 * r_2            ],
        //     [-2. * r_0 * r_3            ,      r_1 * r_1 + r_3 * r_3,  2. * r_0 * r_1            ],
        //     [ 2. * r_0 * r_2            ,-2. * r_0 * r_1            ,       r_1 * r_1 + r_2 * r_2]
        // ];
        let mut mat = [[0.; 3]; 3];
        let d = |i: f32, j: f32| if i == j { 1. } else { 0. };
        let ek = |j: f32, k: f32| -> f32 {
            (0..3).map(|x| ((x as f32 - j) * (j - k) * (k - x as f32)) * (d(j, k) - 1.) / 2.).sum()
        };

        (0..3).for_each(|i| (0..3).for_each(|j| {
                let mut t = ek(i as f32, j as f32);
                if t != 0. { t *= 2. * r[3 - i - j]; }

                mat[i][j] = r[i] * r[j] * (2. - d(i as f32, j as f32)) + r_0 * (r_0 * d(i as f32, j as f32) +
                t) + d(i as f32, j as f32) * (r[i] * r[i] - r[0] * r[0] - r[1] * r[1] - r[2] * r[2]);
            }
        ));

        let view = [
            [mat[0][0], mat[0][1], mat[0][2], x ],
            [mat[1][0], mat[1][1], mat[1][2], y ],
            [mat[2][0], mat[2][1], mat[2][2], z ],
            [0.       , 0.       , 0.       , 1.]
        ];
        data.u_view = view;
        if holding_mouse_button != None {

            first_person.event(&e);
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
// data.u_view = view_mat;
            // data.u_view = first_person.camera(args.ext_dt).orthogonal();
            data.view_pos = first_person.camera(args.ext_dt).position;

            window.encoder.draw(&slice, &pso, &data);
        });

        if let Some(_) = e.resize_args() {
            data.u_proj = get_projection(&window);
            data.out_color = window.output_color.clone();
            data.out_depth = window.output_stencil.clone();
        }
    }
}
