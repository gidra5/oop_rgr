extern crate piston_window;
extern crate vecmath;
extern crate camera_controllers;

#[macro_use]
extern crate gfx;
extern crate shader_version;
#[macro_use] extern crate conrod;

mod load_mesh;
use load_mesh::*;

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

    let mut data = pipe::Data {
            vbuf:           vbuf.clone(),
            u_proj:         get_projection(&window),
            u_view:         vecmath::mat4_id(),
            u_model:        vecmath::mat4_id(),
            view_pos:       [0.; 3],
            light_pos:      [50.; 3],
            light_color:    [1.; 3],
            out_color:      window.output_color.clone(),
            out_depth:      window.output_stencil.clone(),
        };

    while let Some(e) = window.next() {
        if let Some(Button::Mouse(button)) = e.press_args() {
            first_person.event(&e);
            println!("Pressed mouse button '{:?}'", button);
        }
        if let Some(Button::Keyboard(key)) = e.press_args() {
            if key == Key::C {
                println!("Turned capture cursor on");
                capture_cursor = !capture_cursor;
                window.set_capture_cursor(capture_cursor);
            }

            println!("Pressed keyboard key '{:?}'", key);
        };
        if let Some(args) = e.button_args() {
            println!("Scancode {:?}", args.scancode);
        }
        if let Some(button) = e.release_args() {
            match button {
                Button::Keyboard(key) => println!("Released keyboard key '{:?}'", key),
                Button::Mouse(button) => println!("Released mouse button '{:?}'", button),
                Button::Controller(button) => println!("Released controller button '{:?}'", button),
                Button::Hat(hat) => println!("Released controller hat `{:?}`", hat),
            }
        };

        window.draw_3d(&e, |window| {
            let args = e.render_args().unwrap();

            window.encoder.clear(&window.output_color, [0.3, 0.3, 0.3, 1.0]);
            window.encoder.clear_depth(&window.output_stencil, 1.0);

            data.u_view = first_person.camera(args.ext_dt).orthogonal();
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
