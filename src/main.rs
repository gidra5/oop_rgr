extern crate piston_window;
extern crate vecmath;
extern crate camera_controllers;

#[macro_use]
extern crate gfx;
extern crate shader_version;
extern crate quaternion;
// #[macro_use] extern crate conrod;

mod load_mesh;
// use load_mesh::*;

gfx_defines!{
    vertex Vertex2 {
        pos: [i16; 3] = "a_pos",
    }
    vertex Vertex3 {
        pos: [i8; 2] = "a_pos",
    }
    pipeline norm_pipe {
        vbuf:        gfx::VertexBuffer<  Vertex2                  > = (),
        u_proj:      gfx::Global      <  [[f32; 4]; 4]            > = "u_proj",
        u_view:      gfx::Global      <  [[f32; 4]; 4]            > = "u_view",
        out_color:   gfx::RenderTarget<::gfx::format::Rgba16F     > = "frag_color",
        out_depth:   gfx::DepthTarget <  gfx::format::DepthStencil> = gfx::preset::depth::LESS_EQUAL_WRITE,
    }

    pipeline pipe {
        vbuf:        gfx::VertexBuffer  <  Vertex3                  > = (),
        u_proj:      gfx::Global        <  [[f32; 4]; 4]            > = "u_proj",
        u_view:      gfx::Global        <  [[f32; 4]; 4]            > = "u_view",
        view_pos:    gfx::Global        <  [f32; 3]                 > = "view_pos",
        light_pos:   gfx::Global        <  [f32; 3]                 > = "light_pos",
        light_color: gfx::Global        <  [f32; 3]                 > = "light_color",
        depth_tex:   gfx::TextureSampler<   f32                     > = "depth_texture",
        // verticies:   gfx::TextureSampler<  [i8 ; 3]                 > = "vertex_buf",
        // indicies:    gfx::TextureSampler<   u16                     > = "index_buf",
        out_color:   gfx::RenderTarget  <::gfx::format::Srgba8      > = "frag_color",
    }
}

// impl gfx::pso::buffer::Structure<gfx::shade::ConstFormat> for Cube {
//     fn query(_: &str) -> Option<gfx::pso::buffer::Element<gfx::shade::ConstFormat>> {
//         None
//     }
// }

// impl gfx::pso::buffer::Structure<gfx::shade::ConstFormat> for CubeIndex {
//     fn query(_: &str) -> Option<gfx::pso::buffer::Element<gfx::shade::ConstFormat>> {
//         None
//     }
// }
//----------------------------------------

fn main() {
    use piston_window::*;
    use gfx::{
        traits::*,
        memory::{ Bind, Usage },
        format::*,
        texture::*,
        Slice,
    };
    use camera_controllers::CameraPerspective;

    const SENSATIVIY: f32 = -0.001;
    let mut w = 1280;
    let mut h = 720;      

    let opengl = OpenGL::V3_2;

    let mut window: PistonWindow =
        WindowSettings::new("piston: cube", [w, h])
        .exit_on_esc(true)
        .samples(4)
        .graphics_api(opengl)
        .build()
        .unwrap();

    let ref mut factory = window.factory.clone();

    let vertex_data = vec![
        [ 1,  1,  1],
        [ 1,  1, -1],
        [ 1, -1, -1],
        [ 1, -1,  1],
        [-1,  1,  1],
        [-1,  1, -1],
        [-1, -1, -1],
        [-1, -1,  1],
    ];

    let index_data: &[u16] = &[
        0, 1, 2, 0, 2, 3, //right
        6, 5, 4, 7, 6, 4, //left
        0, 1, 4, 1, 4, 5, //top
        0, 3, 7, 0, 7, 4, //front
        5, 2, 1, 6, 2, 5, //back
        2, 3, 7, 6, 7, 2, //bottom
    ];

    let (normals_buf, slice2) = factory.create_vertex_buffer_with_slice
        (vertex_data.iter().map(|x| Vertex2 { pos: *x }).collect::<Vec<_>>().as_slice(), index_data);

    // let norm_pso = factory.create_pipeline_simple(
    //     include_bytes!("../assets/pixel_normal.glslv"),
    //     include_bytes!("../assets/pixel_normal.glslf"),
    //         norm_pipe::new()
    //     ).unwrap();
    let pso = factory.create_pipeline_simple(
        include_bytes!("../assets/combine.glslv"),
        include_bytes!("../assets/combine.glslf"),
        pipe::new()
    ).unwrap();

    let proj = {
        CameraPerspective {
            fov: 60.0, near_clip: 0.1, far_clip: 10.0,
            aspect_ratio: w as f32 / h as f32
        }.projection()
    };

    let mut view_orientation = quaternion::id::<f32>();
    let mut scaling = 4.;
    let kind = gfx::texture::Kind::D2(w as u16, h as u16, AaMode::Single);
    let norm_tex = factory.create_texture(
        kind,
        1, 
        Bind::SHADER_RESOURCE | Bind::RENDER_TARGET, 
        Usage::Data, 
        Some(ChannelType::Float)
    ).unwrap();
    let depth_buf = factory.create_depth_stencil(w as u16, h as u16).unwrap();
    
    let mut data2 = norm_pipe::Data {
        vbuf:           normals_buf,
        u_proj:         proj,
        u_view:         [
            [1., 0.,       0., 0.],
            [0., 1.,       0., 0.],
            [0., 0.,       1., 0.],
            [0., 0., -scaling, 1.]
        ], 
        out_color:      factory.view_texture_as_render_target(&norm_tex, 0, None).unwrap(),
        out_depth:      depth_buf.2,
    };

    let mut data3 = pipe::Data {
        vbuf:           factory.create_vertex_buffer(&[
                Vertex3 { pos: [ 1,  1] },
                Vertex3 { pos: [-1,  1] },
                Vertex3 { pos: [ 1, -1] },
                Vertex3 { pos: [-1, -1] },
                Vertex3 { pos: [-1,  1] },
                Vertex3 { pos: [ 1, -1] }
            ]),
        u_proj:         proj,
        u_view:         [
            [1., 0.,       0., 0.],
            [0., 1.,       0., 0.],
            [0., 0.,       1., 0.],
            [0., 0., -scaling, 1.]
        ], 
        depth_tex:       (
            depth_buf.1, 
            // factory.view_texture_as_shader_resource::<DepthStencil>(&depth_buf.0, (0, 0), Swizzle::new()).unwrap(), 
            factory.create_sampler(SamplerInfo::new(FilterMethod::Bilinear, WrapMode::Clamp))
        ),
        // indicies:        (
        //     factory.create_texture_immutable::<gfx::format::R8>(
        //         gfx::texture::Kind::D2(2, 2, gfx::texture::AaMode::Single),
        //         gfx::texture::Mipmap::Provided,
        //         &[&index_data]).unwrap().1,
        //     factory.create_sampler(SamplerInfo::new(FilterMethod::Bilinear, WrapMode::Clamp))
        // ),
        // verticies:        (
        //     vert_buf.1, 
        //     // factory.view_texture_as_shader_resource::<DepthStencil>(&depth_buf.0, (0, 0), Swizzle::new()).unwrap(), 
        //     factory.create_sampler(SamplerInfo::new(FilterMethod::Bilinear, WrapMode::Clamp))
        // ),
        view_pos:       [0., 0., -scaling],
        light_pos:      [2., 2.,  8.],
        light_color:    [1.; 3],
        out_color:      window.output_color.clone(),
    };

    let mut holding_mouse_button = None;

    while let Some(e) = window.next() {
        if holding_mouse_button != None {
            e.mouse_relative(|d| {
                let right = quaternion::rotate_vector(view_orientation, [1., 0., 0.]);

                // TODO: fix bugs from using reversed 
                // or create other method to invert controls when uppside down

                //calculate if camera is upside down
                // let invert = (r[1] + (1. - r[1]) * (2. * w * w - 1.)).signum();  
                let q_x = quaternion::axis_angle::<f32>(
                    [0., 1., 0.], 
                    // d[0] as f32 * scaling * invert * SENSATIVIY
                    d[0] as f32 * scaling * SENSATIVIY
                );
                let q_y = quaternion::axis_angle::<f32>(
                    right, 
                    d[1] as f32 * scaling * SENSATIVIY
                );
                let q_z = quaternion::rotation_from_to(
                    right, 
                    [right[0], 0., right[2]]
                );

                view_orientation = quaternion::mul(q_x, view_orientation);
                view_orientation = quaternion::mul(q_y, view_orientation);
                view_orientation = quaternion::mul(q_z, view_orientation);
            });
            
            let mut view = [
                [0., 0.,       0., 0.],
                [0., 0.,       0., 0.],
                [0., 0.,       0., 0.],
                [0., 0., -scaling, 1.]
            ];
            let rot = quat_to_mat3(view_orientation);

            (0..3).for_each(|i| (0..3).for_each(|j| {
                view[j][i] = rot[i][j];
                }
            ));

            data2.u_view = view;
            data3.u_view = view;
        }

        e.mouse_scroll(|d| {
            scaling *= 0.95_f64.powf(d[1]) as f32;

            data2.u_view[3][2] = -scaling;
            data3.u_view[3][2] = -scaling;
        });

        if let Some(Button::Mouse(button)) = e.press_args() {
            holding_mouse_button = Some(button);
        }

        if let Some(Button::Mouse(button)) = e.release_args() {
            if button == holding_mouse_button.unwrap() {
                holding_mouse_button = None;
            }
        };

        window.draw_3d(&e, |window| {
            window.encoder.clear(&data2.out_color, [0., 0., 0., 1.]);
            window.encoder.clear_depth(&data2.out_depth, 1.0);

            data3.view_pos  = quaternion::rotate_vector(view_orientation, [0., 0., -scaling]);
            data3.light_pos = quaternion::rotate_vector(view_orientation, vecmath::vec3_scale([2., 2., 8.], scaling));

            // window.encoder.draw(&slice2, &norm_pso, &data2);
            window.encoder.draw(&Slice::new_match_vertex_buffer(&data3.vbuf), &pso, &data3);
        });

        // if let Some(_) = e.resize_args() {
        //     // data.u_proj = get_projection(&window);
            // data3.out_color = window.output_color.clone();
        //     data.out_depth = window.output_stencil.clone();
        // }
    }
}

fn quat_to_mat3((w, r): quaternion::Quaternion<f32>) -> vecmath::Matrix3<f32> {
    let mut mat = [[0.; 3]; 3];

    let del = |i, j| (i == j) as i32 as f32 ;
    let eps = |i, j, k| {
        ((i as i32 - j as i32) * 
         (j as i32 - k as i32) * 
         (k as i32 - i as i32)) as f32 / 2.
    };

    let mut cross_mat = [[0.; 3]; 3];

    (0..3).for_each(|m| (0..3).for_each(|k| {
            cross_mat[m][k] = (0..3).map(|i|
                (0..3).map(|j| del(m, i) * eps(i, j, k) * r[j]).sum::<f32>()
            ).sum::<f32>();
        }
    ));

    (0..3).for_each(|i| (0..3).for_each(|j| {
        mat[i][j] = del(i, j) - 2. *
            (w * cross_mat[i][j] - (0..3).map(|k| cross_mat[i][k] * cross_mat[k][j]).sum::<f32>());
        }
    ));
    
    mat
}
