extern crate piston_window;
extern crate vecmath;
extern crate camera_controllers;

#[macro_use]
extern crate gfx;
extern crate shader_version;
extern crate quaternion;
// #[macro_use] extern crate conrod;
extern crate conrod_core;
extern crate conrod_piston;

extern crate gfx_device_gl;

// use self::piston_window::texture::UpdateTexture;
// use self::piston_window::OpenGL;
// use self::piston_window::{Flip, G2d, G2dTexture, Texture, TextureSettings};

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
        light_pos:   gfx::Global        <  [f32; 3]                 > = "light_pos",
        light_color: gfx::Global        <  [f32; 3]                 > = "light_color",
        norm_tex:    gfx::TextureSampler<  [f32; 4]                 > = "normal_texture",
        depth_tex:   gfx::TextureSampler<  f32                      > = "depth_texture",
        // verticies:   gfx::TextureSampler<  [i8 ; 3]                 > = "vertex_buf",
        // indicies:    gfx::TextureSampler<   u16                     > = "index_buf",
        out_color:   gfx::RenderTarget  <::gfx::format::Srgba8      > = "frag_color",
    }
}

//----------------------------------------

use piston_window::*;
use gfx::{
    traits::*,
    memory::{ Bind, Usage },
    format::*,
    texture::*,
    Slice,
};
use camera_controllers::CameraPerspective;

fn main() {
    const SENSATIVIY: f32 = -0.001; 
    const SPEED: f32 = 0.05;

    let opengl = OpenGL::V3_2;

    let mut window: PistonWindow =
        WindowSettings::new("piston: cube", [1280, 720])
        .exit_on_esc(true)
        .samples(4)
        .graphics_api(opengl)
        .build()
        .unwrap();

    let ref mut factory = window.factory.clone();

    let norm_pso = factory.create_pipeline_simple(
        include_bytes!("../assets/pixel_normal.glslv"),
        include_bytes!("../assets/pixel_normal.glslf"),
            norm_pipe::new()
        ).unwrap();
    let pso = factory.create_pipeline_simple(
        include_bytes!("../assets/combine.glslv"),
        include_bytes!("../assets/combine.glslf"),
        pipe::new()
    ).unwrap();

    let mut view_orientation = quaternion::id::<f32>();
    let mut scaling = 4.;

    let (mut data2, mut data3, slice2) = setup(&window, scaling);

    let mut holding_mouse_button = None;
    let mut mv = [0.; 4]; 

    while let Some(e) = window.next() {
        let front = quaternion::rotate_vector(view_orientation, [0., 0., 1.]);
        let right = quaternion::rotate_vector(view_orientation, [1., 0., 0.]);
        
        if holding_mouse_button != None {
            e.mouse_relative(|d| {
                let q_x = quaternion::axis_angle::<f32>(
                    [0., 1., 0.], 
                    d[0] as f32 * SENSATIVIY
                );
                let q_y = quaternion::axis_angle::<f32>(
                    right, 
                    d[1] as f32 * SENSATIVIY
                );
                let q_z = quaternion::rotation_from_to(
                    right, 
                    [right[0], 0., right[2]]
                );

                view_orientation = quaternion::mul(q_x, view_orientation);
                view_orientation = quaternion::mul(q_y, view_orientation);
                view_orientation = quaternion::mul(q_z, view_orientation);
            });
            
            let rot = quat_to_mat3(view_orientation);

            (0..3).for_each(|i| (0..3).for_each(|j| {
                data2.u_view[j][i] = rot[i][j];
                data3.u_view[j][i] = rot[i][j];
                }
            ));
        }

        e.mouse_scroll(|d| {
            scaling *= 0.95_f64.powf(d[1]) as f32;

            data2.u_view[3][3] = scaling;
            data3.u_view[3][3] = scaling;
        });

        match e.press_args() {
            Some(Button::Mouse(button)) => {
                holding_mouse_button = Some(button);
            },
            Some(Button::Keyboard(key)) => {
                match key {
                    Key::W => mv[0] = SPEED,
                    Key::A => mv[1] = SPEED,
                    Key::S => mv[2] = SPEED,
                    Key::D => mv[3] = SPEED,
                    _ => {}
                }
            },
            _ => {}
        }

        match e.release_args() {
            Some(Button::Mouse(button)) => {
                if button == holding_mouse_button.unwrap() {
                    holding_mouse_button = None;
                }
            },
            Some(Button::Keyboard(key)) => {
                match key {
                    Key::W => mv[0] = 0.,
                    Key::A => mv[1] = 0.,
                    Key::S => mv[2] = 0.,
                    Key::D => mv[3] = 0.,
                    _ => {}
                }
            },
            _ => {}
        }

        window.draw_3d(&e, |window| {
            (0..3).for_each(|i| {
                data2.u_view[3][i] += (mv[0] - mv[2]) * front[i] + (mv[1] - mv[3]) * right[i];
                data3.u_view[3][i] += (mv[0] - mv[2]) * front[i] + (mv[1] - mv[3]) * right[i];
            });

            window.encoder.clear(&data2.out_color, [1., 1., 1., 1.]);
            window.encoder.clear_depth(&data2.out_depth, 1.0);

            window.encoder.draw(&slice2, &norm_pso, &data2);
            window.encoder.draw(&Slice::new_match_vertex_buffer(&data3.vbuf), &pso, &data3);
        });

        if let Some(_) = e.resize_args() {
            resize(&window, &mut data2, &mut data3);
        }
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

fn resize(window: &piston_window::PistonWindow, 
    data2: &mut norm_pipe::Data<gfx_device_gl::Resources>, 
    data3: &mut pipe::Data<gfx_device_gl::Resources>) 
{
    let piston_window::Size {
        width,
        height
    } = window.window.draw_size();

    let ref mut factory = window.factory.clone();

    let norm_tex = factory.create_texture(
        gfx::texture::Kind::D2(width as u16, height as u16, AaMode::Single),
        1, 
        Bind::SHADER_RESOURCE | Bind::RENDER_TARGET, 
        Usage::Data, 
        Some(ChannelType::Float)
    ).unwrap();
    let (_, depth_buf, depth_view) = factory.create_depth_stencil(width as u16, height as u16).unwrap();
    
    data2.u_proj = {
        CameraPerspective {
            fov: 60.0, near_clip: 0.1, far_clip: 10.0,
            aspect_ratio: (width / height) as f32
        }.projection()
    };
    data2.out_color = factory.view_texture_as_render_target(&norm_tex, 0, None).unwrap();
    data2.out_depth = depth_view;

    data3.u_proj = {
        CameraPerspective {
            fov: 60.0, near_clip: 0.1, far_clip: 10.0,
            aspect_ratio: (width / height) as f32
        }.projection()
    };
    data3.depth_tex = (
        depth_buf, 
        factory.create_sampler(SamplerInfo::new(FilterMethod::Bilinear, WrapMode::Clamp))
    );
    data3.out_color = window.output_color.clone();
}

fn setup(
    window: &piston_window::PistonWindow,
    scaling: f32
) 
-> (norm_pipe::Data<gfx_device_gl::Resources>, 
    pipe::Data<gfx_device_gl::Resources>, 
    gfx::Slice<gfx_device_gl::Resources>
) {
    let piston_window::Size {
        width,
        height
    } = window.window.draw_size();

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

    let norm_tex = factory.create_texture(
        gfx::texture::Kind::D2(width as u16, height as u16, AaMode::Single),
        1, 
        Bind::SHADER_RESOURCE | Bind::RENDER_TARGET, 
        Usage::Data, 
        Some(ChannelType::Float)
    ).unwrap();

    let (_, depth_buf, depth_view) = factory.create_depth_stencil(width as u16, height as u16).unwrap();
    
    let data2 = norm_pipe::Data {
        vbuf:           normals_buf,
        u_proj:         {
            CameraPerspective {
                fov: 60.0, near_clip: 0.1, far_clip: 10.0,
                aspect_ratio: (width / height) as f32
            }.projection()
        },
        u_view:         [
            [1., 0.,  0.,      0.],
            [0., 1.,  0.,      0.],
            [0., 0.,  1.,      0.],
            [0., 0., -6., scaling]
        ], 
        out_color:      factory.view_texture_as_render_target(&norm_tex, 0, None).unwrap(),
        out_depth:      depth_view,
    };

    let data3 = pipe::Data {
        vbuf:           factory.create_vertex_buffer(&[
                Vertex3 { pos: [ 1,  1] },
                Vertex3 { pos: [-1,  1] },
                Vertex3 { pos: [ 1, -1] },
                Vertex3 { pos: [-1, -1] },
                Vertex3 { pos: [-1,  1] },
                Vertex3 { pos: [ 1, -1] }
            ]),
        u_proj:         {
            CameraPerspective {
                fov: 60.0, near_clip: 0.1, far_clip: 10.0,
                aspect_ratio: (width / height) as f32
            }.projection()
        },
        u_view:         [
            [1., 0.,  0.,      0.],
            [0., 1.,  0.,      0.],
            [0., 0.,  1.,      0.],
            [0., 0., -6., scaling]
        ], 
        norm_tex:       (
            factory.view_texture_as_shader_resource::<Rgba16F>(&norm_tex, (0, 0), Swizzle::new()).unwrap(), 
            factory.create_sampler(SamplerInfo::new(FilterMethod::Bilinear, WrapMode::Clamp))
        ),
        depth_tex: (
            depth_buf, 
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
        light_pos:      [2., 8.,  8.],
        light_color:    [0xfc as f32 / 255., 0x0f as f32 / 255., 0xc0 as f32 / 255.],
        out_color:      window.output_color.clone(),
    };

    (data2, data3, slice2)
}