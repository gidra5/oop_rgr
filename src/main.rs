extern crate camera_controllers;
extern crate image;
extern crate piston_window;
extern crate vecmath;

#[macro_use]
extern crate gfx;
extern crate quaternion;
extern crate shader_version;
#[macro_use]
extern crate conrod_core;
extern crate conrod_piston;

extern crate gfx_device_gl;

use self::piston_window::texture::UpdateTexture;
use std::path::Path;

gfx_defines! {
    vertex Vertex {
        pos: [i8; 2] = "a_pos",
    }

    pipeline pipe {
        vbuf:            gfx::VertexBuffer  <  Vertex                   > = (),
        u_proj:          gfx::Global        <  [[f32; 4]; 4]            > = "u_proj",
        u_view:          gfx::Global        <  [[f32; 4]; 4]            > = "u_view",
        light_pos:       gfx::Global        <  [f32; 3]                 > = "light_pos",
        light_color:     gfx::Global        <  [f32; 3]                 > = "light_color",
        u_res:           gfx::Global        <  [f32; 2]                 > = "u_resolution",
        sphere_center:   gfx::Global        <  [f32; 3]                 > = "sphere_center",
        plane_center:    gfx::Global        <  [f32; 3]                 > = "plane_center",
        cylinder_center: gfx::Global        <  [f32; 3]                 > = "cylinder_center",
        t:               gfx::Global        <  f32                      > = "t",
        samples:         gfx::Global        <  u32                      > = "samples",
        out_color:       gfx::RenderTarget  <::gfx::format::Srgba8      > = "frag_color",
        skybox:          gfx::TextureSampler<  [f32; 4]                 > = "skybox",
    }
    pipeline blur_pipe {
        vbuf:            gfx::VertexBuffer  <  Vertex                   > = (),
        t:               gfx::Global        <  f32                      > = "t",
        out_color:       gfx::RenderTarget  <::gfx::format::Srgba8      > = "frag_color",
        // prev:            gfx::TextureSampler< [f32; 3]                  > = "prev", // array texture of 8 prev frames
    }
}

widget_ids! {
    pub struct Ids {
        background,

        toggle_light_options,

        text_light_pos,
        xypad_light_pos_dx_dz,
        slider_light_pos_dy,

        text_light_color,
        text_samples,
        slider_samples,
        xypad_light_color_hue_brightness,
        slider_light_color_r,
        slider_light_color_g,
        slider_light_color_b,

        // toggle_change_quality,
        // dialer_change_quality_max_steps,
        // dialer_change_quality_max_dist,
        // dialer_change_quality_min_dist,
        // dialer_change_quality_light_soft_samples,
        toggle_shapes_position,

        text_plane_pos,
        xypad_plane_pos_dx_dz,
        slider_plane_pos_dy,

        text_sphere_pos,
        xypad_sphere_pos_dx_dz,
        slider_sphere_pos_dy,

        text_cylinder_pos,
        xypad_cylinder_pos_dx_dz,
        slider_cylinder_pos_dy,

        toggle_light_follow,
    }
}

//----------------------------------------

use conrod_core::*;
use piston_window::*;

use camera_controllers::CameraPerspective;
use gfx::{traits::*, Slice};

fn main() {
    const SENSATIVIY: f32 = -0.001;
    const SPEED: f32 = 10.;
    let mut view_orientation = quaternion::id::<f32>();
    let mut scaling = 1.;
    let mut holding_mouse_button = None;
    let mut mv = [0.; 6];
    let mut shape_menu = false;
    let mut light_menu = false;
    let mut light_follow = false;

    let opengl = OpenGL::V3_2;

    let mut window: PistonWindow = WindowSettings::new("piston: cube", [1280, 720])
        // let mut window: PistonWindow = WindowSettings::new("piston: cube", [1920, 1080])
        // let mut window: PistonWindow = WindowSettings::new("piston: cube", [2560, 1440])
        .exit_on_esc(true)
        .graphics_api(opengl)
        // .fullscreen(true)
        .build()
        .unwrap();

    let ref mut factory = window.factory.clone();

    let pso = factory
        .create_pipeline_simple(
            include_bytes!("../assets/combine.glslv"),
            include_bytes!("../assets/combine.glslf"),
            pipe::new(),
        )
        .unwrap();

    let mut data = setup(&mut window, scaling);

    let mut texture_context = window.create_texture_context();
    let (mut ui, mut glyph_cache, mut text_texture_cache, mut text_vertex_data) = {
        let piston_window::Size { width, height } = window.window.draw_size();
        // construct our `Ui`.
        let mut ui = conrod_core::UiBuilder::new([width, height])
            // .theme(conrod_example_shared::theme())
            .build();
        // Add a `Font` to the `Ui`'s `font::Map` from file.

        ui.fonts
            .insert_from_file("D:\\Projects\\oop-rgr\\assets\\FiraSans-Regular.ttf")
            .unwrap();
        // Create a texture to use for efficiently caching text on the GPU.
        let text_vertex_data = Vec::new();
        let (glyph_cache, text_texture_cache) = {
            const SCALE_TOLERANCE: f32 = 0.1;
            const POSITION_TOLERANCE: f32 = 0.1;
            let cache = conrod_core::text::GlyphCache::builder()
                .dimensions(width as u32, height as u32)
                .scale_tolerance(SCALE_TOLERANCE)
                .position_tolerance(POSITION_TOLERANCE)
                .build();
            let buffer_len = width as usize * width as usize;
            let init = vec![128; buffer_len];
            let settings = TextureSettings::new();
            let texture = G2dTexture::from_memory_alpha(
                &mut texture_context,
                &init,
                width as u32,
                height as u32,
                &settings,
            )
            .unwrap();
            (cache, texture)
        };
        (ui, glyph_cache, text_texture_cache, text_vertex_data)
    };
    let ids = Ids::new(ui.widget_id_generator());
    let image_map = conrod_core::image::Map::new();

    while let Some(e) = window.next() {
        let right = quaternion::rotate_vector(view_orientation, [1., 0., 0.]);
        let mv_up = [0., 1., 0.];
        let mv_right = [right[0], 0., right[2]];
        let mv_front = {
            let front = quaternion::rotate_vector(view_orientation, [0., 0., 1.]);
            [front[0], 0., front[2]]
        };

        if holding_mouse_button != None {
            e.mouse_relative(|d| {
                let q_x = quaternion::axis_angle::<f32>(
                    [0., 1., 0.],
                    -d[0] as f32 / scaling * SENSATIVIY,
                );
                let q_y = quaternion::axis_angle::<f32>(right, -d[1] as f32 / scaling * SENSATIVIY);
                let q_z = quaternion::rotation_from_to(right, mv_right);

                view_orientation = quaternion::mul(q_x, view_orientation);
                view_orientation = quaternion::mul(q_y, view_orientation);
                view_orientation = quaternion::mul(q_z, view_orientation);
            });

            let rot = quat_to_mat3(view_orientation);

            (0..3).for_each(|i| {
                (0..3).for_each(|j| {
                    data.u_view[i][j] = rot[i][j];
                })
            });
        }

        e.mouse_scroll(|d| {
            scaling *= 0.95_f64.powf(-d[1]) as f32;

            data.u_view[3][3] = scaling;
        });

        match e.press_args() {
            Some(Button::Mouse(button)) => {
                holding_mouse_button = Some(button);
            }
            Some(Button::Keyboard(key)) => match key {
                Key::W => mv[0] = SPEED,
                Key::A => mv[1] = -SPEED,
                Key::S => mv[3] = SPEED,
                Key::D => mv[4] = -SPEED,
                Key::Space => mv[5] = SPEED,
                Key::LShift => mv[2] = SPEED,
                _ => {}
            },
            _ => {}
        }

        match e.release_args() {
            Some(Button::Mouse(button)) => {
                if button == holding_mouse_button.unwrap() {
                    holding_mouse_button = None;
                }
            }
            Some(Button::Keyboard(key)) => match key {
                Key::W => mv[0] = 0.,
                Key::A => mv[1] = 0.,
                Key::S => mv[3] = 0.,
                Key::D => mv[4] = 0.,
                Key::Space => mv[5] = 0.,
                Key::LShift => mv[2] = 0.,
                _ => {}
            },
            _ => {}
        }

        window.draw_3d(&e, |window| {
            window
                .encoder
                .draw(&Slice::new_match_vertex_buffer(&data.vbuf), &pso, &data);
            // window.encoder.draw(&slice2, &norm_pso, &data2);
            // window.encoder.draw(&Slice::new_match_vertex_buffer(&data3.vbuf), &pso, &data3);
        });

        // Convert the src event to a conrod event.
        let size = window.size();
        if let Some(e) = conrod_piston::event::convert(e.clone(), size.width, size.height) {
            ui.handle_event(e);
        }

        e.update(|args| {
            data.t += args.dt as f32;
            let k = 5. * args.dt as f32;
            let mut dir = vecmath::vec3_sub([mv[0], mv[1], mv[2]], [mv[3], mv[4], mv[5]]);
            dir = vecmath::col_mat3_transform([mv_front, mv_right, mv_up], dir);
            dir = vecmath::vec3_scale(dir, args.dt as f32);

            data.u_view[3] = vecmath::vec4_add(data.u_view[3], [dir[0], -dir[1], dir[2], 0.]);

            if light_follow {
                data.light_pos = vecmath::vec3_mul(
                    vecmath::vec3_add(
                        [data.u_view[3][0], data.u_view[3][1], data.u_view[3][2]],
                        quaternion::rotate_vector(view_orientation, [0., 0., -2.]),
                    ),
                    [1., -1., 1.],
                );
            }

            let mut ui = ui.set_widgets();

            {
                const MARGIN: conrod_core::Scalar = 30.0;
                let mut bg_height = 2. * MARGIN + 200.;
                bg_height += if light_menu {
                    2. * (20. + 20. + 20. + 50.)
                } else {
                    0.
                };
                bg_height += if shape_menu {
                    3. * (20. + 20. + 20. + 50.)
                } else {
                    0.
                };

                widget::bordered_rectangle::BorderedRectangle::new([200., bg_height])
                    .top_left_with_margin_on(ui.window, MARGIN)
                    .set(ids.background, &mut ui);

                widget::Text::new("Samples")
                    .mid_top_with_margin_on(ids.background, MARGIN)
                    .down(20.)
                    .set(ids.text_samples, &mut ui);
                for samples in widget::NumberDialer::new(data.samples as f32, 1., 256., 1)
                    .w_h(50., 10.)
                    .x_relative_to(ids.background, 30.)
                    .down(20.)
                    .set(ids.slider_samples, &mut ui)
                {
                    data.samples = samples as u32;
                }

                for toggled in widget::Toggle::new(light_menu)
                    .color(Color::Rgba(0.5, 0.5, 0.5, 0.5))
                    .label("Light options")
                    .mid_top_with_margin_on(ids.background, MARGIN)
                    .w_h(160., 30.)
                    .set(ids.toggle_light_options, &mut ui)
                {
                    light_menu = toggled;
                }
                if light_menu {
                    widget::Text::new("Color")
                        .align_middle_x_of(ids.toggle_light_options)
                        .down(20.)
                        .set(ids.text_light_pos, &mut ui);

                    for g in widget::Slider::new(data.light_color[1], 0., 1.)
                        .w_h(10., 50.)
                        .align_middle_x_of(ids.toggle_light_options)
                        .down(20.)
                        .set(ids.slider_light_color_g, &mut ui)
                    {
                        data.light_color[1] = g;
                    }
                    for b in widget::Slider::new(data.light_color[2], 0., 1.)
                        .w_h(10., 50.)
                        .right(10.)
                        .set(ids.slider_light_color_b, &mut ui)
                    {
                        data.light_color[2] = b;
                    }
                    for r in widget::Slider::new(data.light_color[0], 0., 1.)
                        .w_h(10., 50.)
                        .left_from(ids.slider_light_color_g, 10.)
                        .set(ids.slider_light_color_r, &mut ui)
                    {
                        data.light_color[0] = r;
                    }

                    widget::Text::new("Position")
                        .down(20.)
                        .align_middle_x_of(ids.text_light_pos)
                        .set(ids.text_light_color, &mut ui);
                    for dy in widget::Slider::new(0., -1., 1.)
                        .w_h(10., 50.)
                        .x_relative_to(ids.background, 30.)
                        .down(20.)
                        .set(ids.slider_light_pos_dy, &mut ui)
                    {
                        data.light_pos[1] += k * dy;
                    }
                    for (dx, dz) in widget::XYPad::new(0., -1., 1., 0., -1., 1.)
                        .w_h(50., 50.)
                        .left(10.)
                        .set(ids.xypad_light_pos_dx_dz, &mut ui)
                    {
                        data.light_pos[2] += k * dz;
                        data.light_pos[0] -= k * dx;
                    }
                }

                for toggled in widget::Toggle::new(shape_menu)
                    .color(Color::Rgba(0.5, 0.5, 0.5, 0.5))
                    .label("Shape positioning")
                    .w_h(160., 30.)
                    .down(20.)
                    .align_left_of(ids.toggle_light_options)
                    .set(ids.toggle_shapes_position, &mut ui)
                {
                    shape_menu = toggled;
                }

                if shape_menu {
                    widget::Text::new("Plane position")
                        .align_middle_x_of(ids.toggle_light_options)
                        .down(20.)
                        .set(ids.text_plane_pos, &mut ui);

                    for dy in widget::Slider::new(0., -1., 1.)
                        .w_h(10., 50.)
                        .x_relative_to(ids.background, 30.)
                        .down(20.)
                        .set(ids.slider_plane_pos_dy, &mut ui)
                    {
                        data.plane_center[1] += k * dy;
                    }
                    for (dx, dz) in widget::XYPad::new(0., -1., 1., 0., -1., 1.)
                        .w_h(50., 50.)
                        .left(10.)
                        .set(ids.xypad_plane_pos_dx_dz, &mut ui)
                    {
                        data.plane_center[2] += k * dz;
                        data.plane_center[0] -= k * dx;
                    }
                    widget::Text::new("Sphere position")
                        .align_middle_x_of(ids.toggle_light_options)
                        .down(20.)
                        .set(ids.text_sphere_pos, &mut ui);

                    for dy in widget::Slider::new(0., -1., 1.)
                        .w_h(10., 50.)
                        .x_relative_to(ids.background, 30.)
                        .down(20.)
                        .set(ids.slider_sphere_pos_dy, &mut ui)
                    {
                        data.sphere_center[1] += k * dy;
                    }
                    for (dx, dz) in widget::XYPad::new(0., -1., 1., 0., -1., 1.)
                        .w_h(50., 50.)
                        .left(10.)
                        .set(ids.xypad_sphere_pos_dx_dz, &mut ui)
                    {
                        data.sphere_center[2] += k * dz;
                        data.sphere_center[0] -= k * dx;
                    }

                    widget::Text::new("Cylinder position")
                        .align_middle_x_of(ids.toggle_light_options)
                        .down(20.)
                        .set(ids.text_cylinder_pos, &mut ui);

                    for dy in widget::Slider::new(0., -1., 1.)
                        .w_h(10., 50.)
                        .x_relative_to(ids.background, 30.)
                        .down(20.)
                        .set(ids.slider_cylinder_pos_dy, &mut ui)
                    {
                        data.cylinder_center[1] += k * dy;
                    }
                    for (dx, dz) in widget::XYPad::new(0., -1., 1., 0., -1., 1.)
                        .w_h(50., 50.)
                        .left(10.)
                        .set(ids.xypad_cylinder_pos_dx_dz, &mut ui)
                    {
                        data.cylinder_center[2] += k * dz;
                        data.cylinder_center[0] -= k * dx;
                    }
                }
                for toggled in widget::Toggle::new(light_follow)
                    .color(Color::Rgba(0.5, 0.5, 0.5, 0.5))
                    .label("Light follows you")
                    .w_h(160., 30.)
                    .down(20.)
                    .align_left_of(ids.toggle_light_options)
                    .set(ids.toggle_light_follow, &mut ui)
                {
                    light_follow = toggled;
                }
            }
        });

        window.draw_2d(&e, |context, graphics, device| {
            // A function used for caching glyphs to the texture cache.
            let cache_queued_glyphs = |_graphics: &mut G2d,
                                       cache: &mut G2dTexture,
                                       rect: conrod_core::text::rt::Rect<u32>,
                                       data: &[u8]| {
                let offset = [rect.min.x, rect.min.y];
                let size = [rect.width(), rect.height()];
                let format = piston_window::texture::Format::Rgba8;
                text_vertex_data.clear();
                text_vertex_data.extend(data.iter().flat_map(|&b| vec![255, 255, 255, b]));
                UpdateTexture::update(
                    cache,
                    &mut texture_context,
                    format,
                    &text_vertex_data[..],
                    offset,
                    size,
                )
                .expect("failed to update texture")
            };

            // Specify how to get the drawable texture from the image. In this case, the image
            // *is* the texture.
            fn texture_from_image<T>(img: &T) -> &T {
                img
            }

            // Draw the conrod `render::Primitives`.
            conrod_piston::draw::primitives(
                ui.draw(),
                context,
                graphics,
                &mut text_texture_cache,
                &mut glyph_cache,
                &image_map,
                cache_queued_glyphs,
                texture_from_image,
            );

            texture_context.encoder.flush(device);
        });

        if let Some(_) = e.resize_args() {
            resize(&window, &mut data);
        }
    }
}

fn quat_to_mat3((w, r): quaternion::Quaternion<f32>) -> vecmath::Matrix3<f32> {
    let mut mat = [[0.; 3]; 3];

    let del = |i, j| (i == j) as i32 as f32;
    let eps = |i, j, k| {
        ((i as i32 - j as i32) * (j as i32 - k as i32) * (k as i32 - i as i32)) as f32 / 2.
    };

    let mut cross_mat = [[0.; 3]; 3];

    (0..3).for_each(|m| {
        (0..3).for_each(|k| {
            cross_mat[m][k] = (0..3)
                .map(|i| (0..3).map(|j| del(m, i) * eps(i, j, k) * r[j]).sum::<f32>())
                .sum::<f32>();
        })
    });

    (0..3).for_each(|i| {
        (0..3).for_each(|j| {
            mat[i][j] = del(i, j)
                - 2. * (w * cross_mat[i][j]
                    - (0..3)
                        .map(|k| cross_mat[i][k] * cross_mat[k][j])
                        .sum::<f32>());
        })
    });

    mat
}

fn resize(window: &piston_window::PistonWindow, data: &mut pipe::Data<gfx_device_gl::Resources>) {
    let piston_window::Size { width, height } = window.window.draw_size();

    data.u_proj = {
        CameraPerspective {
            fov: 60.0,
            near_clip: 0.1,
            far_clip: 10.0,
            aspect_ratio: (width / height) as f32,
        }
        .projection()
    };
    data.u_res = [width as f32, height as f32];
    data.out_color = window.output_color.clone();
}

fn setup(
    window: &mut piston_window::PistonWindow,
    scaling: f32,
) -> pipe::Data<gfx_device_gl::Resources> {
    let piston_window::Size { width, height } = window.window.draw_size();

    let ref mut factory = window.factory.clone();
    let skybox = factory
        .create_texture_immutable_u8::< (gfx::format::R8_G8_B8_A8, gfx::format::Unorm)>(
            // gfx::texture::Kind::D2(2048, 2048, gfx::texture::AaMode::Single),
            gfx::texture::Kind::Cube(2048),
            gfx::texture::Mipmap::Provided,
            &[
                &image::load(std::io::Cursor::new(include_bytes!("../assets/right.jpg")), image::ImageFormat::Jpeg).unwrap().to_rgba8(),
                &image::load(std::io::Cursor::new(include_bytes!("../assets/left.jpg")), image::ImageFormat::Jpeg).unwrap().to_rgba8(),
                &image::load(std::io::Cursor::new(include_bytes!("../assets/top.jpg")), image::ImageFormat::Jpeg).unwrap().to_rgba8(),
                &image::load(std::io::Cursor::new(include_bytes!("../assets/bottom.jpg")), image::ImageFormat::Jpeg).unwrap().to_rgba8(),
                &image::load(std::io::Cursor::new(include_bytes!("../assets/back.jpg")), image::ImageFormat::Jpeg).unwrap().to_rgba8(),
                &image::load(std::io::Cursor::new(include_bytes!("../assets/front.jpg")), image::ImageFormat::Jpeg).unwrap().to_rgba8(),
            ],
        )
        .unwrap();
    // let prev = factory
    //     .create_texture::<gfx::format::R32_G32_B32>(
    //         gfx::texture::Kind::D2Array(12, 8, 0, gfx::texture::AaMode::Single),
    //         1,
    //         gfx::memory::Bind::SHADER_RESOURCE,
    //         gfx::memory::Usage::Data,
    //         None,
    //     )
    //     .unwrap();
    let data = pipe::Data {
        vbuf: factory.create_vertex_buffer(&[
            Vertex { pos: [1, 1] },
            Vertex { pos: [-1, 1] },
            Vertex { pos: [1, -1] },
            Vertex { pos: [-1, -1] },
            Vertex { pos: [-1, 1] },
            Vertex { pos: [1, -1] },
        ]),
        u_proj: {
            CameraPerspective {
                fov: 60.0,
                near_clip: 0.1,
                far_clip: 10.0,
                aspect_ratio: (width / height) as f32,
            }
            .projection()
        },
        u_view: [
            [1., 0., 0., 0.],
            [0., 1., 0., 0.],
            [0., 0., 1., 0.],
            [0., 0., -6., scaling],
        ],
        u_res: [width as f32, height as f32],

        light_pos: [2., 0., 1.],
        light_color: [0xfc as f32 / 255., 0x0f as f32 / 255., 0xc0 as f32 / 255.],
        sphere_center: [-1., 0., 0.],
        plane_center: [0., -1., 0.],
        cylinder_center: [1., 0., 4.],

        out_color: window.output_color.clone(),
        skybox: (
            skybox.1,
            factory.create_sampler(gfx::texture::SamplerInfo::new(
                gfx::texture::FilterMethod::Bilinear,
                gfx::texture::WrapMode::Clamp,
            )),
        ),
        t: 0. as f32,
        samples: 4,
    };
    let blur_data = blur_pipe::Data {
        vbuf: factory.create_vertex_buffer(&[
            Vertex { pos: [1, 1] },
            Vertex { pos: [-1, 1] },
            Vertex { pos: [1, -1] },
            Vertex { pos: [-1, -1] },
            Vertex { pos: [-1, 1] },
            Vertex { pos: [1, -1] },
        ]),

        out_color: window.output_color.clone(),
        // prev: (
        //     factory
        //         .view_texture_as_shader_resource::<[f32; 3]>(
        //             &prev,
        //             (0, 1),
        //             gfx::format::Swizzle::new(),
        //         )
        //         .unwrap(),
        //     factory.create_sampler(gfx::texture::SamplerInfo::new(
        //         gfx::texture::FilterMethod::Bilinear,
        //         gfx::texture::WrapMode::Clamp,
        //     )),
        // ),
        t: 0. as f32,
    };

    data
}
