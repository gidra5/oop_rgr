extern crate piston_window;
extern crate vecmath;
extern crate camera_controllers;

#[macro_use]
extern crate gfx;
extern crate shader_version;
extern crate quaternion;
#[macro_use]
extern crate conrod_core;
extern crate conrod_piston;

extern crate gfx_device_gl;

use self::piston_window::texture::UpdateTexture;

gfx_defines!{
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
        out_color:       gfx::RenderTarget  <::gfx::format::Srgba8      > = "frag_color",
    }
}

widget_ids!{
    pub struct Ids {
        background,

        toggle_light_options,

        text_light_pos,
        xypad_light_pos_dx_dz,
        slider_light_pos_dy,

        text_light_color,
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

use piston_window::*;
use conrod_core::*;

use gfx::{
    traits::*,
    Slice,
};
use camera_controllers::CameraPerspective;

fn main() {
    const SENSATIVIY: f32 = -0.001;                     // чутливість камери при повороті
    const SPEED: f32 = 10.;                             // швидкість руху камери 
    let mut view_orientation = quaternion::id::<f32>(); // оріентація камери у виді кватерніона
    let mut scaling = 1.;                               // маштаб камери
    let mut holding_mouse_button = None;                // для того щоб поворот камери працював 
                                                        // лише при затримці клавиші на миші
    let mut mv = [0.; 6];                               // масив змінних що описують які з клавіш руху натиснуті зараз
    let mut shape_menu = false;                         // змінна стану меню - чи розгорнуте меню налаштувать фігур?
    let mut light_menu = false;                         // змінна стану меню - чи розгорнуте меню налаштувать світла?
    let mut light_follow = false;                       // змінна стану меню - чи повинне світло слідувати позаду камери?

    let opengl = OpenGL::V3_2;

    let mut window: PistonWindow = WindowSettings::new("piston: cube", [1280, 720])
        .exit_on_esc(true)
        .graphics_api(opengl)
        .build()
        .unwrap();

    let ref mut factory = window.factory.clone();

    let pso = factory.create_pipeline_simple(
        include_bytes!("../assets/combine.glslv"),
        include_bytes!("../assets/combine.glslf"),
        pipe::new()
    ).unwrap();

    let mut data = setup(&window, scaling);

    let mut texture_context = window.create_texture_context();
    let (mut ui, mut glyph_cache, mut text_texture_cache, mut text_vertex_data) = {
        let piston_window::Size {
            width,
            height
        } = window.window.draw_size();
        
        // construct our `Ui`.
        let mut ui = conrod_core::UiBuilder::new([width, height])
            // .theme(conrod_example_shared::theme())
            .build();
    
        // Add a `Font` to the `Ui`'s `font::Map` from file.
    
        ui.fonts.insert_from_file("D:\\Projects\\oop-rgr\\assets\\FiraSans-Regular.ttf").unwrap();
    
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
            let texture =
                G2dTexture::from_memory_alpha(&mut texture_context, &init, width as u32, height as u32, &settings)
                    .unwrap();
            (cache, texture)
        };
        
        (ui, glyph_cache, text_texture_cache, text_vertex_data)
    };
    let ids = Ids::new(ui.widget_id_generator());
    let image_map = conrod_core::image::Map::new();

    while let Some(e) = window.next() {
        let right = quaternion::rotate_vector(view_orientation, [1., 0., 0.]);

        // напрямки по яких рухатись при натиску клавіш
        let mv_up    = [0., 1., 0.];
        let mv_right = [right[0], 0., right[2]];
        let mv_front = {
            let front = quaternion::rotate_vector(view_orientation, [0., 0., 1.]);
            [front[0], 0., front[2]]
        };

        if holding_mouse_button != None {
            e.mouse_relative(|d| {
                let q_x = quaternion::axis_angle::<f32>(
                    [0., 1., 0.],
                    d[0] as f32 / scaling * SENSATIVIY
                );
                let q_y = quaternion::axis_angle::<f32>(
                    right,
                    d[1] as f32 / scaling * SENSATIVIY
                );
                let q_z = quaternion::rotation_from_to(
                    right,
                    mv_right
                );

                view_orientation = quaternion::mul(q_x, view_orientation);
                view_orientation = quaternion::mul(q_y, view_orientation);
                view_orientation = quaternion::mul(q_z, view_orientation);
            });

            let rot = quat_to_mat3(view_orientation);

            (0..3).for_each(|i| (0..3).for_each(|j| {
                data.u_view[j][i] = rot[i][j];
                }
            ));
        }

        e.mouse_scroll(|d| {
            scaling *= 0.95_f64.powf(-d[1]) as f32;

            data.u_view[3][3] = scaling;
        });

        match e.press_args() {
            Some(Button::Mouse(button)) => {
                holding_mouse_button = Some(button);
            },
            Some(Button::Keyboard(key)) => {
                match key {
                    Key::W     => mv[0] = SPEED,
                    Key::A     => mv[1] = SPEED,
                    Key::S     => mv[2] = SPEED,
                    Key::D     => mv[3] = SPEED,
                    Key::Space => mv[4] = SPEED,
                    Key::LShift => mv[5] = SPEED,
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
                    Key::Space => mv[4] = 0.,
                    Key::LShift => mv[5] = 0.,
                    _ => {}
                }
            },
            _ => {}
        }

        window.draw_3d(&e, |window| {
            window.encoder.draw(&Slice::new_match_vertex_buffer(&data.vbuf), &pso, &data);
        });

        // Convert the src event to a conrod event.
        let size = window.size();
        let (win_w, win_h) = (
            size.width as conrod_core::Scalar,
            size.height as conrod_core::Scalar,
        );
        if let Some(e) = conrod_piston::event::convert(e.clone(), win_w, win_h) {
            ui.handle_event(e);
        }

        e.update(|args| {
            let k = 5. * args.dt as f32;

            (0..3).for_each(|i| {
                data.u_view[3][i] += args.dt as f32 * (
                    (mv[0] - mv[2]) * mv_front[i] +
                    (mv[1] - mv[3]) * mv_right[i] +
                    (mv[5] - mv[4]) * mv_up[i]
                ) ;
            });

            if light_follow {
                data.light_pos = vecmath::vec3_mul(vecmath::vec3_add(
                    [data.u_view[3][0], data.u_view[3][1], data.u_view[3][2]], 
                    quaternion::rotate_vector(view_orientation, [0., 0., -2.])
                ), [1., -1., 1.]);
            }

            let mut ui = ui.set_widgets();

            // розставляємо всі елементи інтерфейсу
            {
                const MARGIN: conrod_core::Scalar = 30.0;
                let mut bg_height = 2. * MARGIN + 130.;
                
                bg_height += if light_menu {
                    2. * (20. + 20. + 20. + 50.)
                } else { 0. };
                bg_height += if shape_menu {
                    3. * (20. + 20. + 20. + 50.)
                } else { 0. };

                widget::bordered_rectangle::BorderedRectangle::new([200., bg_height])
                    .top_left_with_margin_on(ui.window, MARGIN)
                    .set(ids.background, &mut ui);

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
                        .align_right_of(ids.slider_plane_pos_dy)
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
                        .align_right_of(ids.slider_plane_pos_dy)
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

fn resize(
    window: &piston_window::PistonWindow,
    data: &mut pipe::Data<gfx_device_gl::Resources>
)
{
    let piston_window::Size {
        width,
        height
    } = window.window.draw_size();

    data.u_proj = {
        CameraPerspective {
            fov: 60.0, near_clip: 0.1, far_clip: 10.0,
            aspect_ratio: (width / height) as f32
        }.projection()
    };
    data.u_res = [width as f32, height as f32];
    data.out_color = window.output_color.clone();
}

fn setup(window: &piston_window::PistonWindow, scaling: f32) -> pipe::Data<gfx_device_gl::Resources> {
    let piston_window::Size {
        width,
        height
    } = window.window.draw_size();

    let ref mut factory = window.factory.clone();

    let data = pipe::Data {

        // малювати 2 трикутника, що покривають весь екран
        vbuf: factory.create_vertex_buffer(&[
            Vertex { pos: [ 1,  1] },
            Vertex { pos: [-1,  1] },
            Vertex { pos: [ 1, -1] },
            Vertex { pos: [-1, -1] },
            Vertex { pos: [-1,  1] },
            Vertex { pos: [ 1, -1] }
        ]),

        // матриця перспективної проекції з кутом видимости 60 градусів
        u_proj: {
            CameraPerspective {
                fov: 60.0, near_clip: 0.1, far_clip: 10.0,
                aspect_ratio: (width / height) as f32
            }.projection()
        },

        // матриця виду з камери - містить в собі оріентацію, позицію та маштаб камери
        u_view: [
            [1., 0.,  0.,      0.],
            [0., 1.,  0.,      0.],
            [0., 0.,  1.,      0.],
            [0., 0., -6., scaling]
        ],

        // розміри вікна
        u_res:           [width as f32, height as f32],

        // параметри сцени
        light_pos:       [2., 0.,  1.],
        light_color:     [1., 1., 1.],
        sphere_center:   [-1.,  0., 0.],
        plane_center:    [ 0., -1., 0.],
        cylinder_center: [ 1.,  0., 4.],

        // буфер в який записати результат - буфер вікна додатка
        out_color:      window.output_color.clone(),
    };

    data
}