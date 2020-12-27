extern crate piston_window;
extern crate vecmath;
extern crate camera_controllers;

#[macro_use]
extern crate gfx;
extern crate shader_version;
extern crate quaternion;
// #[macro_use] extern crate conrod;
#[macro_use]
extern crate conrod_core;
extern crate conrod_piston;

extern crate gfx_device_gl;

use self::piston_window::texture::UpdateTexture;
// use self::piston_window::OpenGL;
// use self::piston_window::{Flip, G2d, G2dTexture, Texture, TextureSettings};

mod load_mesh;
// use load_mesh::*;

gfx_defines!{
    vertex Vertex {
        pos: [i8; 2] = "a_pos",
    }

    pipeline pipe {
        vbuf:        gfx::VertexBuffer  <  Vertex                   > = (),
        u_proj:      gfx::Global        <  [[f32; 4]; 4]            > = "u_proj",
        u_view:      gfx::Global        <  [[f32; 4]; 4]            > = "u_view",
        light_pos:   gfx::Global        <  [f32; 3]                 > = "light_pos",
        light_color: gfx::Global        <  [f32; 3]                 > = "light_color",
        u_res:       gfx::Global        <  [f32; 2]                 > = "u_resolution",
        out_color:   gfx::RenderTarget  <::gfx::format::Srgba8      > = "frag_color",
    }
}

widget_ids!{
    pub struct Ids {
        button_light,
        button_light_pos,
        button_light_color,
        button_change_quality,
        button_change_quality_max_steps,
        button_change_quality_max_dist,
        button_change_quality_min_dist,
        button_change_quality_light_soft_samples,
        toggle_light_follow,
    }
}

//----------------------------------------

use piston_window::*;
use gfx::{
    traits::*,
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

    let pso = factory.create_pipeline_simple(
        include_bytes!("../assets/combine.glslv"),
        include_bytes!("../assets/combine.glslf"),
        pipe::new()
    ).unwrap();

    let mut view_orientation = quaternion::id::<f32>();
    let mut scaling = 1.;

    let mut data = setup(&window, scaling);
    // let (mut ui, mut image_map, mut texture_context, mut glyph_cache, mut text_texture_cache, mut text_vertex_data)
    //   = setup_ui(&mut window);
    let mut ui = setup_ui(&mut window);

    let mut holding_mouse_button = None;
    let mut mv = [0.; 6];

    while let Some(e) = window.next() {
        let right = quaternion::rotate_vector(view_orientation, [1., 0., 0.]);
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
            (0..3).for_each(|i| {
                data.u_view[3][i] +=
                    (mv[0] - mv[2]) * mv_front[i] +
                    (mv[1] - mv[3]) * mv_right[i] +
                    (mv[5] - mv[4]) * mv_up[i];
            });
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
            // (0..3).for_each(|i| {
            //     data.u_view[3][i] += args.dt as f32 * (
            //         (mv[0] - mv[2]) * mv_front[i] +
            //         (mv[1] - mv[3]) * mv_right[i] +
            //         (mv[5] - mv[4]) * mv_up[i]
            //     ) ;
            // });

            let mut ui = ui.set_widgets();

            // {
            //     use conrod_core::{widget, Colorable, Labelable, Positionable, Sizeable, Widget};
            //     use std::iter::once;

            //     const MARGIN: conrod_core::Scalar = 30.0;
            //     const SHAPE_GAP: conrod_core::Scalar = 50.0;
            //     const TITLE_SIZE: conrod_core::FontSize = 42;
            //     const SUBTITLE_SIZE: conrod_core::FontSize = 32;

            //     // `Canvas` is a widget that provides some basic functionality for laying out children widgets.
            //     // By default, its size is the size of the window. We'll use this as a background for the
            //     // following widgets, as well as a scrollable container for the children widgets.
            //     const TITLE: &'static str = "All Widgets";
            //     widget::Canvas::new()
            //         .pad(MARGIN)
            //         .scroll_kids_vertically()
            //         .set(ids.canvas, ui);

            //     ////////////////
            //     ///// TEXT /////
            //     ////////////////

            //     // We'll demonstrate the `Text` primitive widget by using it to draw a title and an
            //     // introduction to the example.
            //     widget::Text::new(TITLE)
            //         .font_size(TITLE_SIZE)
            //         .mid_top_of(ids.canvas)
            //         .set(ids.title, ui);

            //     const INTRODUCTION: &'static str =
            //         "This example aims to demonstrate all widgets that are provided by conrod.\
            //         \n\nThe widget that you are currently looking at is the Text widget. The Text widget \
            //         is one of several special \"primitive\" widget types which are used to construct \
            //         all other widget types. These types are \"special\" in the sense that conrod knows \
            //         how to render them via `conrod_core::render::Primitive`s.\
            //         \n\nScroll down to see more widgets!";
            //     widget::Text::new(INTRODUCTION)
            //         .padded_w_of(ids.canvas, MARGIN)
            //         .down(60.0)
            //         .align_middle_x_of(ids.canvas)
            //         .center_justify()
            //         .line_spacing(5.0)
            //         .set(ids.introduction, ui);

            //     ////////////////////////////
            //     ///// Lines and Shapes /////
            //     ////////////////////////////

            //     widget::Text::new("Lines and Shapes")
            //         .down(70.0)
            //         .align_middle_x_of(ids.canvas)
            //         .font_size(SUBTITLE_SIZE)
            //         .set(ids.shapes_title, ui);

            //     // Lay out the shapes in two horizontal columns.
            //     //
            //     // TODO: Have conrod provide an auto-flowing, fluid-list widget that is more adaptive for these
            //     // sorts of situations.
            //     widget::Canvas::new()
            //         .down(0.0)
            //         .align_middle_x_of(ids.canvas)
            //         .kid_area_w_of(ids.canvas)
            //         .h(360.0)
            //         .color(conrod_core::color::TRANSPARENT)
            //         .pad(MARGIN)
            //         .flow_down(&[
            //             (ids.shapes_left_col, widget::Canvas::new()),
            //             (ids.shapes_right_col, widget::Canvas::new()),
            //         ])
            //         .set(ids.shapes_canvas, ui);

            //     let shapes_canvas_rect = ui.rect_of(ids.shapes_canvas).unwrap();
            //     let w = shapes_canvas_rect.w();
            //     let h = shapes_canvas_rect.h() * 5.0 / 6.0;
            //     let radius = 10.0;
            //     widget::RoundedRectangle::fill([w, h], radius)
            //         .color(conrod_core::color::CHARCOAL.alpha(0.25))
            //         .middle_of(ids.shapes_canvas)
            //         .set(ids.rounded_rectangle, ui);

            //     let start = [-40.0, -40.0];
            //     let end = [40.0, 40.0];
            //     widget::Line::centred(start, end)
            //         .mid_left_of(ids.shapes_left_col)
            //         .set(ids.line, ui);

            //     let left = [-40.0, -40.0];
            //     let top = [0.0, 40.0];
            //     let right = [40.0, -40.0];
            //     let points = once(left).chain(once(top)).chain(once(right));
            //     widget::PointPath::centred(points)
            //         .right(SHAPE_GAP)
            //         .set(ids.point_path, ui);

            //     widget::Rectangle::fill([80.0, 80.0])
            //         .right(SHAPE_GAP)
            //         .set(ids.rectangle_fill, ui);

            //     widget::Rectangle::outline([80.0, 80.0])
            //         .right(SHAPE_GAP)
            //         .set(ids.rectangle_outline, ui);

            //     let bl = [-40.0, -40.0];
            //     let tl = [-20.0, 40.0];
            //     let tr = [20.0, 40.0];
            //     let br = [40.0, -40.0];
            //     let points = once(bl).chain(once(tl)).chain(once(tr)).chain(once(br));
            //     widget::Polygon::centred_fill(points)
            //         .mid_left_of(ids.shapes_right_col)
            //         .set(ids.trapezoid, ui);

            //     widget::Oval::fill([40.0, 80.0])
            //         .right(SHAPE_GAP + 20.0)
            //         .align_middle_y()
            //         .set(ids.oval_fill, ui);

            //     widget::Oval::outline([80.0, 40.0])
            //         .right(SHAPE_GAP + 20.0)
            //         .align_middle_y()
            //         .set(ids.oval_outline, ui);

            //     widget::Circle::fill(40.0)
            //         .right(SHAPE_GAP)
            //         .align_middle_y()
            //         .set(ids.circle, ui);

            //     /////////////////
            //     ///// Image /////
            //     /////////////////

            //     widget::Text::new("Image")
            //         .down_from(ids.shapes_canvas, MARGIN)
            //         .align_middle_x_of(ids.canvas)
            //         .font_size(SUBTITLE_SIZE)
            //         .set(ids.image_title, ui);

            //     const LOGO_SIDE: conrod_core::Scalar = 144.0;
            //     widget::Image::new(app.rust_logo)
            //         .w_h(LOGO_SIDE, LOGO_SIDE)
            //         .down(60.0)
            //         .align_middle_x_of(ids.canvas)
            //         .set(ids.rust_logo, ui);

            //     /////////////////////////////////
            //     ///// Button, XYPad, Toggle /////
            //     /////////////////////////////////

            //     widget::Text::new("Button, XYPad and Toggle")
            //         .down_from(ids.rust_logo, 60.0)
            //         .align_middle_x_of(ids.canvas)
            //         .font_size(SUBTITLE_SIZE)
            //         .set(ids.button_title, ui);

            //     let ball_x_range = ui.kid_area_of(ids.canvas).unwrap().w();
            //     let ball_y_range = ui.h_of(ui.window).unwrap() * 0.5;
            //     let min_x = -ball_x_range / 3.0;
            //     let max_x = ball_x_range / 3.0;
            //     let min_y = -ball_y_range / 3.0;
            //     let max_y = ball_y_range / 3.0;
            //     let side = 130.0;

            //     for _press in widget::Button::new()
            //         .label("PRESS ME")
            //         .mid_left_with_margin_on(ids.canvas, MARGIN)
            //         .down_from(ids.button_title, 60.0)
            //         .w_h(side, side)
            //         .set(ids.button, ui)
            //     {
            //         let x = rand::random::<conrod_core::Scalar>() * (max_x - min_x) - max_x;
            //         let y = rand::random::<conrod_core::Scalar>() * (max_y - min_y) - max_y;
            //         app.ball_xy = [x, y];
            //     }

            //     for (x, y) in widget::XYPad::new(app.ball_xy[0], min_x, max_x, app.ball_xy[1], min_y, max_y)
            //         .label("BALL XY")
            //         .wh_of(ids.button)
            //         .align_middle_y_of(ids.button)
            //         .align_middle_x_of(ids.canvas)
            //         .parent(ids.canvas)
            //         .set(ids.xy_pad, ui)
            //     {
            //         app.ball_xy = [x, y];
            //     }

            //     let is_white = app.ball_color == conrod_core::color::WHITE;
            //     let label = if is_white { "WHITE" } else { "BLACK" };
            //     for is_white in widget::Toggle::new(is_white)
            //         .label(label)
            //         .label_color(if is_white {
            //             conrod_core::color::WHITE
            //         } else {
            //             conrod_core::color::LIGHT_CHARCOAL
            //         })
            //         .mid_right_with_margin_on(ids.canvas, MARGIN)
            //         .align_middle_y_of(ids.button)
            //         .set(ids.toggle, ui)
            //     {
            //         app.ball_color = if is_white {
            //             conrod_core::color::WHITE
            //         } else {
            //             conrod_core::color::BLACK
            //         };
            //     }

            //     let ball_x = app.ball_xy[0];
            //     let ball_y = app.ball_xy[1] - max_y - side * 0.5 - MARGIN;
            //     widget::Circle::fill(20.0)
            //         .color(app.ball_color)
            //         .x_y_relative_to(ids.xy_pad, ball_x, ball_y)
            //         .set(ids.ball, ui);

            //     //////////////////////////////////
            //     ///// NumberDialer, PlotPath /////
            //     //////////////////////////////////

            //     widget::Text::new("NumberDialer and PlotPath")
            //         .down_from(ids.xy_pad, max_y - min_y + side * 0.5 + MARGIN)
            //         .align_middle_x_of(ids.canvas)
            //         .font_size(SUBTITLE_SIZE)
            //         .set(ids.dialer_title, ui);

            //     // Use a `NumberDialer` widget to adjust the frequency of the sine wave below.
            //     let min = 0.5;
            //     let max = 200.0;
            //     let decimal_precision = 1;
            //     for new_freq in widget::NumberDialer::new(app.sine_frequency, min, max, decimal_precision)
            //         .down(60.0)
            //         .align_middle_x_of(ids.canvas)
            //         .w_h(160.0, 40.0)
            //         .label("F R E Q")
            //         .set(ids.number_dialer, ui)
            //     {
            //         app.sine_frequency = new_freq;
            //     }

            //     // Use the `PlotPath` widget to display a sine wave.
            //     let min_x = 0.0;
            //     let max_x = std::f32::consts::PI * 2.0 * app.sine_frequency;
            //     let min_y = -1.0;
            //     let max_y = 1.0;
            //     widget::PlotPath::new(min_x, max_x, min_y, max_y, f32::sin)
            //         .kid_area_w_of(ids.canvas)
            //         .h(240.0)
            //         .down(60.0)
            //         .align_middle_x_of(ids.canvas)
            //         .set(ids.plot_path, ui);

            //     /////////////////////
            //     ///// Scrollbar /////
            //     /////////////////////

            //     widget::Scrollbar::y_axis(ids.canvas)
            //         .auto_hide(true)
            //         .set(ids.canvas_scrollbar, ui);
            // }
        });

        // window.draw_2d(&e, |context, graphics, device| {
        //     if let Some(primitives) = ui.draw_if_changed() {
        //         // A function used for caching glyphs to the texture cache.
        //         let cache_queued_glyphs = |_graphics: &mut G2d,
        //                                 cache: &mut G2dTexture,
        //                                 rect: conrod_core::text::rt::Rect<u32>,
        //                                 data: &[u8]| {
        //             let offset = [rect.min.x, rect.min.y];
        //             let size = [rect.width(), rect.height()];
        //             let format = piston_window::texture::Format::Rgba8;
        //             text_vertex_data.clear();
        //             text_vertex_data.extend(data.iter().flat_map(|&b| vec![255, 255, 255, b]));
        //             UpdateTexture::update(
        //                 cache,
        //                 &mut texture_context,
        //                 format,
        //                 &text_vertex_data[..],
        //                 offset,
        //                 size,
        //             )
        //             .expect("failed to update texture")
        //         };

        //         // Specify how to get the drawable texture from the image. In this case, the image
        //         // *is* the texture.
        //         fn texture_from_image<T>(img: &T) -> &T {
        //             img
        //         }

        //         // Draw the conrod `render::Primitives`.
        //         conrod_piston::draw::primitives(
        //             primitives,
        //             context,
        //             graphics,
        //             &mut text_texture_cache,
        //             &mut glyph_cache,
        //             &image_map,
        //             cache_queued_glyphs,
        //             texture_from_image,
        //         );

        //         texture_context.encoder.flush(device);
        //     }
        // });

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

fn resize(window: &piston_window::PistonWindow,
    data: &mut pipe::Data<gfx_device_gl::Resources>)
{
    let piston_window::Size {
        width,
        height
    } = window.window.draw_size();

    let ref mut factory = window.factory.clone();

    data.u_proj = {
        CameraPerspective {
            fov: 60.0, near_clip: 0.1, far_clip: 10.0,
            aspect_ratio: (width / height) as f32
        }.projection()
    };
    data.u_res = [width as f32, height as f32];
    data.out_color = window.output_color.clone();
}

fn setup(
    window: &piston_window::PistonWindow,
    scaling: f32
) -> pipe::Data<gfx_device_gl::Resources>
{
    let piston_window::Size {
        width,
        height
    } = window.window.draw_size();

    let ref mut factory = window.factory.clone();

    let data = pipe::Data {
        vbuf:           factory.create_vertex_buffer(&[
                Vertex { pos: [ 1,  1] },
                Vertex { pos: [-1,  1] },
                Vertex { pos: [ 1, -1] },
                Vertex { pos: [-1, -1] },
                Vertex { pos: [-1,  1] },
                Vertex { pos: [ 1, -1] }
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
        u_res:          [width as f32, height as f32],
        light_pos:      [2., 0.,  1.],
        light_color:    [0xfc as f32 / 255., 0x0f as f32 / 255., 0xc0 as f32 / 255.],
        out_color:      window.output_color.clone(),
    };

    data
}

fn setup_ui(window: &mut piston_window::PistonWindow) -> conrod_core::Ui
// (
//     conrod_core::Ui,
//     conrod_core::image::Map<piston_window::Texture<gfx_device_gl::Resources>>,
//     piston_window::TextureContext<gfx_device_gl::Factory, gfx_device_gl::Resources, gfx_device_gl::CommandBuffer>,
//     conrod_core::text::GlyphCache,
//     piston_window::Texture<gfx_device_gl::Resources>,
//     Vec<u8>
// )
{
    let piston_window::Size {
        width,
        height
    } = window.window.draw_size();

    // construct our `Ui`.
    let mut ui = conrod_core::UiBuilder::new([width, height])
        // .theme(conrod_example_shared::theme())
        .build();

    // Add a `Font` to the `Ui`'s `font::Map` from file.

    // ui.fonts.insert_from_file("../assets/FiraSans-regular.ttf").unwrap();

    // Create texture context to perform operations on textures.
    // let mut texture_context = (*window).create_texture_context();

    // Create a texture to use for efficiently caching text on the GPU.
    // let text_vertex_data = Vec::new();
    // let (glyph_cache, text_texture_cache) = {
    //     const SCALE_TOLERANCE: f32 = 0.1;
    //     const POSITION_TOLERANCE: f32 = 0.1;
    //     let cache = conrod_core::text::GlyphCache::builder()
    //         .dimensions(width as u32, height as u32)
    //         .scale_tolerance(SCALE_TOLERANCE)
    //         .position_tolerance(POSITION_TOLERANCE)
    //         .build();
    //     let buffer_len = width as usize * width as usize;
    //     let init = vec![128; buffer_len];
    //     let settings = TextureSettings::new();
    //     let texture =
    //         G2dTexture::from_memory_alpha(&mut texture_context, &init, width as u32, height as u32, &settings)
    //             .unwrap();
    //     (cache, texture)
    // };

    // Instantiate the generated list of widget identifiers.
    let ids = Ids::new(ui.widget_id_generator());

    // Create our `conrod_core::image::Map` which describes each of our widget->image mappings.
    // let image_map = conrod_core::image::Map::new();

    // A demonstration of some state that we'd like to control with the App.
    // let mut app = conrod_example_shared::DemoApp::new(rust_logo);

    // (ui, image_map, texture_context, glyph_cache, text_texture_cache, text_vertex_data)
    ui
}