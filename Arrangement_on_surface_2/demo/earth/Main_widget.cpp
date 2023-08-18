// Copyright(c) 2012, 2020  Tel - Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Main_widget.h"

#include <cmath>
#include <iostream>
#include <string>

#include <QMouseEvent>

#include "Aos.h"
#include "Aos_triangulator.h"
#include "Camera_manip_rot.h"
#include "Camera_manip_rot_bpa.h"
#include "Camera_manip_zoom.h"
#include "GUI_country_pick_handler.h"
#include "Kml_reader.h"
#include "Message_manager.h"
#include "Shapefile.h"
#include "Timer.h"
#include "Tools.h"
#include "Verification.h"



Main_widget::~Main_widget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}

void Main_widget::hightlight_country(const std::string& country_name)
{
    static std::string  prev_picked_country;

    if (!prev_picked_country.empty())
    {
      // dim the previous country color
      auto& prev_country = m_gr_country_triangles[prev_picked_country];
      auto color = prev_country->get_color();
      color *= m_dimming_factor;
      color.setW(1);
      prev_country->set_color(color);
    }

    if (!country_name.empty())
    {
      // highlight the current country color
      auto& curr_country = m_gr_country_triangles[country_name];
      auto color = curr_country->get_color();
      color /= m_dimming_factor;
      color.setW(1);
      curr_country->set_color(color);
      qDebug() << "SELECTED COUNTRY: " << country_name;
    }

    prev_picked_country = country_name;
}
void Main_widget::mousePressEvent(QMouseEvent* e)
{
  // forward the event to the camera manipulators
  m_camera_manip_rot->mousePressEvent(e);
  m_camera_manip_zoom->mousePressEvent(e);
  m_pick_handler->mousePressEvent(e);
}
void Main_widget::mouseMoveEvent(QMouseEvent* e)
{
  // forward the event to the camera manipulator
  m_camera_manip_rot->mouseMoveEvent(e);
  m_camera_manip_zoom->mouseMoveEvent(e);
  m_pick_handler->mouseMoveEvent(e);
}
void Main_widget::mouseReleaseEvent(QMouseEvent* e)
{
  // forward the event to the camera manipulator
  m_camera_manip_rot->mouseReleaseEvent(e);
  m_camera_manip_zoom->mouseReleaseEvent(e);
  m_pick_handler->mouseReleaseEvent(e);
}
void Main_widget::timerEvent(QTimerEvent*)
{
  update();
}


void Main_widget::keyPressEvent(QKeyEvent* event)
{
}




void Main_widget::initializeGL()
{
  m_pick_handler = std::make_unique<GUI_country_pick_handler>(*this);

  m_mouse_pos = QVector3D(0, -1, 0);
  m_gr_mouse_vertex = std::make_unique<SingleVertex>(m_mouse_pos);

  // triangulation
  {
    qDebug() << "loading arr..";
    //auto arrh = Aos::construct(m_countries);
    m_arrh = Aos::load_arr("C:/work/gsoc2023/ne_110m_admin_0_countries.json");
    if (m_arrh == nullptr)
    {
      qDebug() << "** FAILED TO LOAD THE ARRANGEMENT!!!";
      exit(1);
    }

    qDebug() << "generating triangles..";
    //auto triangle_points = Aos::get_triangles(arrh);
    //auto triangle_points = Aos_triangulator::get_all(arrh);
    //auto country_triangles_map = Aos::get_triangles_by_country(m_arrh);
    auto country_triangles_map = Aos_triangulator::get_by_country(m_arrh);
    //auto color_map = Aos::get_color_mapping(m_arrh);
    //qDebug() << "color map size = " << color_map.size();
    qDebug() << "num countries = " << country_triangles_map.size();
    auto rndm = [] {return rand() / double(RAND_MAX); };
    for (auto& [country_name, triangle_points] : country_triangles_map)
    {
      auto country_triangles = std::make_unique<Triangles>(triangle_points);
      auto color = QVector4D(rndm(), rndm(), rndm(), 1);
      auto m = std::max(color.x(), std::max(color.y(), color.z()));
      color /= m;
      color *= m_dimming_factor;
      color.setW(1);
      country_triangles->set_color(color);
      //country_triangles->set_color(colors[color_map[country_name]]);
      m_gr_country_triangles.emplace(country_name, std::move(country_triangles));
    }
    
    //qDebug() << "num triangles = " << triangle_points.size() / 3;
    //m_gr_all_triangles = std::make_unique<Triangles>(triangle_points);
  }


  initializeOpenGLFunctions();

  init_camera();
  init_geometry();
  init_shader_programs();

  m_current_approx_error = 0.001;
  init_country_borders(m_current_approx_error);

  glClearColor(0, 0, 0, 1);
  glEnable(GL_DEPTH_TEST);  // Enable depth buffer
  //glEnable(GL_CULL_FACE); // Enable back face culling

  // Use QBasicTimer because its faster than QTimer
  m_timer.start(12, this);
}


void Main_widget::init_camera()
{
  m_camera.set_pos(0, 0, 3);
  m_camera_manip_rot = std::make_unique<Camera_manip_rot>(m_camera);
  //m_camera_manip_rot = std::make_unique<Camera_manip_rot_bpa>(m_camera);
  m_camera_manip_zoom = std::make_unique<Camera_manip_zoom>(m_camera);
  
  // this makes z-axes point upwards!
  m_model.rotate(-90, 1, 0, 0);

  // register the zoom-changed function
  Message_manager::add("zoom_changed", [&] 
    { 
      qDebug() << "ZOOM CHANGED!!!"; 
      //const auto error = compute_backprojected_error(0.5);
      //qDebug() << "new error = " << error;
      m_update_approx_error = true;
      //qDebug() << "re-initializing the country borders..";
      //init_country_borders(error);
    });
}
void Main_widget::init_geometry()
{
  // SPHERE
  int num_slices, num_stacks;
  num_slices = num_stacks = 64;
  float r = 1;
  m_gr_sphere = std::make_unique<Sphere>(num_slices, num_stacks, r);
  const float c = 0.8;
  m_gr_sphere->set_color(c, c, c, 1);

  // IDENTIFICATION CURVE
  const double error = 0.001;
  auto approx_ident_curve = Aos::get_approx_identification_curve(error);
  m_gr_identification_curve = std::make_unique<Line_strips>(approx_ident_curve);

  const float axes_length = 2;
  m_gr_world_coord_axes = std::make_unique<World_coord_axes>(axes_length);
}
void Main_widget::init_shader_programs()
{
  Shader_program::set_shader_path("shaders/");
  m_sp_smooth.init_with_vs_fs("smooth");;
  m_sp_per_vertex_color.init_with_vs_fs("per_vertex_color");
  m_sp_arc.init_with_vs_fs("arc");
}

void Main_widget::init_country_borders(float error)
{
  // this part does the same as the code below but using arrangement!
  // NOTE: the old code interferes with some logic (NEEDS REFACTORING!!!)
  m_gr_all_country_borders.reset(nullptr);
  qDebug() << "approximating the arcs of each edge of all faces..";
  auto all_approx_arcs = Aos::get_approx_arcs_from_faces_edges(m_arrh, error);
  m_gr_all_country_borders = std::make_unique<Line_strips>(all_approx_arcs);
}


float Main_widget::compute_backprojected_error(float pixel_error)
{
  // compute the back-projected error
  QRect vp(0, 0, m_vp_width, m_vp_height);
  auto proj = m_camera.get_projection_matrix();
  auto view = m_camera.get_view_matrix();
  QMatrix4x4 model;
  auto model_view = view * model;

  QVector3D p0(m_vp_width / 2, m_vp_height / 2, 0);
  QVector3D p1(p0.x() + pixel_error, p0.y(), 0);
  auto wp0 = p0.unproject(model_view, proj, vp);
  auto wp1 = p1.unproject(model_view, proj, vp);
  const float z_near = m_camera.get_z_near();
  const float r = 1.f; // sphere radius
  const QVector3D origin(0, 0, 0);
  const float dist_to_cam = m_camera.get_pos().distanceToPoint(origin);

  float d = dist_to_cam - r;
  float err = wp0.distanceToPoint(wp1) * (d / z_near);
  //find_minimum_projected_error_on_sphere(err);
  return err;
}


void Main_widget::resizeGL(int w, int h)
{
  m_camera_manip_rot->resizeGL(w, h);
  m_pick_handler->resizeGL(w, h);
  
  m_vp_width = w;
  m_vp_height = h;

  // Reset projection
  qreal aspect = qreal(w) / qreal(h ? h : 1);
  const qreal z_near = 0.1, z_far = 100.0, fov = 45.0;
  m_camera.perspective(fov, aspect, z_near, z_far);

  // signal to look into the approximation error
  m_update_approx_error = true;
}

template<typename T>
void draw_safe(T& ptr)
{
  if (ptr)
    ptr->draw();
}

void Main_widget::paintGL()
{
  if (m_update_approx_error)
  {
    const auto error = compute_backprojected_error(0.5);
    qDebug() << "new approx error = " << error;
    qDebug() << "current error = " << m_current_approx_error;
    if(error < m_current_approx_error)
    {
      init_country_borders(error);
      m_current_approx_error = error;
    }
    m_update_approx_error = false;
  }

  const auto view = m_camera.get_view_matrix();
  const auto projection = m_camera.get_projection_matrix();
  const auto model_view = view * m_model;
  const auto mvp = projection * model_view;
  const auto normal_matrix = model_view.normalMatrix();

  // compute the cutting plane
  // remember that we are passing the local vertex positions of the sphere 
  // between the vertex and fragment shader stages, so we need to convert
  // the camera-pos in world coords to sphere's local coords!
  auto c = m_model.inverted() * m_camera.get_pos();
  const auto d = c.length();
  const auto r = 1.0f;
  const auto sin_alpha = r / d;
  const auto n = (c / d); // plane unit normal vector
  const auto cos_beta = sin_alpha;
  const auto p = (r * cos_beta) * n;
  QVector4D plane(n.x(), n.y(), n.z(), -QVector3D::dotProduct(p, n));

  // Clear color and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // SMOTH RENDERING
  {
    glEnable(GL_DEPTH_TEST);

    auto& sp = m_sp_smooth;
    sp.use();

    // SPHERE
    {
      sp.set_uniform("u_mvp", mvp);
      sp.set_uniform("u_normal_matrix", normal_matrix);
      auto sphere_color = QVector4D(167, 205, 242, 255) / 255;
      sp.set_uniform("u_color", sphere_color);
      sp.set_uniform("u_plane", QVector4D(0, 0, 0, 0));
      //sp.set_uniform("u_color", m_sphere->get_color());

      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      m_gr_sphere->draw();
    }

    // DRAW SOLID FACES
    {
      glDisable(GL_DEPTH_TEST);
      //auto face_color = QVector4D(241, 141, 0, 255) / 255;
      //sp.set_uniform("u_color", face_color);
      sp.set_uniform("u_plane", plane);
      //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      //m_gr_all_triangles->draw();
      for (auto& [country_name, country] : m_gr_country_triangles)
      {
        sp.set_uniform("u_color", country->get_color());
        country->draw();
      }

      //sp.set_uniform("u_color", QVector4D(0,0,0,1));
      //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      //m_gr_all_triangles->draw();
    }

    sp.unuse();
  }

  // WORLD COORDINATE AXES
  {
    auto& sp = m_sp_per_vertex_color;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    glEnable(GL_DEPTH_TEST);
    m_gr_world_coord_axes->draw();

    sp.unuse();
  }

  // VERTICES & GEODESIC ARCS
  {
    glDisable(GL_DEPTH_TEST);

    auto& sp = m_sp_arc;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    const QVector4D arc_color(0, 0.5, 1, 1);
    glLineWidth(5);
    sp.set_uniform("u_plane", plane);

    // IDENTIFICATION CURVE
    sp.set_uniform("u_color", QVector4D(0, 1, 1, 1));
    m_gr_identification_curve->draw();

    // draw all countries' borders 
    float a = 0.0;
    sp.set_uniform("u_color", QVector4D(a, a, a, 1));
    m_gr_all_country_borders->draw();

    // MOUSE VERTEX
    {
      glPointSize(5);
      sp.set_uniform("u_color", QVector4D(1, 0, 0, 1));
      //auto pos = m_mouse_vertex->get_pos();
      //pos.setX(pos.x() + 0.01);
      //m_mouse_vertex->set_pos(pos);
      m_gr_mouse_vertex->set_pos(m_mouse_pos);
      draw_safe(m_gr_mouse_vertex);
    }

    sp.unuse();
  }
}
