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
#include "Kml_reader.h"
#include "Message_manager.h"
#include "Shapefile.h"
#include "Timer.h"
#include "Tools.h"
#include "Verification.h"


namespace
{
  // used when dimming / highlighting selected countries
  float s_dimming_factor = 0.4;
}

Main_widget::~Main_widget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}

void Main_widget::handle_country_picking(QMouseEvent* e)
{
  //// handle country selection
  //if (e->button() == Qt::RightButton)
  //{
  //  auto p = e->pos();
  //  QVector3D  sp0(p.x(), m_vp_height - p.y(), 0);
  //  QVector3D  sp1(p.x(), m_vp_height - p.y(), 1);

  //  auto proj = m_camera.get_projection_matrix();
  //  auto view = m_camera.get_view_matrix();
  //  auto model_view = view * m_model;
  //  QRect viewport(0, 0, m_vp_width, m_vp_height);
  //  auto wp0 = sp0.unproject(model_view, proj, viewport);
  //  auto wp1 = sp1.unproject(model_view, proj, viewport);

  //  // ASSERTION!!!
  //  m_mouse_pos = wp0;

  //  // define a ray from the camera pos to the world-point
  //  //auto o = m_camera.get_pos();
  //  //auto u = wp - o;
  //  auto o = wp0;
  //  auto u = wp1 - wp0;

  //  // solve the quadratic equation to check for intersection of ray with sphere
  //  auto a = QVector3D::dotProduct(u, u);
  //  auto b = 2 * QVector3D::dotProduct(u, o);
  //  auto c = QVector3D::dotProduct(o, o) - 1;
  //  auto d = b * b - 4 * a * c;

  //  float ti = -1;
  //  if (abs(d) < std::numeric_limits<float>::epsilon())
  //  {
  //    // single intersection
  //    ti = -b / (2 * a);
  //  }
  //  else
  //  {
  //    if (d < 0)
  //    {
  //      // no intersection
  //      return;
  //    }
  //    else
  //    {
  //      // two intersections
  //      auto sd = sqrt(d);
  //      auto t1 = (-b - sd) / (2 * a);
  //      auto t2 = (-b + sd) / (2 * a);
  //      if (t1 > 0 && t2 > 0)
  //        ti = std::min(t1, t2);
  //      else if (t1 > 0)
  //        ti = t1;
  //      else
  //        ti = t2;
  //    }
  //  }

  //  m_mouse_pos = o + ti * u;
  //  static std::string prev_picked_country;
  //  auto picked_country = Aos::locate_country(m_arrh, m_mouse_pos);

  //  if (!prev_picked_country.empty())
  //  {
  //    // dim the previous country color
  //    auto& prev_country = m_country_triangles[prev_picked_country];
  //    auto color = prev_country->get_color();
  //    color *= s_dimming_factor;
  //    color.setW(1);
  //    prev_country->set_color(color);
  //  }

  //  if (!picked_country.empty())
  //  {
  //    // highlight the current country color
  //    auto& curr_country = m_country_triangles[picked_country];
  //    auto color = curr_country->get_color();
  //    color /= s_dimming_factor;
  //    color.setW(1);
  //    curr_country->set_color(color);
  //    qDebug() << "SELECTED COUNTRY: " << picked_country;
  //  }

  //  prev_picked_country = picked_country;
  //}
}
void Main_widget::hightlight_country(const std::string& country_name)
{
    static std::string  prev_picked_country;

    if (!prev_picked_country.empty())
    {
      // dim the previous country color
      auto& prev_country = m_country_triangles[prev_picked_country];
      auto color = prev_country->get_color();
      color *= s_dimming_factor;
      color.setW(1);
      prev_country->set_color(color);
    }

    if (!country_name.empty())
    {
      // highlight the current country color
      auto& curr_country = m_country_triangles[country_name];
      auto color = curr_country->get_color();
      color /= s_dimming_factor;
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
  //handle_country_picking(e);
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
  switch (event->key())
  {
  case Qt::Key_Q:
  {
    auto num_arcs = m_country_borders[m_selected_country_index]->get_num_line_strips();
    if (++m_selected_arc_index == num_arcs)
      m_selected_arc_index--;
    qDebug() << "---------------------------------------";
    qDebug() << "selected arc index = " << m_selected_arc_index;
    
    const auto& arc = m_selected_country_arcs[m_selected_arc_index];
    std::cout << arc.from << "  TO  " << arc.to << std::endl;
  }
    break;
  case Qt::Key_A:
  {
    auto num_arcs = m_country_borders[m_selected_country_index]->get_num_line_strips();
    if (--m_selected_arc_index < 0)
      m_selected_arc_index = 0;
    std::cout << "selected arc = " << m_selected_arc_index << std::endl;
  }
    break;

  case Qt::Key_Up:
    m_selected_country_index++;
    if (m_selected_country_index == m_country_names.size())
      m_selected_country_index--;
    std::cout << m_selected_country_index << ": " 
              << m_country_names[m_selected_country_index] << std::endl;
    
    m_selected_arc_index = 0;
    m_selected_country = &m_countries[m_selected_country_index];
    m_selected_country_nodes = m_selected_country->get_all_nodes();
    m_selected_country_arcs = m_selected_country->get_all_arcs();

    {
      auto num_arcs = m_country_borders[m_selected_country_index]->get_num_line_strips();
      
      qDebug() << "num KML arcs = " << m_selected_country_arcs.size();
      qDebug() << "num arcs = " << num_arcs;
    }
    break;

  case Qt::Key_Down:
    m_selected_country_index--;
    if (m_selected_country_index < 0)
      m_selected_country_index = 0;
    std::cout << m_selected_country_index << ": " 
              << m_country_names[m_selected_country_index] << std::endl;
    
    m_selected_arc_index = 0;
    m_selected_country = &m_countries[m_selected_country_index];
    m_selected_country_nodes = m_selected_country->get_all_nodes();
    break;
  }
}


void Main_widget::init_problematic_nodes()
{
  Kml::Nodes prob_nodes = {
    {23.8058134294668,8.66631887454253},
    {24.1940677211877,8.7286964724039 },
    {24.5673690121521,8.22918793378547},
    {23.8869795808607,8.61972971293307}
  };
  std::vector<QVector3D> prob_vertices;
  for (const auto& node : prob_nodes)
    prob_vertices.push_back(node.get_coords_3f());
  m_problematic_vertices = std::make_unique<Vertices>(prob_vertices);
}

#include "GUI_country_pick_handler.h"

void Main_widget::initializeGL()
{
  m_pick_handler = std::make_unique<GUI_country_pick_handler>(*this);

  // verify that the node (180.0, -84.71338) in Antarctica is redundant
  //Verification::verify_antarctica_node_is_redundant();

  //init_problematic_nodes();
  m_mouse_pos = QVector3D(0, -1, 0);
  m_mouse_vertex = std::make_unique<SingleVertex>(m_mouse_pos);


  std::string data_path = "C:/work/gsoc2023/data/";
  //std::string shape_file_path = data_path + "ne_110m_admin_0_countries/";
  //auto shape_file_name = shape_file_path + "ne_110m_admin_0_countries.shp";
  //Shapefile::read(shape_file_name);

  //const auto file_name = data_path + "world_countries.kml";
  //const auto file_name = data_path + "ne_110m_admin_0_countries.kml";
  const auto file_name = data_path + "ne_110m_admin_0_countries_africa.kml";
  //const auto file_name = data_path + "ne_110m_admin_0_countries_equatorial_guinea.kml";
  m_countries = Kml::read(file_name);
  
  // find the country with the least number of nodes
  if(0)
  {
    std::string smallest;
    int min_num_nodes = std::numeric_limits<int>::max();
    for (auto& p : m_countries)
    {
      int num_nodes = p.get_all_nodes_count();
      if (min_num_nodes > num_nodes)
      {
        min_num_nodes = num_nodes;
        smallest = p.name;
      }
      qDebug() << p.name << " = " << p.get_all_nodes_count();
    }
    qDebug() << "smallest = " << smallest;
    exit(0);
  }
  
  auto dup_nodes = Kml::get_duplicates(m_countries);
  //auto all_nodes = Kml::generate_ids(m_countries);
  qDebug() << "*** KML number of polygons = " << 
                                      Kml::get_number_of_polygons(m_countries);
  if(0)
  {
    auto created_faces = Aos::find_new_faces(m_countries);
    m_new_faces = std::make_unique<Line_strips>(created_faces);
  }

  // SAVING ARR
  if(0)
  {
    std::string dest_path = "C:/work/gsoc2023/";
    //std::string file_name = "ne_110m_admin_0_countries.json";
    //std::string file_name = "ne_110m_admin_0_countries_africa_1.json";
    std::string file_name = "ne_110m_admin_0_countries_equatorial_guinea.json";
    auto full_path = dest_path + file_name;
    Aos::save_arr(m_countries, full_path);
    qDebug() << "done saving!";
    exit(0);
  }

  // triangulation
  {
    qDebug() << "constructiong arr..";
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
    //QVector4D colors[] = {
    //  QVector4D(1,0,0,1),
    //  QVector4D(0,1,0,1),
    //  QVector4D(0,0,1,1),
    //  QVector4D(1,1,0,1),
    //  QVector4D(1,0,1,1)
    //};
    for (auto& [country_name, triangle_points] : country_triangles_map)
    {
      auto country_triangles = std::make_unique<Triangles>(triangle_points);
      auto color = QVector4D(rndm(), rndm(), rndm(), 1);
      auto m = std::max(color.x(), std::max(color.y(), color.z()));
      color /= m;
      color *= s_dimming_factor;
      color.setW(1);
      country_triangles->set_color(color);
      //country_triangles->set_color(colors[color_map[country_name]]);
      m_country_triangles.emplace(country_name, std::move(country_triangles));
    }
    
    //qDebug() << "num triangles = " << triangle_points.size() / 3;
    //m_all_triangles = std::make_unique<Triangles>(triangle_points);
  }

  
  // initialize rendering of DUPLICATE VERTICES
  if(1)
  {
    qDebug() << "identifying duplicate nodes";
    std::vector<QVector3D> vertices;
    for (const auto& node : dup_nodes)
      vertices.push_back(node.get_coords_3f());
  
    m_vertices = std::make_unique<Vertices>(vertices);
  }
  else
  {
    // check the arrangement constructed from the GIS data-set
    auto created_vertices = Aos::ext_check(m_countries);
    //auto created_vertices = Aos::ext_check_id_based(m_countries);
    m_vertices = std::make_unique<Vertices>(created_vertices);
  }


  initializeOpenGLFunctions();

  init_camera();
  init_geometry();
  init_shader_programs();

  m_current_approx_error = 0.001;
  init_country_borders(m_current_approx_error);
  init_country_selection();

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
  m_sphere = std::make_unique<Sphere>(num_slices, num_stacks, r);
  const float c = 0.8;
  m_sphere->set_color(c, c, c, 1);

  // IDENTIFICATION CURVE
  const double error = 0.001;
  auto approx_ident_curve = Aos::get_approx_identification_curve(error);
  m_identification_curve = std::make_unique<Line_strips>(approx_ident_curve);

  const float axes_length = 2;
  m_world_coord_axes = std::make_unique<World_coord_axes>(axes_length);
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
  // TO-DO: move this code to resizeGL (when viewport is initialized)
  // has to be defined after camera has been defined:
  // because we want to compute the error based on camera parameters!
  //const double error = 0.001; // calculate this from cam parameters!
  //auto lsa = Aos::get_approx_arcs(countries, error);
  //auto lsa = Aos::get_approx_arcs(error);
  //m_geodesic_arcs = std::make_unique<Line_strips>(lsa);
  m_country_borders.clear();
  for (const auto& country : m_countries)
  {
    m_country_names.push_back(country.name);
    auto approx_arcs = Aos::get_approx_arcs(country, error);
    auto country_border = std::make_unique<Line_strips>(approx_arcs);
    m_country_borders.push_back(std::move(country_border));
  }
}
void Main_widget::init_country_selection()
{
  m_selected_country_index = 0;
  //m_selected_country_index = 159; // ANTARCTICA
  m_selected_country = &m_countries[m_selected_country_index];
  m_selected_country_nodes = m_selected_country->get_all_nodes();
  m_selected_country_arcs = m_selected_country->get_all_arcs();
  m_selected_arc_index = 0;
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
  // SPHERE
  {
    glEnable(GL_DEPTH_TEST);

    auto& sp = m_sp_smooth;
    sp.use();
    sp.set_uniform("u_mvp", mvp);
    sp.set_uniform("u_normal_matrix", normal_matrix);
    auto sphere_color = QVector4D(167, 205, 242, 255) / 255;
    sp.set_uniform("u_color", sphere_color);
    sp.set_uniform("u_plane", QVector4D(0,0,0,0));
    //sp.set_uniform("u_color", m_sphere->get_color());
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    m_sphere->draw();

    // DRAW SOLID FACES
    if(1)
    {
      glDisable(GL_DEPTH_TEST);
      //auto face_color = QVector4D(241, 141, 0, 255) / 255;
      //sp.set_uniform("u_color", face_color);
      sp.set_uniform("u_plane", plane);
      //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      //m_all_triangles->draw();
      for (auto& [country_name, country] : m_country_triangles)
      {
        sp.set_uniform("u_color", country->get_color());
        country->draw();
      }

      //sp.set_uniform("u_color", QVector4D(0,0,0,1));
      //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      //m_all_triangles->draw();
    }

    sp.unuse();
  }

  // WORLD COORDINATE AXES
  {
    auto& sp = m_sp_per_vertex_color;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    glEnable(GL_DEPTH_TEST);
    m_world_coord_axes->draw();

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
    m_identification_curve->draw(m_selected_arc_index);

    // draw all countries 
    //float a = 0.0;
    //sp.set_uniform("u_color", QVector4D(a, a, a, 1));
    //for(auto& country_border : m_country_borders)
    //  country_border->draw();

    //// draw the SELECTED COUNTRY in BLUE
    //auto& selected_country = m_country_borders[m_selected_country_index];
    //sp.set_uniform("u_color", QVector4D(0, .6, 1, 1));
    //selected_country->draw();

    //// draw the CURRENT ARC of the selected country in YELLOW
    //sp.set_uniform("u_color", QVector4D(1, 1, 0, 1));
    //selected_country->draw(m_selected_arc_index);

    //const QVector4D vertex_color(1, 0, 0, 1);
    //sp.set_uniform("u_color", vertex_color);
    //glPointSize(3);
    ////m_vertices->draw();

    //sp.set_uniform("u_color", QVector4D(0,1,0,1));
    //glPointSize(2);
    ////m_problematic_vertices->draw();
    //draw_safe(m_problematic_vertices);

    //// NEW FACES in RED
    //sp.set_uniform("u_color", QVector4D(1, 0, 0, 1));
    ////m_new_faces->draw();

    // MOUSE VERTEX
    {
      glPointSize(5);
      sp.set_uniform("u_color", QVector4D(1, 0, 0, 1));
      //auto pos = m_mouse_vertex->get_pos();
      //pos.setX(pos.x() + 0.01);
      //m_mouse_vertex->set_pos(pos);
      m_mouse_vertex->set_pos(m_mouse_pos);
      draw_safe(m_mouse_vertex);
    }

    sp.unuse();
  }
}
