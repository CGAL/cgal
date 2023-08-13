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
#include "Camera_manip_rot.h"
#include "Camera_manip_rot_bpa.h"
#include "Camera_manip_zoom.h"
#include "Kml_reader.h"
#include "Message_manager.h"
#include "Shapefile.h"
#include "Timer.h"
#include "Tools.h"


Main_widget::~Main_widget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}

void Main_widget::mousePressEvent(QMouseEvent* e)
{
  // forward the event to the camera manipulators
  m_camera_manip_rot->mousePressEvent(e);
  m_camera_manip_zoom->mousePressEvent(e);
}
void Main_widget::mouseMoveEvent(QMouseEvent* e)
{
  // forward the event to the camera manipulator
  m_camera_manip_rot->mouseMoveEvent(e);
  m_camera_manip_zoom->mouseMoveEvent(e);
}
void Main_widget::mouseReleaseEvent(QMouseEvent* e)
{
  // forward the event to the camera manipulator
  m_camera_manip_rot->mouseReleaseEvent(e);
  m_camera_manip_zoom->mouseReleaseEvent(e);
}
void Main_widget::timerEvent(QTimerEvent*)
{
  update();
}

//class GUI_event_handler
//{
//public:
//  void keyPressEvent(QKeyEvent* event);
//
//protected:
//  virtual void key_press_event(QKeyEvent* event);
//};

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


void Main_widget::verify_antarctica_node_is_redundant()
{
  Kml::Node n1(178.277211542064, -84.4725179992025),
    n2(180.0, -84.71338),
    n3(-179.942499356179, -84.7214433735525);

  // 1) check if it is collinear with its neighboring nodes:
  // all of the vectors in 3D must lie in the same plane
  auto v1 = n1.get_coords_3f();
  auto v2 = n2.get_coords_3f();
  auto v3 = n3.get_coords_3f();
  auto n = QVector3D::crossProduct(v1, v3);
  n.normalize();
  std::cout << "*** DOT PRODUCT = " << QVector3D::dotProduct(n, v2) << std::endl;

  // 2) check if it is between its neighbors (check if r,s > 0)
  auto det = [](float ax, float ay, float bx, float by) { return ax * by - ay * bx; };
  auto D = det(v1.x(), v1.y(), v3.x(), v3.y());
  auto Dr = det(v2.x(), v2.y(), v3.x(), v3.y());
  auto Ds = det(v1.x(), v1.y(), v2.x(), v2.y());
  auto r = Dr / D;
  auto s = Ds / D;
  std::cout << "r = " << r << "\ns=" << s << std::endl;
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



std::unique_ptr<Line_strips>   new_faces;

#include "Triangles.h"
std::unique_ptr<Triangles>  g_all_triangles;
std::vector<std::unique_ptr<Triangles>>  g_country_triangles;


void Main_widget::initializeGL()
{
  // verify that the node (180.0, -84.71338) in Antarctica is redundant
  //verify_antarctica_node_is_redundant();

  //init_problematic_nodes();


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
    new_faces = std::make_unique<Line_strips>(created_faces);
  }

  // SAVING ARR
  if(0)
  {
    //Aos::save_arr(m_countries, "C:/work/gsoc2023/ne_110m_admin_0_countries_africa_1.json");
    Aos::save_arr(m_countries, "C:/work/gsoc2023/ne_110m_admin_0_countries_equatorial_guinea.json");
    qDebug() << "done saving!";
    exit(0);
    //Aos::save_arr(m_countries, "C:/work/gsoc2023/ne_110m_admin_0_countries.json");
  }

  // triangulation
  {
    qDebug() << "constructiong arr..";
    //auto arrh = Aos::construct(m_countries);
    auto arrh = Aos::load_arr("C:/work/gsoc2023/ne_110m_admin_0_countries.json");
    if (arrh == nullptr)
    {
      qDebug() << "** FAILED TO LOAD THE ARRANGEMENT!!!";
      exit(1);
    }

    qDebug() << "generating triangles..";
    //auto triangle_points = Aos::get_triangles(arrh);
    auto country_triangles_map = Aos::get_triangles_by_country(arrh);
    qDebug() << "num countries = " << country_triangles_map.size();
    auto rndm = [] {return rand() / double(RAND_MAX); };
    for (auto& [country_name, triangle_points] : country_triangles_map)
    {
      auto country_triangles = std::make_unique<Triangles>(triangle_points);
      country_triangles->set_color(QVector4D(rndm(), rndm(), rndm(), 1));
      g_country_triangles.push_back(std::move(country_triangles));
    }
    //qDebug() << "num triangles = " << triangle_points.size() / 3;
    //g_all_triangles = std::make_unique<Triangles>(triangle_points);
    
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
void Main_widget::find_minimum_projected_error_on_sphere(float we)
{
  QRect vp(0, 0, m_vp_width, m_vp_height);
  auto proj = m_camera.get_projection_matrix();
  auto view = m_camera.get_view_matrix();
  QMatrix4x4 model;
  auto model_view = view * model;

  float max_err = 0;
  float max_theta = -1;
  float max_phi = -1;

  int num_divs = 200;
  const float dtheta = M_PI_2 / num_divs;
  const float dphi = M_PI_2 / num_divs;

  const float r1 = 1.f;
  const float r2 = r1 - we;
  for (int i = 0; i <= num_divs; i++)
  {
    const float theta = dtheta * i;
    const float cos_theta = std::cos(theta);
    const float sin_theta = std::sin(theta);

    for (int j = 0; j <= num_divs; j++)
    {
      QVector3D p1, p2;
      const float phi = dphi * j;
      const float cos_phi = std::cos(phi);
      const float sin_phi = std::sin(phi);

      // p1
      const float r1xz = r1 * sin_phi;
      p1.setY(r1 * cos_phi);
      p1.setX(r1xz * cos_theta);
      p1.setZ(r1xz * sin_theta);

      // p2
      const float r2xz = r2 * sin_phi;
      p2.setY(r2 * cos_phi);
      p2.setX(r2xz * cos_theta);
      p2.setZ(r2xz * sin_theta);

      auto wp1 = p1.project(model_view, proj, vp);
      auto wp2 = p2.project(model_view, proj, vp);

      const auto pe = wp1.distanceToPoint(wp2);
      if (max_err < pe)
      {
        max_err = pe;
        max_theta = theta;
        max_phi = phi;
      }
    }
  }

  std::cout << "max err = " << max_err << std::endl;
  std::cout << "max phi = " << max_phi * 180 / M_PI << std::endl;
  std::cout << "max theta = " << max_theta * 180 / M_PI << std::endl;

  auto wp1 = QVector3D(0, r1, 0).project(model_view, proj, vp);
  auto wp2 = QVector3D(0, r2, 0).project(model_view, proj, vp);
  auto pe = wp1.distanceToPoint(wp2);
  std::cout << "polar err = " << pe << std::endl;

  wp1 = QVector3D(r1, 0, 0).project(model_view, proj, vp);
  wp2 = QVector3D(r2, 0, 0).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "x-axis err = " << pe << std::endl;

  wp1 = QVector3D(0, 0, 1).project(model_view, proj, vp);
  wp2 = QVector3D(we, 0, 1).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "nearest proj err = " << pe << std::endl;

  wp1 = QVector3D(0, 0, -1).project(model_view, proj, vp);
  wp2 = QVector3D(we, 0, -1).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "farthest proj err = " << pe << std::endl;

  // project the origin on the screen (to check if it projects to the mid-vp)
  //std::cout << QVector3D(0, 0, 0).project(model_view, proj, vp) << std::endl;
}

void Main_widget::resizeGL(int w, int h)
{
  m_camera_manip_rot->resizeGL(w, h);
  
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

  QMatrix4x4 model;
  model.rotate(-90, 1,0,0); // this makes z-axes point upwards!
  const auto view = m_camera.get_view_matrix();
  const auto projection = m_camera.get_projection_matrix();
  const auto mvp = projection * view * model;

  // compute the cutting plane
// remember that we are passing the local vertex positions of the sphere 
// between the vertex and fragment shader stages, so we need to convert
// the camera-pos in world coords to sphere's local coords!
  auto c = model.inverted() * m_camera.get_pos();
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
    auto sphere_color = QVector4D(167, 205, 242, 255) / 255;
    sp.set_uniform("u_color", sphere_color);
    sp.set_uniform("u_plane", QVector4D(0,0,0,0));
    //sp.set_uniform("u_color", m_sphere->get_color());
    
    //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    m_sphere->draw();

    // DRAW SOLID FACES
    if(1)
    {
      glDisable(GL_DEPTH_TEST);
      auto face_color = QVector4D(241, 141, 0, 255) / 255;
      sp.set_uniform("u_color", face_color);
      sp.set_uniform("u_plane", plane);
      //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      //g_all_triangles->draw();
      for (auto& country : g_country_triangles)
      {
        sp.set_uniform("u_color", country->get_color());
        country->draw();
      }

      //sp.set_uniform("u_color", QVector4D(0,0,0,1));
      //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      //g_all_triangles->draw();
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
    float a = 0.0;
    sp.set_uniform("u_color", QVector4D(a, a, a, 1));
    for(auto& country_border : m_country_borders)
      country_border->draw();

    // draw the SELECTED COUNTRY in BLUE
    auto& selected_countru = m_country_borders[m_selected_country_index];
    sp.set_uniform("u_color", QVector4D(0, .6, 1, 1));
    selected_countru->draw();

    // draw the CURRENT ARC of the selected country in YELLOW
    sp.set_uniform("u_color", QVector4D(1, 1, 0, 1));
    selected_countru->draw(m_selected_arc_index);

    const QVector4D vertex_color(1, 0, 0, 1);
    sp.set_uniform("u_color", vertex_color);
    glPointSize(3);
    //m_vertices->draw();

    sp.set_uniform("u_color", QVector4D(0,1,0,1));
    glPointSize(2);
    //m_problematic_vertices->draw();
    draw_safe(m_problematic_vertices);

    // NEW FACES in RED
    sp.set_uniform("u_color", QVector4D(1, 0, 0, 1));
    //new_faces->draw();


    sp.unuse();
  }
}
