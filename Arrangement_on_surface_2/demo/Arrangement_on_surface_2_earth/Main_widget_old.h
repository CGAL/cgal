// Copyright (c) 2023, 2024  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef MAIN_WIDGET_H
#define MAIN_WIDGET_H

#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

#include <functional>
#include <memory>

#include <qopenglwidget.h>

#include "Aos.h"
#include "Camera.h"
#include "Camera_manip.h"
#include "Common_defs.h"
#include "GUI_event_handler.h"
#include "Kml_reader.h"
#include "Line_strips.h"
#include "Shader_program.h"
#include "SingleVertex.h"
#include "Sphere.h"
#include "Triangles.h"
#include "Vertices.h"
#include "World_coordinate_axes.h"

class Main_widget : public QOpenGLWidget, protected OpenGLFunctionsBase {
  Q_OBJECT

public:
  using QOpenGLWidget::QOpenGLWidget;
  ~Main_widget();

  auto& get_camera() { return m_camera; }
  auto& get_model_matrix() { return m_model; }
  auto& get_arr_handle() { return m_arrh; }

  void set_mouse_pos(const QVector3D mouse_pos) { m_mouse_pos = mouse_pos; }

  void hightlight_country(const std::string& country_name);

protected:
  void mousePressEvent(QMouseEvent* e) override;
  void mouseMoveEvent(QMouseEvent* e) override;
  void mouseReleaseEvent(QMouseEvent* e) override;
  void timerEvent(QTimerEvent* e) override;
  void keyPressEvent(QKeyEvent* event) override;

  void initializeGL() override;
  void resizeGL(int w, int h) override;
  void paintGL() override;


  void init_camera();
  void init_geometry();
  void init_shader_programs();

  void init_country_borders(float error);
  void init_country_selection();

  void handle_country_picking(QMouseEvent* e);

  // This is called when the required approximation of the arcs is below the
  // currently required one defined by the zoom level and window size. If you
  // zoom-in or increase the window-size this can be called. But once a minimum
  // approximation error is needed, it will stay there until further change.
  // SEE the definition of "m_current_approx_error" member variable below!
  float compute_backprojected_error(float pixel_error);


  // init problematic vertices: these are the vertices incident to deg-4 vertex
  void init_problematic_nodes();

private:
  // ARRANGEMENT
  Aos::Arr_handle m_arrh;
  std::unique_ptr<Line_strips> m_gr_all_approx_arcs;

  // GUI: event handler for picking with right mouse button
  std::unique_ptr<GUI_event_handler> m_pick_handler;

  // Objects in the scene
  std::unique_ptr<Sphere> m_sphere;
  std::unique_ptr<World_coord_axes> m_world_coord_axes;
  std::unique_ptr<Line_strips> m_geodesic_arcs;
  std::unique_ptr<Vertices> m_vertices, m_problematic_vertices;
  std::unique_ptr<Line_strips> m_identification_curve;

  // New faces not in the KML-file but created during arr-construction.
  // This is used to identify the Caspian Sea!
  std::unique_ptr<Line_strips>   m_new_faces;

  // These are used to highlight the picked position by right-mouse click
  QVector3D m_mouse_pos;
  std::unique_ptr<SingleVertex> m_mouse_vertex;

  // COUNTRY DATA
  Kml::Placemarks m_countries;
  std::vector<std::string> m_country_names;
  std::vector<std::unique_ptr<Line_strips>> m_country_borders;

  // boundary-arcs by country
  int m_selected_country_index;
  int m_selected_arc_index;
  Kml::Nodes m_selected_country_nodes;
  Kml::Arcs m_selected_country_arcs;
  Kml::Placemark* m_selected_country;

  // TRIANGLES for rendering the countries in solid
  std::unique_ptr<Triangles>  m_all_triangles;
  std::map<std::string, std::unique_ptr<Triangles>>  m_country_triangles;

  // Shaders
  Shader_program m_sp_smooth;
  Shader_program m_sp_per_vertex_color;
  Shader_program m_sp_arc;

  // Camera & controls
  Camera  m_camera;
  std::unique_ptr<GUI_event_handler> m_camera_manip_rot;
  std::unique_ptr<Camera_manip> m_camera_manip_zoom;
  QMatrix4x4  m_model;

  // view-port
  int m_vp_width = 0;
  int m_vp_height = 0;

  // After zooming in or making the viewport larger, the approximation-error
  // needs to be updated and checked against the old value. If a lower approxi-
  // mation error is needed the necessary graphics-side updates need to be made
  // INSIDE the paintGL (or wherever the OpenGL context is active)!
  bool m_update_approx_error = false;
  float m_current_approx_error;

  // Timer for continuous screen-updates
  QBasicTimer m_timer;
};

#endif
