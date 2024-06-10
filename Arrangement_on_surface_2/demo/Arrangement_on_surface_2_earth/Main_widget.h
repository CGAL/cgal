// Copyright(c) 2023, 2024  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef MAIN_WIDGET_H
#define MAIN_WIDGET_H

#include <functional>
#include <memory>

#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

#include <qopenglwidget.h>

#include "Aos.h"
#include "Camera.h"
#include "Camera_manip.h"
#include "Common_defs.h"
#include "GUI_event_handler.h"
#include "Kml_reader.h"
#include "Line_strips.h"
#include "Shader_program.h"
#include "Single_vertex.h"
#include "Sphere.h"
#include "Triangles.h"
#include "Vertices.h"
#include "World_coordinate_axes.h"

class Main_widget : public QOpenGLWidget, protected OpenGLFunctionsBase {
  Q_OBJECT

public:
  using QOpenGLWidget::QOpenGLWidget;

  Main_widget(const QString& file_name);
  ~Main_widget();

  auto& get_camera() { return m_camera; }
  auto& get_model_matrix() { return m_model; }
  auto& get_arr_handle() { return m_arrh; }

  void set_mouse_pos(const QVector3D mouse_pos);

  void hightlight_country(const std::string& country_name);

protected:
  void mousePressEvent(QMouseEvent* e) override;
  void mouseMoveEvent(QMouseEvent* e) override;
  void mouseReleaseEvent(QMouseEvent* e) override;
  // void timerEvent(QTimerEvent* e) override;
  void keyPressEvent(QKeyEvent* event) override;

  void initializeGL() override;
  void resizeGL(int w, int h) override;
  void paintGL() override;


  void init_camera();
  void init_geometry();
  void init_shader_programs();

  void init_country_borders(float error);

  // This is called when the required approximation of the arcs is below the
  // currently required one defined by the zoom level and window size. If you
  // zoom-in or increase the window-size this can be called. But once a minimum
  // approximation error is needed, it will stay there until further change.
  // SEE the definition of "m_current_approx_error" member variable below!
  float compute_backprojected_error(float pixel_error);

private:
  // COUNTRY ARRANGEMENT SPECIFIC DATA
  QString m_file_name;
  Aos::Arr_handle m_arrh;
  std::unique_ptr<Line_strips> m_gr_all_country_borders;

  // used when dimming / highlighting selected countries
  const float m_dimming_factor = 0.4f;

  // GUI: event handler for picking with right mouse button
  std::unique_ptr<GUI_event_handler> m_pick_handler;

  // These are used to highlight the picked position by right-mouse click
  std::unique_ptr<Single_vertex> m_gr_mouse_vertex;


  // TRIANGLES for rendering the countries in solid
  std::unique_ptr<Triangles> m_gr_all_triangles;
  std::map<std::string, std::unique_ptr<Triangles>> m_gr_country_triangles;


  // -------------------------------
  // --> COMMON SETUP FOR ALL SCENES

  // Basic objects in the scene
  std::unique_ptr<Sphere> m_gr_sphere;
  std::unique_ptr<World_coord_axes> m_gr_world_coord_axes;
  std::unique_ptr<Line_strips> m_gr_identification_curve;

  // Shaders
  Shader_program m_sp_smooth;
  Shader_program m_sp_per_vertex_color;
  Shader_program m_sp_arc;

  // Camera & controls
  Camera  m_camera;
  std::unique_ptr<GUI_event_handler> m_camera_manip_rot;
  std::unique_ptr<GUI_event_handler> m_camera_manip_zoom;
  QMatrix4x4 m_model;

  // view-port
  int m_vp_width = 0;
  int m_vp_height = 0;

  // After zooming in or making the viewport larger, the approximation-error
  // needs to be updated and checked against the old value. If a lower approxi-
  // mation error is needed the necessary graphics-side updates need to be made
  // INSIDE the paintGL (or wherever the OpenGL context is active)!
  bool m_update_approx_error = false;
  float m_current_approx_error;
  std::size_t m_num_uniform_points;

  // Timer for continuous screen-updates
  // QBasicTimer m_timer;

  // <-- COMMON SETUP FOR ALL SCENES
  // -------------------------------
};

#endif
