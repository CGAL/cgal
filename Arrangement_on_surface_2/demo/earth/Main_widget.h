// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#ifndef MAIN_WIDGET_H
#define MAIN_WIDGET_H

#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

#include <memory>

#include <qopenglwidget.h>

#include "Camera.h"
#include "Common_defs.h"
#include "Geodesic_arcs.h"
#include "Shader_program.h"
#include "Sphere.h"
#include "World_coordinate_axes.h"


class Main_widget : public QOpenGLWidget, protected OpenGLFunctionsBase
{
  Q_OBJECT

public:
  using QOpenGLWidget::QOpenGLWidget;
  ~Main_widget();

protected:
  void set_mouse_button_pressed_flag(QMouseEvent* e, bool flag);
  void mousePressEvent(QMouseEvent* e) override;
  void mouseMoveEvent(QMouseEvent* e) override;
  void mouseReleaseEvent(QMouseEvent* e) override;
  void timerEvent(QTimerEvent* e) override;

  void initializeGL() override;
  void resizeGL(int w, int h) override;
  void paintGL() override;

    
  void init_camera();
  void init_geometry();
  
  void init_shader_programs();
  void init_sp_smooth();
  void init_sp_per_vertex_color();
  void init_sp_arc();

private:
  // Objects in the scene
  std::unique_ptr<Sphere>           m_sphere;
  std::unique_ptr<World_coord_axes> m_world_coord_axes;
  std::unique_ptr<Geodesic_arcs>    m_geodesic_arcs;

  // Shaders
  Shader_program  m_sp_smooth;
  Shader_program  m_sp_per_vertex_color;
  Shader_program  m_sp_arc;
  
  // Camera & controls
  Camera m_camera;
  bool m_left_mouse_button_down = false;
  bool m_middle_mouse_button_down = false;
  QVector2D m_last_mouse_pos;
  float m_theta = 0, m_phi = 0;
  int   m_vp_width = 0, m_vp_height = 0;

  // Timer for continuous screen-updates
  QBasicTimer m_timer;
};

#endif
